"""TranSiGen-style Sci-Plex 3 pseudobulk adaptation.

This is a minimal, auditable adaptation for the local benchmark: basal/control
expression, compound fingerprint, log dose, and cell-line one-hot features are
used to predict treated-control response vectors. It is not a full reproduction
of the original TranSiGen training recipe.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn


class _ResponseMLP(nn.Module):
    def __init__(self, input_dim: int, output_dim: int, hidden_dim: int = 256, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class TranSiGenAdapter:
    name = "transigen_adapted_sciplex3"
    requires_raw_single_cell = False
    supports_pseudobulk_response = True
    supports_drug_structure = True
    supports_dose = True
    supports_cell_context = True

    def __init__(self, input_dim: int, output_dim: int, hidden_dim: int = 256, dropout: float = 0.1, seed: int = 0):
        torch.manual_seed(seed)
        self.config = {
            "input_dim": input_dim,
            "output_dim": output_dim,
            "hidden_dim": hidden_dim,
            "dropout": dropout,
            "seed": seed,
        }
        self.model = _ResponseMLP(input_dim, output_dim, hidden_dim=hidden_dim, dropout=dropout)
        self.x_mean: np.ndarray | None = None
        self.x_std: np.ndarray | None = None
        self.y_mean: np.ndarray | None = None
        self.y_std: np.ndarray | None = None

    def _standardize_x(self, x: np.ndarray, fit: bool = False) -> np.ndarray:
        if fit:
            self.x_mean = x.mean(axis=0, keepdims=True)
            self.x_std = x.std(axis=0, keepdims=True)
            self.x_std[self.x_std < 1e-6] = 1.0
        if self.x_mean is None or self.x_std is None:
            raise RuntimeError("Feature scaler is not fitted.")
        return ((x - self.x_mean) / self.x_std).astype(np.float32)

    def _standardize_y(self, y: np.ndarray, fit: bool = False) -> np.ndarray:
        if fit:
            self.y_mean = y.mean(axis=0, keepdims=True)
            self.y_std = y.std(axis=0, keepdims=True)
            self.y_std[self.y_std < 1e-6] = 1.0
        if self.y_mean is None or self.y_std is None:
            raise RuntimeError("Target scaler is not fitted.")
        return ((y - self.y_mean) / self.y_std).astype(np.float32)

    def fit(self, train_data: dict[str, np.ndarray], val_data: dict[str, np.ndarray] | None = None, config: dict[str, Any] | None = None) -> None:
        cfg = {
            "epochs": 20,
            "batch_size": 64,
            "learning_rate": 1e-3,
            "weight_decay": 1e-5,
            "patience": 5,
        }
        if config:
            cfg.update(config)
        x_train = self._standardize_x(train_data["x"], fit=True)
        y_train = self._standardize_y(train_data["y"], fit=True)
        if val_data is not None and len(val_data["x"]):
            x_val = self._standardize_x(val_data["x"], fit=False)
            y_val = self._standardize_y(val_data["y"], fit=False)
        else:
            x_val = y_val = None

        optimizer = torch.optim.AdamW(self.model.parameters(), lr=cfg["learning_rate"], weight_decay=cfg["weight_decay"])
        loss_fn = nn.MSELoss()
        n = x_train.shape[0]
        batch_size = int(cfg["batch_size"])
        best_state = None
        best_loss = float("inf")
        bad_epochs = 0
        rng = np.random.default_rng(int(self.config["seed"]))

        for _epoch in range(int(cfg["epochs"])):
            order = rng.permutation(n)
            self.model.train()
            for start in range(0, n, batch_size):
                idx = order[start : start + batch_size]
                xb = torch.from_numpy(x_train[idx])
                yb = torch.from_numpy(y_train[idx])
                optimizer.zero_grad()
                loss = loss_fn(self.model(xb), yb)
                loss.backward()
                optimizer.step()

            if x_val is not None:
                self.model.eval()
                with torch.no_grad():
                    val_loss = float(loss_fn(self.model(torch.from_numpy(x_val)), torch.from_numpy(y_val)).item())
                if val_loss < best_loss:
                    best_loss = val_loss
                    best_state = {k: v.detach().clone() for k, v in self.model.state_dict().items()}
                    bad_epochs = 0
                else:
                    bad_epochs += 1
                    if bad_epochs >= int(cfg["patience"]):
                        break
        if best_state is not None:
            self.model.load_state_dict(best_state)

    def predict(self, test_data: dict[str, np.ndarray]) -> np.ndarray:
        x = self._standardize_x(test_data["x"], fit=False)
        self.model.eval()
        with torch.no_grad():
            pred_scaled = self.model(torch.from_numpy(x)).numpy()
        if self.y_mean is None or self.y_std is None:
            raise RuntimeError("Target scaler is not fitted.")
        return (pred_scaled * self.y_std + self.y_mean).astype(np.float32)

    def save(self, path: Path) -> None:
        path.mkdir(parents=True, exist_ok=True)
        torch.save(self.model.state_dict(), path / "model.pt")
        np.savez_compressed(path / "scalers.npz", x_mean=self.x_mean, x_std=self.x_std, y_mean=self.y_mean, y_std=self.y_std)
        (path / "config.json").write_text(json.dumps(self.config, indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "TranSiGenAdapter":
        cfg = json.loads((path / "config.json").read_text(encoding="utf-8"))
        obj = cls(**cfg)
        obj.model.load_state_dict(torch.load(path / "model.pt", map_location="cpu"))
        scalers = np.load(path / "scalers.npz")
        obj.x_mean = scalers["x_mean"]
        obj.x_std = scalers["x_std"]
        obj.y_mean = scalers["y_mean"]
        obj.y_std = scalers["y_std"]
        return obj
