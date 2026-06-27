"""PRnet-style Sci-Plex 3 pseudobulk adapter.

This adapter uses the official PRnet PGM module from ``external/PRnet`` when
available, but trains it against the local Sci-Plex 3 pseudobulk response
benchmark. It is a smoke-test bridge, not a full PRnet reproduction.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


ROOT = Path(__file__).resolve().parents[3]
PRNET_ROOT = ROOT / "external/PRnet"
if str(PRNET_ROOT) not in sys.path:
    sys.path.insert(0, str(PRNET_ROOT))

try:
    from models.PRnet import PGM  # type: ignore
except Exception as exc:  # pragma: no cover - exercised by smoke script.
    PGM = None
    PRNET_IMPORT_ERROR = exc
else:
    PRNET_IMPORT_ERROR = None


class PRnetAdapter:
    name = "prnet_adapted_sciplex3"
    requires_raw_single_cell = False
    supports_pseudobulk_response = True
    supports_drug_structure = True
    supports_dose = True
    supports_cell_context = False

    def __init__(
        self,
        x_dimension: int,
        drug_dimension: int = 1024,
        hidden_layer_sizes: list[int] | None = None,
        z_dimension: int = 64,
        adaptor_layer_sizes: list[int] | None = None,
        comb_dimension: int = 64,
        dr_rate: float = 0.05,
        seed: int = 0,
    ):
        if PGM is None:
            raise RuntimeError(f"Could not import official PRnet PGM from {PRNET_ROOT}: {PRNET_IMPORT_ERROR}")
        torch.manual_seed(seed)
        self.config = {
            "x_dimension": x_dimension,
            "drug_dimension": drug_dimension,
            "hidden_layer_sizes": hidden_layer_sizes or [128],
            "z_dimension": z_dimension,
            "adaptor_layer_sizes": adaptor_layer_sizes or [128],
            "comb_dimension": comb_dimension,
            "dr_rate": dr_rate,
            "seed": seed,
        }
        self.model = PGM(
            x_dim=x_dimension,
            c_dim=comb_dimension,
            n_dim=10,
            hidden_layer_sizes=self.config["hidden_layer_sizes"],
            z_dimension=z_dimension,
            adaptor_layer_sizes=self.config["adaptor_layer_sizes"],
            comb_adapt_dim=drug_dimension,
            dr_rate=dr_rate,
        )
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.model.to(self.device)

    def _loader(self, control: np.ndarray, treated: np.ndarray, drug_dose: np.ndarray, batch_size: int, shuffle: bool, seed: int):
        n = control.shape[0]
        rng = np.random.default_rng(seed)
        order = np.arange(n)
        if shuffle:
            order = rng.permutation(order)
        for start in range(0, n, batch_size):
            idx = order[start : start + batch_size]
            yield (
                torch.from_numpy(control[idx]).to(self.device, dtype=torch.float32),
                torch.from_numpy(treated[idx]).to(self.device, dtype=torch.float32),
                torch.from_numpy(drug_dose[idx]).to(self.device, dtype=torch.float32),
            )

    def fit(
        self,
        train_data: dict[str, np.ndarray],
        val_data: dict[str, np.ndarray],
        config: dict[str, Any] | None = None,
    ) -> None:
        cfg = {
            "epochs": 2,
            "batch_size": 128,
            "learning_rate": 1e-3,
            "weight_decay": 1e-8,
            "patience": 2,
        }
        if config:
            cfg.update(config)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=cfg["learning_rate"], weight_decay=cfg["weight_decay"])
        loss_fn = nn.GaussianNLLLoss()
        best_state = None
        best_loss = float("inf")
        bad_epochs = 0
        seed = int(self.config["seed"])

        for epoch in range(int(cfg["epochs"])):
            self.model.train()
            for control, treated, drug_dose in self._loader(
                train_data["control"], train_data["treated"], train_data["drug_dose"], int(cfg["batch_size"]), True, seed + epoch
            ):
                noise = torch.randn(control.size(0), 10, device=self.device)
                recon = self.model(control, drug_dose, noise)
                dim = recon.size(1) // 2
                mean = recon[:, :dim]
                var = F.softplus(recon[:, dim:]) + 1e-4
                loss = loss_fn(mean, treated, var)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            val_loss = self._validation_loss(val_data, int(cfg["batch_size"]))
            if val_loss < best_loss:
                best_loss = val_loss
                best_state = {k: v.detach().cpu().clone() for k, v in self.model.state_dict().items()}
                bad_epochs = 0
            else:
                bad_epochs += 1
                if bad_epochs >= int(cfg["patience"]):
                    break
        if best_state is not None:
            self.model.load_state_dict(best_state)

    def _validation_loss(self, val_data: dict[str, np.ndarray], batch_size: int) -> float:
        if len(val_data["control"]) == 0:
            return float("inf")
        loss_fn = nn.GaussianNLLLoss(reduction="mean")
        losses = []
        self.model.eval()
        with torch.no_grad():
            for control, treated, drug_dose in self._loader(val_data["control"], val_data["treated"], val_data["drug_dose"], batch_size, False, int(self.config["seed"])):
                noise = torch.randn(control.size(0), 10, device=self.device)
                recon = self.model(control, drug_dose, noise)
                dim = recon.size(1) // 2
                mean = recon[:, :dim]
                var = F.softplus(recon[:, dim:]) + 1e-4
                losses.append(float(loss_fn(mean, treated, var).item()))
        return float(np.mean(losses)) if losses else float("inf")

    def predict_treated(self, test_data: dict[str, np.ndarray], batch_size: int = 128) -> np.ndarray:
        preds = []
        self.model.eval()
        with torch.no_grad():
            for control, _treated, drug_dose in self._loader(test_data["control"], test_data["treated"], test_data["drug_dose"], batch_size, False, int(self.config["seed"])):
                noise = torch.randn(control.size(0), 10, device=self.device)
                recon = self.model(control, drug_dose, noise)
                dim = recon.size(1) // 2
                preds.append(recon[:, :dim].detach().cpu().numpy())
        return np.vstack(preds).astype(np.float32)
