"""Frozen single-cell foundation-model embedding adapters.

The adapters require precomputed cell-state embeddings. They do not download
or run foundation models by default, which keeps benchmark split logic
separate from model-weight management.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from src.models.base import ModelCapabilities


class FrozenFoundationRidgeAdapter:
    name = "foundation_frozen_ridge"
    capabilities = ModelCapabilities(
        requires_raw_single_cell=False,
        supports_pseudobulk_response=True,
        supports_drug_structure=True,
        supports_dose=True,
        supports_cell_context=True,
    )

    def __init__(self, embedding_path: Path | None = None) -> None:
        if embedding_path is None or not embedding_path.exists():
            raise FileNotFoundError("Precomputed frozen foundation embeddings are required.")
        self.embedding_path = embedding_path

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        raise NotImplementedError("Frozen embedding ridge head is not configured for this dataset yet.")

    def predict(self, test_data: Any) -> Any:
        raise NotImplementedError("Frozen embedding ridge head is not configured for this dataset yet.")

    def save(self, path: Path) -> None:
        path.write_text(f"embedding_path={self.embedding_path}\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "FrozenFoundationRidgeAdapter":
        embedding_path = Path(path.read_text(encoding="utf-8").split("=", 1)[1].strip())
        return cls(embedding_path)


class FrozenFoundationMLPAdapter(FrozenFoundationRidgeAdapter):
    name = "foundation_frozen_mlp"

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        raise NotImplementedError("Frozen embedding MLP head is not configured for this dataset yet.")
