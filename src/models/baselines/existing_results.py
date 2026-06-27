"""Adapter marker for models whose results are produced by existing scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from src.models.base import ModelCapabilities


class ExistingResultsAdapter:
    name = "existing_results"
    capabilities = ModelCapabilities(
        requires_raw_single_cell=False,
        supports_pseudobulk_response=True,
        supports_drug_structure=True,
        supports_dose=True,
        supports_cell_context=True,
    )

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        return None

    def predict(self, test_data: Any) -> Any:
        raise NotImplementedError("Existing results are loaded from CSV summaries, not recomputed here.")

    def save(self, path: Path) -> None:
        path.write_text("Existing-results adapter; no model state.\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "ExistingResultsAdapter":
        return cls()
