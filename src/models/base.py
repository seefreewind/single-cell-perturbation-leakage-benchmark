"""Shared model adapter contract for leakage-aware perturbation benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol


@dataclass(frozen=True)
class ModelCapabilities:
    requires_raw_single_cell: bool
    supports_pseudobulk_response: bool
    supports_drug_structure: bool
    supports_dose: bool
    supports_cell_context: bool


@dataclass
class PredictionManifestRow:
    record_id: str
    dataset: str
    split_type: str
    seed: int
    model_name: str
    true_response_vector_path: str | None = None
    true_response_hash: str | None = None
    pred_response_vector_path: str | None = None
    pred_response_array: str | None = None
    status: str = "pending"
    runtime_seconds: float | None = None
    notes: str = ""


class PerturbationModel(Protocol):
    name: str
    capabilities: ModelCapabilities

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        """Fit the model using only training data and optional train-internal validation data."""

    def predict(self, test_data: Any) -> Any:
        """Return predictions aligned to the test_data record order."""

    def save(self, path: Path) -> None:
        """Save model state or an adapter manifest."""

    @classmethod
    def load(cls, path: Path) -> "PerturbationModel":
        """Load a model state or adapter manifest."""
