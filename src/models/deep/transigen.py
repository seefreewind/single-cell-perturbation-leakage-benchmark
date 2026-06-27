"""TranSiGen adapter stub for chemical-induced transcriptional response prediction."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from src.models.base import ModelCapabilities


class TranSiGenAdapter:
    name = "transigen"
    capabilities = ModelCapabilities(
        requires_raw_single_cell=False,
        supports_pseudobulk_response=True,
        supports_drug_structure=True,
        supports_dose=True,
        supports_cell_context=True,
    )

    def __init__(self) -> None:
        try:
            import transigen  # noqa: F401
        except Exception as exc:  # pragma: no cover - environment dependent
            raise ImportError("TranSiGen is not installed in the current environment.") from exc

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        raise NotImplementedError("External TranSiGen training wrapper has not been configured.")

    def predict(self, test_data: Any) -> Any:
        raise NotImplementedError("External TranSiGen prediction wrapper has not been configured.")

    def save(self, path: Path) -> None:
        path.write_text("TranSiGen adapter placeholder; no trained state available.\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "TranSiGenAdapter":
        return cls()
