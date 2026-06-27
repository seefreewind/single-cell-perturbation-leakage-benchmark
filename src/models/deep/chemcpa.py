"""chemCPA/CPA adapter stub."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from src.models.base import ModelCapabilities


class ChemCPAAdapter:
    name = "chemcpa_or_cpa"
    capabilities = ModelCapabilities(
        requires_raw_single_cell=True,
        supports_pseudobulk_response=False,
        supports_drug_structure=True,
        supports_dose=True,
        supports_cell_context=True,
    )

    def __init__(self) -> None:
        try:
            import chemCPA  # noqa: F401
        except Exception:
            try:
                import cpa  # noqa: F401
            except Exception as exc:  # pragma: no cover - environment dependent
                raise ImportError("Neither chemCPA nor CPA is installed in the current environment.") from exc

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        raise NotImplementedError("External chemCPA/CPA training wrapper has not been configured.")

    def predict(self, test_data: Any) -> Any:
        raise NotImplementedError("External chemCPA/CPA prediction wrapper has not been configured.")

    def save(self, path: Path) -> None:
        path.write_text("chemCPA/CPA adapter placeholder; no trained state available.\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "ChemCPAAdapter":
        return cls()
