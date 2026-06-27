"""PRnet adapter stub.

This file defines the interface expected by the benchmark. It intentionally
does not fabricate PRnet predictions when the external implementation or
required input representation is unavailable.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from src.models.base import ModelCapabilities


class PRNetAdapter:
    name = "prnet"
    capabilities = ModelCapabilities(
        requires_raw_single_cell=True,
        supports_pseudobulk_response=False,
        supports_drug_structure=True,
        supports_dose=False,
        supports_cell_context=True,
    )

    def __init__(self) -> None:
        try:
            import prnet  # noqa: F401
        except Exception as exc:  # pragma: no cover - environment dependent
            raise ImportError("PRnet is not installed in the current environment.") from exc

    def fit(self, train_data: Any, val_data: Any | None = None, config: dict[str, Any] | None = None) -> None:
        raise NotImplementedError("External PRnet training wrapper has not been configured.")

    def predict(self, test_data: Any) -> Any:
        raise NotImplementedError("External PRnet prediction wrapper has not been configured.")

    def save(self, path: Path) -> None:
        path.write_text("PRnet adapter placeholder; no trained state available.\n", encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> "PRNetAdapter":
        return cls()
