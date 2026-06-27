"""Smoke-run the current deep-model-panel assembly workflow."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.reporting.assemble_deep_model_panel import main


if __name__ == "__main__":
    main()
