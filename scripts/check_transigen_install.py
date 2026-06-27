"""Check local TranSiGen repository and adaptation runtime."""

from __future__ import annotations

import importlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "external/TranSiGen"


def import_status(name: str) -> str:
    try:
        mod = importlib.import_module(name)
        return f"OK {getattr(mod, '__version__', 'unknown')}"
    except Exception as exc:
        return f"FAIL {type(exc).__name__}: {exc}"


def main() -> None:
    print(f"repo_exists={REPO.exists()}")
    if (REPO / ".git").exists():
        commit = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
        print(f"repo_commit={commit}")
    for name in ["torch", "numpy", "pandas", "sklearn", "anndata", "rdkit"]:
        print(f"{name}={import_status(name)}")


if __name__ == "__main__":
    main()
