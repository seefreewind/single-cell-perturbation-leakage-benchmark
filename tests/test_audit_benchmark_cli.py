from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd


def test_audit_benchmark_cli_on_example(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    records = root / "metadata/examples/splits/candidate_splits.csv"
    outdir = tmp_path / "audit"
    cmd = [
        "python3",
        "-m",
        "audit_benchmark",
        "--records",
        str(records),
        "--split-manifest",
        str(records),
        "--split-column",
        "split_cell_scaffold_heldout",
        "--assignment-column",
        "split_cell_scaffold_heldout",
        "--output-dir",
        str(outdir),
    ]
    subprocess.run(cmd, cwd=root, check=True)

    expected = [
        "record_level_audit.csv",
        "split_level_summary.csv",
        "similarity_matrix_metadata.json",
        "integrity_report.json",
        "matched_random_manifest.csv",
        "audit_report.html",
    ]
    for name in expected:
        assert (outdir / name).exists(), name

    record = pd.read_csv(outdir / "record_level_audit.csv")
    summary = pd.read_csv(outdir / "split_level_summary.csv")
    integrity = json.loads((outdir / "integrity_report.json").read_text())
    metadata = json.loads((outdir / "similarity_matrix_metadata.json").read_text())

    assert len(record) == 1
    assert int(summary.loc[0, "n_test"]) == 1
    assert int(summary.loc[0, "n_excluded"]) == 3
    assert integrity["passed"] is True
    assert metadata["rdkit_available"] is True
