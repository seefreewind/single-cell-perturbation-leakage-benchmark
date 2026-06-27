"""Quality-control checks for the leakage-aware deep model panel."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.metrics.vector_metrics import rmse, rowwise_pearson, topk_direction_agreement, topk_overlap  # noqa: E402


def assert_close(name: str, value: float, expected: float, tol: float = 1e-8) -> None:
    if not np.isclose(value, expected, atol=tol, equal_nan=True):
        raise AssertionError(f"{name}: expected {expected}, observed {value}")


def metric_tests() -> list[str]:
    y = np.array([[1, 2, 3], [1, 1, 1], [1, -1, 0]], dtype=float)
    pred = np.array([[1, 2, 3], [2, 2, 2], [-1, 1, 0]], dtype=float)
    corr = rowwise_pearson(y, pred)
    err = rmse(y, pred)
    assert_close("rowwise pearson identical row", corr[0], 1.0)
    if not np.isnan(corr[1]):
        raise AssertionError("constant row should produce nan Pearson")
    assert_close("rowwise pearson inverse row", corr[2], -1.0)
    assert_close("rmse identical row", err[0], 0.0)
    overlap = topk_overlap(y, pred, 2)
    direction = topk_direction_agreement(y, pred, 2)
    if overlap.shape[0] != 3 or direction.shape[0] != 3:
        raise AssertionError("top-k metrics returned wrong number of rows")
    return ["metric reproducibility toy tests passed"]


def split_integrity_tests() -> list[str]:
    messages = []
    op = pd.read_csv(ROOT / "results/deep_model_panel/all_model_leakage_audit_long.csv")
    scaffold = op[(op["dataset"] == "openproblems_neurips2023") & (op["split_type"] == "scaffold_heldout")]
    joint = op[(op["dataset"] == "openproblems_neurips2023") & (op["split_type"] == "joint_cell_scaffold_heldout")]
    moa = op[(op["dataset"] == "openproblems_neurips2023") & (op["split_type"] == "moa_heldout")]
    sci_scaffold = op[(op["dataset"] == "sciplex3_24h_top2000") & (op["split_type"] == "scaffold_heldout")]
    sci_joint = op[(op["dataset"] == "sciplex3_24h_top2000") & (op["split_type"] == "joint_cell_scaffold_heldout")]

    if not scaffold.empty:
        if scaffold["test_drug_same_in_train_frac"].fillna(scaffold.get("same_drug_in_train", 0)).max() > 0:
            raise AssertionError("OpenProblems scaffold-held-out contains same-drug leakage")
        if scaffold["test_scaffold_same_in_train_frac"].fillna(scaffold.get("same_scaffold_in_train", 0)).max() > 0:
            raise AssertionError("OpenProblems scaffold-held-out contains same-scaffold leakage")
    if not joint.empty:
        for col in ["test_drug_same_in_train_frac", "test_scaffold_same_in_train_frac", "test_cell_same_in_train_frac"]:
            if col in joint.columns and joint[col].fillna(0).max() > 0:
                raise AssertionError(f"OpenProblems joint-held-out has nonzero {col}")
    if not moa.empty and "same_moa_in_train" in moa.columns:
        if moa["same_moa_in_train"].fillna(0).max() > 0:
            raise AssertionError("OpenProblems MoA-held-out has same-MoA overlap")
    for df, name in [(sci_scaffold, "Sci-Plex scaffold"), (sci_joint, "Sci-Plex joint")]:
        if not df.empty and "same_drug_cell_in_train" in df.columns:
            if df["same_drug_cell_in_train"].fillna(0).max() > 0:
                raise AssertionError(f"{name} has same drug-cell overlap")
    messages.append("split integrity checks passed for available leakage columns")
    return messages


def output_completeness_tests() -> list[str]:
    required = [
        "all_model_metrics_long.csv",
        "all_model_leakage_audit_long.csv",
        "model_rank_stability.csv",
        "random_to_strict_contrasts.csv",
        "bootstrap_ci.csv",
    ]
    for filename in required:
        path = ROOT / "results/deep_model_panel" / filename
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"missing or empty output: {path}")
    metrics = pd.read_csv(ROOT / "results/deep_model_panel/all_model_metrics_long.csv")
    planned = metrics[metrics["status"].str.contains("not_run", na=False)]
    if planned.empty:
        raise AssertionError("planned model status rows are missing")
    return ["output completeness checks passed"]


def no_test_tuning_check() -> list[str]:
    config_text = (ROOT / "configs/deep_model_panel.yaml").read_text(encoding="utf-8")
    if "test_split_used_for_tuning: false" not in config_text:
        raise AssertionError("no-test-tuning flag is not explicit in config")
    if "validation_source: train_split_internal_only" not in config_text:
        raise AssertionError("validation source is not train-internal")
    return ["no-test-tuning config check passed"]


def main() -> None:
    messages = []
    messages.extend(metric_tests())
    messages.extend(split_integrity_tests())
    messages.extend(output_completeness_tests())
    messages.extend(no_test_tuning_check())
    out = ROOT / "results/deep_model_panel/logs/qc_deep_model_panel.log"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(messages) + "\n", encoding="utf-8")
    print("\n".join(messages))


if __name__ == "__main__":
    main()
