#!/usr/bin/env python3
"""Create a BMC-specific, provenance-tracked master result source.

This script consumes only completed fixed-seed result tables.  It does not
train models, alter manifests, or manually enter numerical results.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TCBB = ROOT / "results/tcbb_master"
OUT = ROOT / "results/bmc_master"
N_BOOT = 10_000


def qstats(series: pd.Series) -> dict[str, float]:
    x = pd.to_numeric(series, errors="coerce").dropna()
    return {
        "mean": float(x.mean()),
        "sd": float(x.std(ddof=1)),
        "median": float(x.median()),
        "iqr": float(x.quantile(0.75) - x.quantile(0.25)),
        "n_seeds": int(x.nunique() if x.name == "seed" else len(x)),
    }


def performance_summary(performance: pd.DataFrame) -> pd.DataFrame:
    fields = ["dataset", "split", "split_label", "model", "metric", "model_status", "source_file", "evidence_role"]
    rows = []
    for key, group in performance.groupby(fields, dropna=False):
        values = group.value
        rows.append({
            **dict(zip(fields, key)),
            "mean": values.mean(), "sd": values.std(ddof=1), "median": values.median(),
            "iqr": values.quantile(.75) - values.quantile(.25), "n_seeds": group.seed.nunique(),
            "mean_n_train": group.n_train.mean(), "mean_n_test": group.n_test.mean(),
        })
    return pd.DataFrame(rows).sort_values(["dataset", "model", "metric", "split"]).reset_index(drop=True)


def bootstrap(values: pd.Series, seed: int) -> tuple[float, float]:
    x = values.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    samples = rng.choice(x, size=(N_BOOT, len(x)), replace=True).mean(axis=1)
    return float(np.quantile(samples, .025)), float(np.quantile(samples, .975))


def controls() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    op_path = ROOT / "results/openproblems_train_size_matched_100/train_size_matched_random_metrics_by_seed.csv"
    op = pd.read_csv(op_path)
    op = op.rename(columns={"control": "split", "value": "value"})
    op["dataset"] = "OpenProblems"
    op["split_label"] = op["split"].str.replace("_", " ", regex=False)
    op["model"] = "ridge_cell_drug_fp"
    op["metric"] = "rowwise_pearson"
    op["model_status"] = "completed_baseline"
    op["source_file"] = str(op_path.relative_to(ROOT))
    op["evidence_role"] = "primary_probe_control"
    op = op[["seed", "split", "model", "n_train", "n_test", "value", "dataset", "metric", "model_status", "source_file", "split_label", "evidence_role"]]

    sci_path = ROOT / "results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30/metrics_long.csv"
    sci = pd.read_csv(sci_path)
    sci = sci[sci.status.eq("completed")].copy()
    sci["dataset"] = "Sci-Plex 3"
    sci["split"] = sci["split_type"]
    sci["split_label"] = sci["split"].str.replace("_", " ", regex=False)
    sci["model"] = "control_conditioned_nonlinear_benchmark"
    sci["metric"] = "rowwise_pearson"
    sci["value"] = sci["mean_all_gene_pearson"]
    sci["model_status"] = "benchmark_compatible_adaptation"
    sci["source_file"] = str(sci_path.relative_to(ROOT))
    sci["evidence_role"] = "supporting_probe_control"
    sci = sci[["seed", "split", "model", "n_train", "n_test", "value", "dataset", "metric", "model_status", "source_file", "split_label", "evidence_role"]]

    audit_rows = []
    op_audit = pd.read_csv(ROOT / "results/openproblems_train_size_matched_100/train_size_matched_random_leakage_audit.csv")
    op_audit["dataset"] = "OpenProblems"
    op_audit["model"] = "ridge_cell_drug_fp"
    audit_rows.append(op_audit)
    sci_audit = pd.read_csv(ROOT / "results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30/leakage_audit.csv")
    sci_audit["dataset"] = "Sci-Plex 3"
    sci_audit["model"] = "control_conditioned_nonlinear_benchmark"
    audit_rows.append(sci_audit)
    return pd.concat([op, sci], ignore_index=True), pd.concat(audit_rows, ignore_index=True), pd.DataFrame()


def control_contrasts(performance: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    configurations = [
        ("OpenProblems", "ridge_cell_drug_fp", "split_cell_scaffold_heldout", ["test_size_matched_random", "train_test_size_matched_random", "composition_train_test_size_matched_random"]),
        ("Sci-Plex 3", "control_conditioned_nonlinear_benchmark", "split_cell_scaffold_heldout", ["test_size_matched_random", "train_test_size_matched_random"]),
    ]
    summary, details = [], []
    for idx, (dataset, model, strict, control_names) in enumerate(configurations):
        subset = performance[(performance.dataset.eq(dataset)) & (performance.model.eq(model)) & (performance.metric.eq("rowwise_pearson"))]
        strict_rows = subset[subset.split.eq(strict)][["seed", "value"]].rename(columns={"value": "strict_value"})
        for control in control_names:
            left = subset[subset.split.eq(control)][["seed", "value", "n_train", "n_test"]].rename(columns={"value": "control_value"})
            paired = left.merge(strict_rows, on="seed", validate="one_to_one")
            paired["difference"] = paired.control_value - paired.strict_value
            low, high = bootstrap(paired.difference, 20260722 + idx + len(summary))
            summary.append({
                "dataset": dataset, "model": model, "metric": "rowwise_pearson",
                "comparison": f"{control} minus joint held-out", "n_paired_seeds": len(paired),
                "mean_difference": paired.difference.mean(), "bootstrap_95ci_low": low,
                "bootstrap_95ci_high": high, "inference_unit": "seed",
                "source": "same_seed_fixed_random_manifest_control",
            })
            paired["dataset"] = dataset
            paired["model"] = model
            paired["control"] = control
            details.append(paired)
    return pd.concat(details, ignore_index=True), pd.DataFrame(summary)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Copy read-only legacy masters whose content remains an input to the BMC revision.
    for path in TCBB.glob("*.csv"):
        if path.name not in {"master_performance_seed_level.csv", "master_performance_summary.csv", "master_contrast_summary.csv", "master_contrast_seed_level.csv"}:
            shutil.copy2(path, OUT / path.name)

    performance = pd.read_csv(TCBB / "master_performance_seed_level.csv")
    # Use one public-facing name consistently; saved source paths retain provenance.
    performance.loc[performance.model.eq("transigen_adapted_sciplex3"), "model"] = "control_conditioned_nonlinear_benchmark"
    performance.loc[performance.model.eq("control_conditioned_nonlinear_benchmark"), "model_status"] = "benchmark_compatible_adaptation"
    performance.loc[performance.model.eq("control_conditioned_nonlinear_benchmark"), "evidence_role"] = "supporting_probe"
    added, audit, _ = controls()
    performance = pd.concat([performance, added], ignore_index=True)
    performance.to_csv(OUT / "master_performance_seed_level.csv", index=False)
    performance_summary(performance).to_csv(OUT / "master_performance_summary.csv", index=False)

    old_contrast = pd.read_csv(TCBB / "master_contrast_summary.csv")
    # Legacy TCBB matched controls used a weaker design (test-size or selected
    # test composition only). The BMC master retains only the replacement
    # train-size-aware controls built above.
    old_contrast = old_contrast[~old_contrast["comparison"].str.contains("matched_random", na=False)].copy()
    old_contrast.loc[old_contrast.model.eq("transigen_adapted_sciplex3"), "model"] = "control_conditioned_nonlinear_benchmark"
    detail, control_summary = control_contrasts(performance)
    pd.concat([old_contrast, control_summary], ignore_index=True).to_csv(OUT / "master_contrast_summary.csv", index=False)
    old_detail = pd.read_csv(TCBB / "master_contrast_seed_level.csv")
    if "comparison" in old_detail.columns:
        old_detail = old_detail[~old_detail["comparison"].str.contains("matched_random", na=False)].copy()
    old_detail.loc[old_detail.model.eq("transigen_adapted_sciplex3"), "model"] = "control_conditioned_nonlinear_benchmark"
    pd.concat([old_detail, detail], ignore_index=True).to_csv(OUT / "master_contrast_seed_level.csv", index=False)
    audit.to_csv(OUT / "master_train_size_matched_control_audit_seed_level.csv", index=False)
    # These aliases are the BMC canonical matched-control tables. They replace
    # the older test-only control names copied from the TCBB provenance folder.
    audit.to_csv(OUT / "master_matched_random_audit.csv", index=False)
    detail.to_csv(OUT / "master_matched_random_seed_level.csv", index=False)

    provenance = [
        "# BMC Master Result Provenance",
        "",
        "This directory is the single numerical source for the BMC Bioinformatics revision. It was generated from immutable fixed-manifest, completed seed-level outputs only; no model training occurs in this build step.",
        "",
        "## Added controls",
        "",
        "- OpenProblems: 100 fixed seeds for test-size-matched, train-and-test-size-matched, and cell-context plus DE-gene-count-quartile-matched random controls of the cell-plus-drug-fingerprint ridge.",
        "- Sci-Plex 3: 30 fixed seeds for test-size-matched and train-and-test-size-matched random controls of the control-conditioned nonlinear benchmark. Training uses only sampled random training records; validation is a deterministic partition of that training set.",
        "- Joint manifests supply target counts only for these controls. Their test, training, and excluded record IDs are never reused by a random-control split.",
        "",
        "## Primary source files",
        "",
        "- `results/tcbb_master/`: pre-existing fixed-manifest baseline and audit source tables.",
        "- `results/openproblems_train_size_matched_100/`: new OpenProblems controls and reconstructed audit fields.",
        "- `results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30/`: new Sci-Plex nonlinear controls, prediction manifests, audit fields, and training-loss traces.",
        "",
        "All contrasts are paired by the same manifest seed and use seed-level bootstrap resampling (10,000 resamples).",
    ]
    (OUT / "BMC_result_provenance.md").write_text("\n".join(provenance) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
