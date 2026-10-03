"""Summarize seed-paired nonlinear matched-random controls for Sci-Plex 3."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30"
PRIMARY = ROOT / "results/deep_model_panel/sciplex3/transigen_main30_metrics_long.csv"


def bootstrap(values: np.ndarray, seed: int = 20260722, n_boot: int = 10_000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)])
    return tuple(np.quantile(estimates, [0.025, 0.975]).tolist())


def main() -> None:
    control = pd.read_csv(OUT / "metrics_long.csv")
    control = control[control.status.eq("completed")].copy()
    base = pd.read_csv(PRIMARY)
    joint = base[
        base.split_type.eq("cell_scaffold_heldout") & base.status.eq("completed")
    ][["seed", "mean_all_gene_pearson", "n_train", "n_test"]].rename(
        columns={"mean_all_gene_pearson": "joint_pearson", "n_train": "joint_n_train", "n_test": "joint_n_test"}
    )
    summary = control.groupby("split_type").agg(
        n_seeds=("seed", "nunique"),
        mean_pearson=("mean_all_gene_pearson", "mean"),
        sd_pearson=("mean_all_gene_pearson", "std"),
        median_pearson=("mean_all_gene_pearson", "median"),
        mean_n_train=("n_train", "mean"),
        mean_n_test=("n_test", "mean"),
    ).reset_index()
    summary.to_csv(OUT / "summary.csv", index=False)

    rows = []
    paired = []
    for control_name, group in control.groupby("split_type"):
        merged = group.merge(joint, on="seed", validate="one_to_one")
        merged["difference"] = merged["mean_all_gene_pearson"] - merged["joint_pearson"]
        low, high = bootstrap(merged.difference.to_numpy(), seed=20260722 + len(rows))
        rows.append({
            "model_name": "control_conditioned_nonlinear_benchmark",
            "comparison": f"{control_name} minus joint held-out",
            "n_paired_seeds": len(merged),
            "mean_difference": merged.difference.mean(),
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
            "inference_unit": "seed",
        })
        merged.insert(1, "control", control_name)
        paired.append(merged)
    pd.DataFrame(rows).to_csv(OUT / "paired_contrasts.csv", index=False)
    pd.concat(paired, ignore_index=True).to_csv(OUT / "paired_seed_level.csv", index=False)

    audit = pd.read_csv(OUT / "leakage_audit.csv")
    fields = ["same_drug_overlap", "same_scaffold_overlap", "same_cell_line_overlap", "same_numeric_dose_overlap"]
    audit.groupby("split_type")[fields].mean().reset_index().to_csv(OUT / "leakage_audit_summary.csv", index=False)


if __name__ == "__main__":
    main()
