"""Robustness summaries for TranSiGen-style Sci-Plex 3 adaptation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"
TRANSIGEN = SCI / "transigen_metrics_long.csv"
BASELINE = ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_summary_by_seed.csv"

SPLIT_MAP = {
    "random": "split_random",
    "cell_heldout": "split_cell_heldout",
    "scaffold_heldout": "split_scaffold_heldout",
    "cell_scaffold_heldout": "split_cell_scaffold_heldout",
    "joint_cell_scaffold_heldout": "split_cell_scaffold_heldout",
}


def iqr(x: pd.Series) -> float:
    return float(x.quantile(0.75) - x.quantile(0.25))


def ci(values: np.ndarray, n_boot: int = 10000, seed: int = 1729) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def main() -> None:
    trans = pd.read_csv(TRANSIGEN)
    trans = trans[trans["status"].eq("completed")].copy()
    trans["split_key"] = trans["split_type"].replace({"joint_cell_scaffold_heldout": "cell_scaffold_heldout"})
    trans_seed = trans[
        [
            "dataset",
            "seed",
            "split_key",
            "split_type",
            "model_name",
            "mean_all_gene_pearson",
            "mean_all_gene_rmse",
            "mean_top50_overlap",
            "mean_top100_overlap",
            "mean_top200_overlap",
            "mean_top50_direction_agreement",
            "mean_top100_direction_agreement",
            "mean_top200_direction_agreement",
            "n_train",
            "n_val",
            "n_test",
            "runtime_seconds",
        ]
    ].copy()
    trans_seed.to_csv(SCI / "transigen_seed_level_metrics.csv", index=False)

    summary = (
        trans_seed.groupby(["dataset", "model_name", "split_key"], as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_pearson=("mean_all_gene_pearson", "mean"),
            sd_pearson=("mean_all_gene_pearson", "std"),
            median_pearson=("mean_all_gene_pearson", "median"),
            iqr_pearson=("mean_all_gene_pearson", iqr),
            mean_rmse=("mean_all_gene_rmse", "mean"),
            sd_rmse=("mean_all_gene_rmse", "std"),
            median_rmse=("mean_all_gene_rmse", "median"),
            iqr_rmse=("mean_all_gene_rmse", iqr),
            mean_top50_overlap=("mean_top50_overlap", "mean"),
            mean_top100_overlap=("mean_top100_overlap", "mean"),
            mean_top200_overlap=("mean_top200_overlap", "mean"),
        )
    )
    summary.to_csv(SCI / "transigen_summary_stats.csv", index=False)

    wide = trans_seed.pivot(index="seed", columns="split_key", values="mean_all_gene_pearson")
    contrast_rows = []
    for a, b in [
        ("random", "scaffold_heldout"),
        ("random", "cell_scaffold_heldout"),
        ("scaffold_heldout", "cell_scaffold_heldout"),
    ]:
        if a in wide.columns and b in wide.columns:
            diff = wide[a] - wide[b]
            lo, hi = ci(diff.to_numpy())
            contrast_rows.append(
                {
                    "dataset": "sciplex3_24h_top2000",
                    "model_name": "transigen_adapted_sciplex3",
                    "contrast": f"{a}_minus_{b}",
                    "n_seeds": int(diff.notna().sum()),
                    "mean_difference": float(diff.mean()),
                    "sd_difference": float(diff.std(ddof=1)),
                    "median_difference": float(diff.median()),
                    "iqr_difference": iqr(diff),
                    "ci_low": lo,
                    "ci_high": hi,
                }
            )
    paired_contrasts = pd.DataFrame(contrast_rows)
    paired_contrasts.to_csv(SCI / "transigen_paired_contrasts.csv", index=False)
    paired_contrasts.to_csv(SCI / "transigen_bootstrap_ci.csv", index=False)

    base = pd.read_csv(BASELINE)
    base_keep = base[
        base["baseline"].isin(["global_train_mean", "ridge_drug_fp", "ridge_cell_dose_drug_fp", "svd50_ridge_cell_dose_drug_fp"])
    ].copy()
    base_keep["split_key"] = base_keep["split_label"].replace({"joint": "cell_scaffold_heldout"})
    comparison_rows = []
    for baseline, g in base_keep.groupby("baseline"):
        bwide = g.pivot_table(index="seed", columns="split_key", values="mean_rowwise_pearson", aggfunc="mean")
        for split in sorted(set(wide.columns).intersection(set(bwide.columns))):
            merged = pd.DataFrame({"transigen": wide[split], "baseline": bwide[split]}).dropna()
            if merged.empty:
                continue
            diff = merged["transigen"] - merged["baseline"]
            lo, hi = ci(diff.to_numpy())
            comparison_rows.append(
                {
                    "dataset": "sciplex3_24h_top2000",
                    "split_key": split,
                    "model_name": "transigen_adapted_sciplex3",
                    "baseline": baseline,
                    "n_seeds": int(len(diff)),
                    "mean_transigen": float(merged["transigen"].mean()),
                    "mean_baseline": float(merged["baseline"].mean()),
                    "mean_difference_transigen_minus_baseline": float(diff.mean()),
                    "sd_difference": float(diff.std(ddof=1)),
                    "ci_low": lo,
                    "ci_high": hi,
                    "notes": "SVD-ridge comparison appears only if available in the Sci-Plex 3 baseline table.",
                }
            )
    pd.DataFrame(comparison_rows).to_csv(SCI / "transigen_vs_baseline_paired_differences.csv", index=False)


if __name__ == "__main__":
    main()
