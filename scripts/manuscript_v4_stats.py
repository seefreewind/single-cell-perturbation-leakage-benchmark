#!/usr/bin/env python3
"""Generate v4 manuscript statistical add-ons."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(".")
OUT = ROOT / "results/manuscript_v4_stats"
OUT.mkdir(parents=True, exist_ok=True)


def paired_bootstrap(values: pd.DataFrame, left: str, right: str, metric: str, n_boot: int = 10000, seed: int = 123) -> dict:
    wide = values.pivot(index="seed", columns="split", values=metric).dropna(subset=[left, right])
    diffs = wide[left].to_numpy() - wide[right].to_numpy()
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        boot[i] = rng.choice(diffs, size=len(diffs), replace=True).mean()
    return {
        "metric": metric,
        "contrast": f"{left} minus {right}",
        "n_seeds": len(diffs),
        "mean_difference": float(diffs.mean()),
        "ci_low": float(np.quantile(boot, 0.025)),
        "ci_high": float(np.quantile(boot, 0.975)),
    }


def openproblems_performance_ci() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "results/openproblems_multiseed_100_de/all_baseline_summary.csv")
    df = df[df["baseline"].eq("ridge_cell_drug_fp")].copy()
    rows = []
    for metric in ["mean_rowwise_pearson", "mean_de_rowwise_pearson"]:
        for right in ["split_scaffold_heldout", "split_cell_scaffold_heldout"]:
            rows.append(paired_bootstrap(df[["seed", "split", metric]], "split_random", right, metric))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "openproblems_ridge_performance_delta_ci.csv", index=False)
    return out


def sciplex_performance_ci() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_summary_by_seed.csv")
    df = df[df["baseline"].eq("ridge_cell_dose_drug_fp")].copy()
    rows = []
    for right in ["split_scaffold_heldout", "split_cell_scaffold_heldout"]:
        rows.append(paired_bootstrap(df[["seed", "split", "mean_rowwise_pearson"]], "split_random", right, "mean_rowwise_pearson"))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "sciplex3_ridge_performance_delta_ci.csv", index=False)
    return out


def rank_correlation_stats() -> pd.DataFrame:
    sci = pd.read_csv(ROOT / "results/sciplex3_24h_top2000_multiseed_30/split_rank_correlations.csv")
    op = pd.read_csv(ROOT / "results/ranking_instability_100/split_rank_correlations.csv")
    rows = []
    for name, df in [("OpenProblems", op), ("Sci-Plex 3", sci)]:
        for split, sub in df.groupby("comparison_split"):
            vals = sub["spearman_rank_correlation"].dropna().to_numpy()
            if len(vals) == 0:
                continue
            rng = np.random.default_rng(456)
            boot = np.array([rng.choice(vals, size=len(vals), replace=True).mean() for _ in range(10000)])
            rows.append(
                {
                    "dataset": name,
                    "comparison_split": split,
                    "n_seeds": len(vals),
                    "mean_spearman": float(vals.mean()),
                    "sd_spearman": float(vals.std(ddof=1)),
                    "ci_low": float(np.quantile(boot, 0.025)),
                    "ci_high": float(np.quantile(boot, 0.975)),
                    "negative_seed_fraction": float((vals < 0).mean()),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "rank_correlation_ci_and_negative_fraction.csv", index=False)
    return out


def tanimoto_performance_correlations() -> pd.DataFrame:
    op = pd.read_csv(ROOT / "results/openproblems_multiseed_100_de/all_drug_similarity_summary.csv")
    op = op[op["baseline"].eq("ridge_cell_drug_fp")].copy()
    op["dataset"] = "OpenProblems"
    sci_sim = pd.read_csv(ROOT / "results/sciplex3_24h_top2000_multiseed_30/similarity_per_row_all.csv")
    sci_perf = pd.read_csv(ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_per_row_all.csv")
    sci = sci_perf[sci_perf["baseline"].eq("ridge_cell_dose_drug_fp")].merge(
        sci_sim[["seed", "split", "sample_id", "max_train_tanimoto"]],
        on=["seed", "split", "sample_id"],
        how="left",
    )
    sci_rows = []
    for (seed, split), sub in sci.groupby(["seed", "split"]):
        data = sub[["max_train_tanimoto", "rowwise_pearson", "rmse"]].dropna()
        def corr(a: str, b: str) -> float:
            if data.shape[0] < 3 or data[a].nunique() < 2 or data[b].nunique() < 2:
                return np.nan
            return float(data[a].corr(data[b]))
        sci_rows.append(
            {
                "dataset": "Sci-Plex 3",
                "seed": seed,
                "split": split,
                "baseline": "ridge_cell_dose_drug_fp",
                "n": len(data),
                "mean_max_train_tanimoto": data["max_train_tanimoto"].mean(),
                "corr_tanimoto_pearson": corr("max_train_tanimoto", "rowwise_pearson"),
                "corr_tanimoto_rmse": corr("max_train_tanimoto", "rmse"),
            }
        )
    sci_df = pd.DataFrame(sci_rows)
    op = op.rename(columns={"corr_tanimoto_pearson": "corr_tanimoto_pearson", "corr_tanimoto_rmse": "corr_tanimoto_rmse"})
    common = ["dataset", "seed", "split", "baseline", "n", "mean_max_train_tanimoto", "corr_tanimoto_pearson", "corr_tanimoto_rmse"]
    all_rows = pd.concat([op[common], sci_df[common]], ignore_index=True)
    summary = (
        all_rows.groupby(["dataset", "split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            estimable_pearson_seeds=("corr_tanimoto_pearson", lambda s: int(s.notna().sum())),
            mean_corr_tanimoto_pearson=("corr_tanimoto_pearson", "mean"),
            sd_corr_tanimoto_pearson=("corr_tanimoto_pearson", "std"),
            mean_corr_tanimoto_rmse=("corr_tanimoto_rmse", "mean"),
            sd_corr_tanimoto_rmse=("corr_tanimoto_rmse", "std"),
        )
        .reset_index()
    )
    all_rows.to_csv(OUT / "tanimoto_performance_correlations_by_seed.csv", index=False)
    summary.to_csv(OUT / "tanimoto_performance_correlation_summary.csv", index=False)
    return summary


def main() -> None:
    print(openproblems_performance_ci().to_string(index=False))
    print(sciplex_performance_ci().to_string(index=False))
    print(rank_correlation_stats().to_string(index=False))
    print(tanimoto_performance_correlations().to_string(index=False))


if __name__ == "__main__":
    main()
