#!/usr/bin/env python3
"""Bootstrap seed-level confidence intervals for drug-similarity audit metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def percentile_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int, alpha: float) -> tuple[float, float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan, np.nan
    boot = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        boot[i] = rng.choice(values, size=values.size, replace=True).mean()
    low = float(np.quantile(boot, alpha / 2))
    high = float(np.quantile(boot, 1 - alpha / 2))
    return float(values.mean()), low, high


def paired_difference_ci(
    df: pd.DataFrame,
    split_a: str,
    split_b: str,
    baseline: str,
    metric: str,
    rng: np.random.Generator,
    n_boot: int,
    alpha: float,
) -> dict:
    wide = (
        df.loc[df["baseline"].eq(baseline) & df["split"].isin([split_a, split_b]), ["seed", "split", metric]]
        .pivot(index="seed", columns="split", values=metric)
        .dropna()
    )
    if wide.empty:
        return {
            "contrast": f"{split_a} minus {split_b}",
            "baseline": baseline,
            "metric": metric,
            "n_seeds": 0,
            "mean_difference": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
        }
    diffs = (wide[split_a] - wide[split_b]).to_numpy(dtype=np.float64)
    mean, low, high = percentile_ci(diffs, rng, n_boot, alpha)
    return {
        "contrast": f"{split_a} minus {split_b}",
        "baseline": baseline,
        "metric": metric,
        "n_seeds": int(diffs.size),
        "mean_difference": mean,
        "ci_low": low,
        "ci_high": high,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/drug_similarity_bootstrap_ci"))
    parser.add_argument("--baseline", default="ridge_cell_drug_fp")
    parser.add_argument("--n-boot", type=int, default=10000)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260622)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    df = pd.read_csv(args.summary)

    metrics = [
        "mean_max_train_tanimoto",
        "corr_tanimoto_pearson",
        "corr_tanimoto_rmse",
    ]
    rows = []
    for (split, baseline), group in df.groupby(["split", "baseline"]):
        if baseline != args.baseline:
            continue
        for metric in metrics:
            mean, low, high = percentile_ci(group[metric].to_numpy(dtype=np.float64), rng, args.n_boot, args.alpha)
            rows.append(
                {
                    "split": split,
                    "baseline": baseline,
                    "metric": metric,
                    "n_seeds": int(group[metric].notna().sum()),
                    "mean": mean,
                    "ci_low": low,
                    "ci_high": high,
                }
            )

    ci = pd.DataFrame(rows)
    ci.to_csv(args.outdir / "similarity_metric_seed_bootstrap_ci.csv", index=False)

    contrast_rows = []
    for metric in metrics:
        contrast_rows.append(
            paired_difference_ci(
                df,
                "split_random",
                "split_scaffold_heldout",
                args.baseline,
                metric,
                rng,
                args.n_boot,
                args.alpha,
            )
        )
        contrast_rows.append(
            paired_difference_ci(
                df,
                "split_random",
                "split_cell_scaffold_heldout",
                args.baseline,
                metric,
                rng,
                args.n_boot,
                args.alpha,
            )
        )
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(args.outdir / "similarity_metric_seed_bootstrap_contrasts.csv", index=False)
    print(ci.to_string(index=False))
    print(contrasts.to_string(index=False))


if __name__ == "__main__":
    main()
