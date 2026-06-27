#!/usr/bin/env python3
"""Plot direct ridge versus SVD-ridge under repeated splits."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


SPLIT_LABELS = {
    "split_random": "Random",
    "split_cell_heldout": "Cell held-out",
    "split_scaffold_heldout": "Scaffold held-out",
    "split_cell_scaffold_heldout": "Cell + scaffold held-out",
}

BASELINE_LABELS = {
    "global_train_mean": "Global mean",
    "ridge_cell_drug_fp": "Ridge: cell + drug FP",
    "svd50_ridge_cell_drug_fp": "SVD50 + ridge: cell + drug FP",
    "nearest_drug_any_cell": "Nearest drug",
}


def load_summary(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        frame = pd.read_csv(path)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def plot_metric(df: pd.DataFrame, metric: str, ylabel: str, outfile_stem: Path) -> None:
    plot_df = df[df["baseline"].isin(BASELINE_LABELS)].copy()
    plot_df["split_label"] = pd.Categorical(
        plot_df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )
    plot_df["baseline_label"] = pd.Categorical(
        plot_df["baseline"].map(BASELINE_LABELS),
        categories=list(BASELINE_LABELS.values()),
        ordered=True,
    )

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(13, 6))
    sns.pointplot(
        data=plot_df,
        x="split_label",
        y=metric,
        hue="baseline_label",
        errorbar="sd",
        dodge=0.42,
        markers="o",
        linestyles="-",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title("Stronger baselines still degrade under strict held-out validation")
    ax.tick_params(axis="x", rotation=18)
    ax.legend(title="", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(outfile_stem.with_suffix(".png"), dpi=300)
    fig.savefig(outfile_stem.with_suffix(".pdf"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summaries", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems_svd_ridge"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = load_summary(args.summaries)
    plot_metric(
        df,
        "mean_rowwise_pearson",
        "Mean row-wise Pearson across seeds",
        args.outdir / "svd_ridge_comparison_all_genes",
    )
    plot_metric(
        df,
        "mean_de_rowwise_pearson",
        "DE-gene row-wise Pearson across seeds",
        args.outdir / "svd_ridge_comparison_de_genes",
    )
    df.to_csv(args.outdir / "svd_ridge_comparison_plot_data.csv", index=False)
    print(args.outdir / "svd_ridge_comparison_all_genes.png")
    print(args.outdir / "svd_ridge_comparison_de_genes.png")


if __name__ == "__main__":
    main()
