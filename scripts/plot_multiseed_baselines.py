#!/usr/bin/env python3
"""Plot baseline performance across repeated split seeds."""

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
    "cell_context_mean": "Cell mean",
    "drug_mean": "Drug mean",
    "ridge_cell": "Ridge: cell",
    "ridge_drug_fp": "Ridge: drug FP",
    "ridge_cell_drug_fp": "Ridge: cell + drug FP",
    "nearest_drug_any_cell": "Nearest drug",
    "nearest_drug_same_cell": "Nearest drug, same cell",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems_multiseed"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.summary)
    df["split_label"] = pd.Categorical(
        df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )
    df["baseline_label"] = df["baseline"].map(BASELINE_LABELS).fillna(df["baseline"])

    keep = [
        "Global mean",
        "Ridge: cell + drug FP",
        "Nearest drug",
        "Nearest drug, same cell",
    ]
    plot_df = df[df["baseline_label"].isin(keep)].copy()

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.pointplot(
        data=plot_df,
        x="split_label",
        y="mean_rowwise_pearson",
        hue="baseline_label",
        errorbar="sd",
        dodge=0.35,
        markers="o",
        linestyles="-",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Mean row-wise Pearson across seeds")
    ax.set_title("Baseline performance declines under leakage-resistant validation")
    ax.tick_params(axis="x", rotation=18)
    ax.legend(title="", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(args.outdir / "multiseed_baseline_decay_mean_pearson.png", dpi=300)
    fig.savefig(args.outdir / "multiseed_baseline_decay_mean_pearson.pdf")
    print(args.outdir / "multiseed_baseline_decay_mean_pearson.png")


if __name__ == "__main__":
    main()
