#!/usr/bin/env python3
"""Plot simple baseline performance across leakage-resistant split levels."""

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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.summary)
    df["split_label"] = pd.Categorical(
        df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )
    df["baseline_label"] = df["baseline"].str.replace("_", " ").str.title()

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(11, 5.5))
    sns.lineplot(
        data=df,
        x="split_label",
        y="mean_rowwise_pearson",
        hue="baseline_label",
        marker="o",
        linewidth=2,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Mean row-wise Pearson")
    ax.set_title("Simple baseline performance under increasingly strict splits")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(title="", frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(args.outdir / "simple_baseline_decay_mean_pearson.png", dpi=300)
    fig.savefig(args.outdir / "simple_baseline_decay_mean_pearson.pdf")
    print(args.outdir / "simple_baseline_decay_mean_pearson.png")


if __name__ == "__main__":
    main()

