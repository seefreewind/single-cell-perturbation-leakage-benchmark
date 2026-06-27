#!/usr/bin/env python3
"""Plot drug similarity audit results."""

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
    parser.add_argument("--per-row", type=Path, required=True)
    parser.add_argument("--baseline", default="ridge_cell_drug_fp")
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.per_row)
    df = df[df["baseline"].eq(args.baseline)].copy()
    df["split_label"] = pd.Categorical(
        df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    sns.boxplot(data=df, x="split_label", y="max_train_tanimoto", ax=axes[0], color="#8fb9aa")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Nearest training-drug Tanimoto")
    axes[0].tick_params(axis="x", rotation=20)

    plot_df = df[df["split"].isin(["split_random", "split_scaffold_heldout", "split_cell_scaffold_heldout"])]
    sns.scatterplot(
        data=plot_df,
        x="max_train_tanimoto",
        y="rowwise_pearson",
        hue="split_label",
        alpha=0.8,
        ax=axes[1],
    )
    axes[1].set_xlabel("Nearest training-drug Tanimoto")
    axes[1].set_ylabel("Row-wise Pearson")
    axes[1].legend(title="", frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(args.outdir / "drug_similarity_audit_ridge_cell_drug_fp.png", dpi=300)
    fig.savefig(args.outdir / "drug_similarity_audit_ridge_cell_drug_fp.pdf")
    print(args.outdir / "drug_similarity_audit_ridge_cell_drug_fp.png")


if __name__ == "__main__":
    main()

