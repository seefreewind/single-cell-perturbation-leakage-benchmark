#!/usr/bin/env python3
"""Plot top-k DE overlap and direction agreement."""

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
    "nearest_drug_any_cell": "Nearest drug",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--top-k", type=int, default=100)
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems_topk"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.summary)
    df = df[df["top_k"].eq(args.top_k)].copy()
    df["split_label"] = pd.Categorical(
        df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )
    df["baseline_label"] = df["baseline"].map(BASELINE_LABELS).fillna(df["baseline"])
    df = df[df["baseline_label"].isin(BASELINE_LABELS.values())].copy()

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), sharex=True)
    sns.pointplot(
        data=df,
        x="split_label",
        y="mean_topk_overlap",
        hue="baseline_label",
        errorbar="sd",
        dodge=0.35,
        markers="o",
        linestyles="-",
        ax=axes[0],
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel(f"Top-{args.top_k} overlap")
    axes[0].tick_params(axis="x", rotation=18)
    axes[0].legend_.remove()

    sns.pointplot(
        data=df,
        x="split_label",
        y="mean_direction_agreement",
        hue="baseline_label",
        errorbar="sd",
        dodge=0.35,
        markers="o",
        linestyles="-",
        ax=axes[1],
    )
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Direction agreement among shared genes")
    axes[1].tick_params(axis="x", rotation=18)
    axes[1].legend(title="", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.suptitle(f"Top-{args.top_k} differential-response gene recovery")
    fig.tight_layout()
    fig.savefig(args.outdir / f"top{args.top_k}_de_overlap_direction.png", dpi=300)
    fig.savefig(args.outdir / f"top{args.top_k}_de_overlap_direction.pdf")
    print(args.outdir / f"top{args.top_k}_de_overlap_direction.png")


if __name__ == "__main__":
    main()

