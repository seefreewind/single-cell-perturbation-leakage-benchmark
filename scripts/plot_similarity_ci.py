#!/usr/bin/env python3
"""Plot bootstrap confidence intervals for nearest training-drug similarity."""

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
    parser.add_argument("--ci", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("figures/openproblems_similarity_ci"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.ci)
    df = df[df["metric"].eq("mean_max_train_tanimoto")].copy()
    df["split_label"] = pd.Categorical(
        df["split"].map(SPLIT_LABELS),
        categories=list(SPLIT_LABELS.values()),
        ordered=True,
    )
    df = df.sort_values("split_label")

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    sns.pointplot(data=df, x="split_label", y="mean", color="#2f6f73", markers="o", linestyles="", ax=ax)
    for i, row in enumerate(df.itertuples(index=False)):
        ax.plot([i, i], [row.ci_low, row.ci_high], color="#1f3f41", linewidth=2)
        ax.plot([i - 0.08, i + 0.08], [row.ci_low, row.ci_low], color="#1f3f41", linewidth=2)
        ax.plot([i - 0.08, i + 0.08], [row.ci_high, row.ci_high], color="#1f3f41", linewidth=2)
    ax.set_xlabel("")
    ax.set_ylabel("Nearest training-drug Tanimoto")
    ax.set_title("Scaffold-strict splits remove chemical-neighbor shortcuts")
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis="x", rotation=18)
    fig.tight_layout()
    fig.savefig(args.outdir / "nearest_train_drug_tanimoto_bootstrap_ci.png", dpi=300)
    fig.savefig(args.outdir / "nearest_train_drug_tanimoto_bootstrap_ci.pdf")
    print(args.outdir / "nearest_train_drug_tanimoto_bootstrap_ci.png")


if __name__ == "__main__":
    main()
