#!/usr/bin/env python3
"""Audit whether model rankings from random splits transfer to strict splits."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_SPLIT_ORDER = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]


def rank_within_seed(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    ranked = df.copy()
    ranked["rank"] = ranked.groupby(["seed", "split"])[metric].rank(
        ascending=False, method="min"
    )
    ranked["n_ranked"] = ranked.groupby(["seed", "split"])["baseline"].transform("nunique")
    return ranked


def winner_frequency(ranked: pd.DataFrame) -> pd.DataFrame:
    winners = ranked[ranked["rank"].eq(1)].copy()
    return (
        winners.groupby(["split", "baseline"])
        .agg(
            n_wins=("seed", "nunique"),
            mean_score=("score", "mean"),
            sd_score=("score", "std"),
        )
        .reset_index()
        .sort_values(["split", "n_wins", "mean_score"], ascending=[True, False, False])
    )


def mean_rank_table(ranked: pd.DataFrame) -> pd.DataFrame:
    return (
        ranked.groupby(["split", "baseline"])
        .agg(
            mean_rank=("rank", "mean"),
            sd_rank=("rank", "std"),
            mean_score=("score", "mean"),
            sd_score=("score", "std"),
            n_seeds=("seed", "nunique"),
        )
        .reset_index()
        .sort_values(["split", "mean_rank", "mean_score"], ascending=[True, True, False])
    )


def rank_transfer(ranked: pd.DataFrame, reference_split: str) -> pd.DataFrame:
    ref = ranked[ranked["split"].eq(reference_split)][
        ["seed", "baseline", "rank", "score"]
    ].rename(columns={"rank": "reference_rank", "score": "reference_score"})
    rows = ranked[~ranked["split"].eq(reference_split)].merge(
        ref, on=["seed", "baseline"], how="inner"
    )
    rows["rank_shift_vs_reference"] = rows["rank"] - rows["reference_rank"]
    rows["score_delta_vs_reference"] = rows["score"] - rows["reference_score"]
    return rows


def split_rank_correlations(ranked: pd.DataFrame, reference_split: str) -> pd.DataFrame:
    rows = []
    for seed, seed_df in ranked.groupby("seed"):
        ref = seed_df[seed_df["split"].eq(reference_split)][["baseline", "rank"]].rename(
            columns={"rank": "reference_rank"}
        )
        for split, split_df in seed_df.groupby("split"):
            if split == reference_split:
                continue
            merged = split_df[["baseline", "rank"]].merge(ref, on="baseline", how="inner")
            if len(merged) < 2:
                corr = float("nan")
            else:
                corr = merged["rank"].corr(merged["reference_rank"], method="spearman")
            rows.append(
                {
                    "seed": seed,
                    "reference_split": reference_split,
                    "comparison_split": split,
                    "n_baselines": len(merged),
                    "spearman_rank_correlation": corr,
                }
            )
    return pd.DataFrame(rows)


def write_markdown(
    outpath: Path,
    metric: str,
    mean_ranks: pd.DataFrame,
    winners: pd.DataFrame,
    transfer: pd.DataFrame,
    correlations: pd.DataFrame,
) -> None:
    strict_splits = [
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]
    top_lines = []
    for split in DEFAULT_SPLIT_ORDER:
        table = mean_ranks[mean_ranks["split"].eq(split)].head(3)
        if table.empty:
            continue
        names = ", ".join(
            f"{r.baseline} (rank {r.mean_rank:.2f}, score {r.mean_score:.3f})"
            for r in table.itertuples()
        )
        top_lines.append(f"- {split}: {names}")

    corr_lines = []
    corr_summary = (
        correlations.groupby("comparison_split")
        .agg(
            mean_spearman=("spearman_rank_correlation", "mean"),
            sd_spearman=("spearman_rank_correlation", "std"),
            n_seeds=("seed", "nunique"),
        )
        .reset_index()
    )
    for row in corr_summary.itertuples():
        corr_lines.append(
            f"- random vs {row.comparison_split}: Spearman rho = "
            f"{row.mean_spearman:.3f} +/- {row.sd_spearman:.3f} across {row.n_seeds} seeds"
        )

    transfer_lines = []
    for split in strict_splits:
        rows = (
            transfer[transfer["split"].eq(split)]
            .groupby("baseline")
            .agg(
                mean_rank_shift=("rank_shift_vs_reference", "mean"),
                mean_score_delta=("score_delta_vs_reference", "mean"),
                n_seeds=("seed", "nunique"),
            )
            .reset_index()
            .sort_values("mean_score_delta")
        )
        if rows.empty:
            continue
        worst = rows.head(3)
        vals = ", ".join(
            f"{r.baseline} (Delta score {r.mean_score_delta:.3f}, Delta rank {r.mean_rank_shift:.2f})"
            for r in worst.itertuples()
        )
        transfer_lines.append(f"- {split}: {vals}")

    text = "\n".join(
        [
            "# Ranking instability audit",
            "",
            f"Metric: `{metric}`. Higher values are ranked better within each seed and split.",
            "",
            "## Mean top-ranked baselines",
            *top_lines,
            "",
            "## Rank transfer from random split",
            *corr_lines,
            "",
            "## Largest score declines relative to random split",
            *transfer_lines,
            "",
            "Interpretation: random-split rankings should not be treated as a proxy for scaffold-strict or joint cell+scaffold generalization when rank correlations are weak and the largest score declines concentrate in drug-fingerprint-aware models.",
            "",
        ]
    )
    outpath.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/ranking_instability"))
    parser.add_argument("--metric", default="mean_de_rowwise_pearson")
    parser.add_argument("--reference-split", default="split_random")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.summary)
    required = {"seed", "split", "baseline", args.metric}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    df = df.dropna(subset=[args.metric]).copy()
    df["score"] = df[args.metric].astype(float)
    ranked = rank_within_seed(df, args.metric)
    mean_ranks = mean_rank_table(ranked)
    winners = winner_frequency(ranked)
    transfer = rank_transfer(ranked, args.reference_split)
    correlations = split_rank_correlations(ranked, args.reference_split)

    ranked.to_csv(args.outdir / "baseline_ranks_by_seed.csv", index=False)
    mean_ranks.to_csv(args.outdir / "mean_rank_summary.csv", index=False)
    winners.to_csv(args.outdir / "winner_frequency.csv", index=False)
    transfer.to_csv(args.outdir / "rank_transfer_from_random.csv", index=False)
    correlations.to_csv(args.outdir / "split_rank_correlations.csv", index=False)
    write_markdown(
        args.outdir / "ranking_instability_summary.md",
        args.metric,
        mean_ranks,
        winners,
        transfer,
        correlations,
    )

    print(mean_ranks.to_string(index=False))
    print()
    print(correlations.groupby("comparison_split")["spearman_rank_correlation"].agg(["mean", "std", "count"]).to_string())


if __name__ == "__main__":
    main()
