#!/usr/bin/env python3
"""Run multi-seed Sci-Plex 3 split, baseline, leakage, and ranking audits."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent))

from drug_similarity_audit import fp_from_smiles, max_tanimoto  # noqa: E402
from evaluate_sciplex3_response_baselines import mean_predict, ridge_predict  # noqa: E402
from evaluate_ridge_openproblems import rmse, rowwise_pearson  # noqa: E402
from make_split_plan import make_splits  # noqa: E402


SPLIT_COLS = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]

MEAN_SPECS = [
    ("global_train_mean", None),
    ("cell_context_mean", "cell_context"),
    ("drug_mean", "drug_name"),
    ("cell_dose_mean", "cell_dose_key"),
]

RIDGE_SPECS = [
    ("ridge_cell_dose", "cell+dose"),
    ("ridge_drug_fp", "drug_fp"),
    ("ridge_cell_dose_drug_fp", "cell+dose+drug_fp"),
]


def split_label(split_col: str) -> str:
    return split_col.replace("split_", "")


def summarize_overlap(splits: pd.DataFrame, seed: int) -> list[dict]:
    rows = []
    for split_col in SPLIT_COLS:
        train = splits[splits[split_col].eq("train")]
        test = splits[splits[split_col].eq("test")]
        train_drugs = set(train["drug_name"].dropna())
        train_scaffolds = set(train["scaffold"].dropna())
        train_cells = set(train["cell_context"].dropna())
        train_pairs = set(zip(train["drug_name"], train["cell_context"], strict=False))
        test_pairs = list(zip(test["drug_name"], test["cell_context"], strict=False))
        rows.append(
            {
                "seed": seed,
                "split": split_col,
                "n_train": len(train),
                "n_test": len(test),
                "n_excluded": int(splits[split_col].eq("excluded").sum()),
                "same_drug_in_train": test["drug_name"].isin(train_drugs).mean(),
                "same_scaffold_in_train": test["scaffold"].isin(train_scaffolds).mean(),
                "same_cell_in_train": test["cell_context"].isin(train_cells).mean(),
                "same_drug_cell_in_train": np.mean([pair in train_pairs for pair in test_pairs]) if len(test_pairs) else np.nan,
            }
        )
    return rows


def summarize_similarity(splits: pd.DataFrame, seed: int, fp_map: dict[str, object]) -> tuple[list[dict], pd.DataFrame]:
    rows = []
    per_row = []
    for split_col in SPLIT_COLS:
        train_drugs = sorted(splits.loc[splits[split_col].eq("train"), "drug_name"].dropna().unique())
        train_fps = [fp_map[d] for d in train_drugs if fp_map.get(d) is not None]
        test = splits.loc[splits[split_col].eq("test"), ["sample_id", "drug_name", "cell_context", "scaffold"]].copy()
        for row in test.itertuples(index=False):
            sim = max_tanimoto(fp_map.get(row.drug_name), train_fps)
            per_row.append(
                {
                    "seed": seed,
                    "split": split_col,
                    "sample_id": row.sample_id,
                    "drug_name": row.drug_name,
                    "cell_context": row.cell_context,
                    "scaffold": row.scaffold,
                    "max_train_tanimoto": sim,
                }
            )
        rows.append(
            {
                "seed": seed,
                "split": split_col,
                "mean_max_train_tanimoto": float(np.nanmean([r["max_train_tanimoto"] for r in per_row if r["seed"] == seed and r["split"] == split_col])),
                "n_test": len(test),
            }
        )
    return rows, pd.DataFrame(per_row)


def evaluate_seed(
    seed: int,
    splits: pd.DataFrame,
    meta_base: pd.DataFrame,
    y: np.ndarray,
    outdir: Path,
    alpha: float,
    n_bits: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    seed_dir = outdir / f"seed_{seed:03d}" / "baselines"
    seed_dir.mkdir(parents=True, exist_ok=True)
    meta = meta_base.merge(
        splits[["sample_id", "drug_name", "cell_context", "dose_value", "smiles", "scaffold", *SPLIT_COLS]],
        on="sample_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_split"),
    )
    meta["cell_dose_key"] = meta["cell_context"].astype(str) + "||" + meta["dose_value"].astype(str)

    summary_rows = []
    per_row = []
    for split_col in SPLIT_COLS:
        train_mask = meta[split_col].eq("train").to_numpy()
        test_mask = meta[split_col].eq("test").to_numpy()
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_y = y[train_mask]
        test_y = y[test_mask]
        if len(train_meta) == 0 or len(test_meta) == 0:
            continue

        predictions = []
        for baseline, key in MEAN_SPECS:
            predictions.append((baseline, mean_predict(train_meta, train_y, test_meta, key)))
        for baseline, feature_set in RIDGE_SPECS:
            predictions.append((baseline, ridge_predict(train_meta, test_meta, train_y, feature_set, alpha, n_bits)))

        for baseline, pred in predictions:
            corr = rowwise_pearson(test_y, pred)
            err = rmse(test_y, pred)
            summary_rows.append(
                {
                    "seed": seed,
                    "split": split_col,
                    "split_label": split_label(split_col),
                    "baseline": baseline,
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "mean_rowwise_pearson": float(np.nanmean(corr)),
                    "median_rowwise_pearson": float(np.nanmedian(corr)),
                    "mean_rmse": float(np.nanmean(err)),
                    "median_rmse": float(np.nanmedian(err)),
                    "alpha": alpha if baseline.startswith("ridge") else np.nan,
                    "n_bits": n_bits if "drug_fp" in baseline else np.nan,
                }
            )
            per_row.extend(
                {
                    "seed": seed,
                    "split": split_col,
                    "baseline": baseline,
                    "sample_id": sample_id,
                    "rowwise_pearson": c,
                    "rmse": e,
                }
                for sample_id, c, e in zip(test_meta["sample_id"], corr, err, strict=False)
            )

    summary = pd.DataFrame(summary_rows)
    per_row_df = pd.DataFrame(per_row)
    summary.to_csv(seed_dir / "baseline_summary.csv", index=False)
    per_row_df.to_csv(seed_dir / "baseline_per_row.csv", index=False)
    return summary, per_row_df


def ranking_outputs(summary: pd.DataFrame, outdir: Path) -> None:
    ranked = summary.copy()
    ranked["rank"] = ranked.groupby(["seed", "split"])["mean_rowwise_pearson"].rank(ascending=False, method="min")
    ranked.to_csv(outdir / "baseline_ranks_by_seed.csv", index=False)
    mean_rank = (
        ranked.groupby(["split", "baseline"])
        .agg(
            mean_rank=("rank", "mean"),
            sd_rank=("rank", "std"),
            mean_score=("mean_rowwise_pearson", "mean"),
            sd_score=("mean_rowwise_pearson", "std"),
            n_seeds=("seed", "nunique"),
        )
        .reset_index()
        .sort_values(["split", "mean_rank", "mean_score"], ascending=[True, True, False])
    )
    mean_rank.to_csv(outdir / "mean_rank_summary.csv", index=False)

    correlations = []
    for seed, seed_df in ranked.groupby("seed"):
        ref = seed_df[seed_df["split"].eq("split_random")][["baseline", "rank"]].rename(columns={"rank": "random_rank"})
        for split, split_df in seed_df.groupby("split"):
            if split == "split_random":
                continue
            merged = split_df[["baseline", "rank"]].merge(ref, on="baseline", how="inner")
            correlations.append(
                {
                    "seed": seed,
                    "comparison_split": split,
                    "spearman_rank_correlation": merged["rank"].corr(merged["random_rank"], method="spearman"),
                    "n_baselines": len(merged),
                }
            )
    corr = pd.DataFrame(correlations)
    corr.to_csv(outdir / "split_rank_correlations.csv", index=False)
    (
        corr.groupby("comparison_split")
        .agg(
            mean_spearman=("spearman_rank_correlation", "mean"),
            sd_spearman=("spearman_rank_correlation", "std"),
            n_seeds=("seed", "nunique"),
        )
        .reset_index()
        .to_csv(outdir / "rank_correlation_summary.csv", index=False)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=Path("data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad"))
    parser.add_argument("--metadata", type=Path, default=Path("metadata/sciplex3_24h_top2000_split_metadata.csv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/sciplex3_24h_top2000_multiseed_30"))
    parser.add_argument("--n-seeds", type=int, default=30)
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    meta_base = adata.obs.copy().reset_index(names="sample_id")
    y = np.asarray(adata.X, dtype=np.float64)
    drug_table = (
        pd.read_csv(args.metadata)[["drug_name", "smiles"]]
        .drop_duplicates("drug_name")
        .assign(fp=lambda d: d["smiles"].map(fp_from_smiles))
    )
    fp_map = dict(zip(drug_table["drug_name"], drug_table["fp"], strict=False))

    summaries = []
    per_rows = []
    split_summaries = []
    overlap_rows = []
    similarity_rows = []
    similarity_per_rows = []

    for seed in range(args.seed_start, args.seed_start + args.n_seeds):
        seed_root = args.outdir / f"seed_{seed:03d}"
        split_dir = seed_root / "splits"
        make_splits(args.metadata, split_dir, args.test_fraction, seed)
        splits = pd.read_csv(split_dir / "candidate_splits.csv")
        split_summary = pd.read_csv(split_dir / "candidate_split_summary.csv")
        split_summary.insert(0, "seed", seed)
        split_summaries.append(split_summary)

        overlap_rows.extend(summarize_overlap(splits, seed))
        sim_rows, sim_per = summarize_similarity(splits, seed, fp_map)
        similarity_rows.extend(sim_rows)
        similarity_per_rows.append(sim_per)

        summary, per_row = evaluate_seed(seed, splits, meta_base, y, args.outdir, args.alpha, args.n_bits)
        summaries.append(summary)
        per_rows.append(per_row)
        print(f"finished seed {seed:03d}")

    all_summary = pd.concat(summaries, ignore_index=True)
    all_per_row = pd.concat(per_rows, ignore_index=True)
    split_summary_all = pd.concat(split_summaries, ignore_index=True)
    overlap = pd.DataFrame(overlap_rows)
    similarity = pd.DataFrame(similarity_rows)
    similarity_per = pd.concat(similarity_per_rows, ignore_index=True)

    all_summary.to_csv(args.outdir / "baseline_summary_by_seed.csv", index=False)
    all_per_row.to_csv(args.outdir / "baseline_per_row_all.csv", index=False)
    split_summary_all.to_csv(args.outdir / "split_summary_by_seed.csv", index=False)
    overlap.to_csv(args.outdir / "leakage_overlap_by_seed.csv", index=False)
    similarity.to_csv(args.outdir / "similarity_by_seed.csv", index=False)
    similarity_per.to_csv(args.outdir / "similarity_per_row_all.csv", index=False)

    (
        all_summary.groupby(["split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_test=("n_test", "mean"),
            mean_rowwise_pearson=("mean_rowwise_pearson", "mean"),
            sd_rowwise_pearson=("mean_rowwise_pearson", "std"),
            mean_rmse=("mean_rmse", "mean"),
            sd_rmse=("mean_rmse", "std"),
        )
        .reset_index()
        .to_csv(args.outdir / "baseline_summary_aggregate.csv", index=False)
    )
    (
        split_summary_all.groupby("split")
        .agg(
            mean_train_records=("train_records", "mean"),
            sd_train_records=("train_records", "std"),
            mean_test_records=("test_records", "mean"),
            sd_test_records=("test_records", "std"),
            mean_excluded_records=("excluded_records", "mean"),
            sd_excluded_records=("excluded_records", "std"),
        )
        .reset_index()
        .to_csv(args.outdir / "split_summary_aggregate.csv", index=False)
    )
    (
        overlap.groupby("split")
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_test=("n_test", "mean"),
            same_drug_in_train=("same_drug_in_train", "mean"),
            same_scaffold_in_train=("same_scaffold_in_train", "mean"),
            same_cell_in_train=("same_cell_in_train", "mean"),
            same_drug_cell_in_train=("same_drug_cell_in_train", "mean"),
        )
        .reset_index()
        .to_csv(args.outdir / "leakage_overlap_summary.csv", index=False)
    )
    (
        similarity.groupby("split")
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_test=("n_test", "mean"),
            mean_max_train_tanimoto=("mean_max_train_tanimoto", "mean"),
            sd_max_train_tanimoto=("mean_max_train_tanimoto", "std"),
        )
        .reset_index()
        .to_csv(args.outdir / "similarity_summary.csv", index=False)
    )
    ranking_outputs(all_summary, args.outdir)

    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
