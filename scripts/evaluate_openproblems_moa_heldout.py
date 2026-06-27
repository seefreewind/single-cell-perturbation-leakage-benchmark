#!/usr/bin/env python3
"""Evaluate OpenProblems baselines under a MoA-held-out secondary split."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

from evaluate_ridge_openproblems import (
    build_features,
    fit_predict,
    load_matrix,
    masked_rmse,
    masked_rowwise_pearson,
    rmse,
    rowwise_pearson,
)
from evaluate_simple_baselines_openproblems import (
    grouped_mean,
    load_combined,
    predict_grouped,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/openproblems_moa_heldout_100"


def normalize_name(value: object) -> str:
    return str(value).strip().lower()


def load_moa_map() -> pd.DataFrame:
    moa = pd.read_csv(ROOT / "data/raw/openproblems_neurips2023/moa_annotations.csv")
    moa["drug_key"] = moa["sm_name"].map(normalize_name)
    return moa.dropna(subset=["moa"]).drop_duplicates("drug_key")[["drug_key", "moa"]]


def build_seed_split(seed: int, holdout_fraction: float = 0.2) -> pd.DataFrame:
    source = pd.read_csv(ROOT / "results/openproblems_multiseed_100_de/seed_001/splits/candidate_splits.csv")
    source = source.copy()
    source["drug_key"] = source["drug_name"].map(normalize_name)
    source = source.merge(load_moa_map(), on="drug_key", how="left")
    moas = np.array(sorted(source.loc[source["condition"].eq("treated"), "moa"].dropna().unique()))
    rng = np.random.default_rng(seed + 20260624)
    n_holdout = max(1, int(round(holdout_fraction * len(moas))))
    heldout = set(rng.choice(moas, size=n_holdout, replace=False))
    source["split_moa_heldout"] = np.where(source["moa"].isin(heldout), "test", "train")
    source["heldout_moa"] = source["moa"].where(source["moa"].isin(heldout), "")
    return source


def evaluate_simple(seed: int, splits: pd.DataFrame, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix_meta, x, de_mask = load_combined(
        ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad",
        ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad",
        "clipped_sign_log10_pval",
        "is_de",
    )
    meta = matrix_meta.merge(
        splits[["sample_id", "drug_name", "cell_context", "condition", "split_moa_heldout"]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    treated = meta["condition"].eq("treated").to_numpy()
    meta = meta.loc[treated].reset_index(drop=True)
    x = x[treated]
    de_mask = de_mask[treated]

    train_mask = meta["split_moa_heldout"].eq("train").to_numpy()
    test_mask = meta["split_moa_heldout"].eq("test").to_numpy()
    train_meta = meta.loc[train_mask].reset_index(drop=True)
    test_meta = meta.loc[test_mask].reset_index(drop=True)
    train_x = x[train_mask]
    test_x = x[test_mask]
    test_de_mask = de_mask[test_mask]
    global_mean = train_x.mean(axis=0)
    baselines = {
        "global_train_mean": np.tile(global_mean, (len(test_meta), 1)),
        "cell_context_mean": predict_grouped(
            test_meta, grouped_mean(train_meta, train_x, "cell_context"), "cell_context", global_mean
        ),
        "drug_mean": predict_grouped(test_meta, grouped_mean(train_meta, train_x, "drug_name"), "drug_name", global_mean),
    }
    rows = []
    per_rows = []
    for baseline, pred in baselines.items():
        corr = rowwise_pearson(test_x, pred)
        err = rmse(test_x, pred)
        de_corr = masked_rowwise_pearson(test_x, pred, test_de_mask, 10)
        de_err = masked_rmse(test_x, pred, test_de_mask, 10)
        rows.append(
            {
                "seed": seed,
                "split": "split_moa_heldout",
                "baseline": baseline,
                "n_train": len(train_meta),
                "n_test": len(test_meta),
                "mean_rowwise_pearson": float(np.nanmean(corr)),
                "mean_rmse": float(np.nanmean(err)),
                "mean_de_rowwise_pearson": float(np.nanmean(de_corr)),
                "mean_de_rmse": float(np.nanmean(de_err)),
                "n_de_evaluable": int(np.isfinite(de_corr).sum()),
            }
        )
        per_rows.extend(
            {
                "seed": seed,
                "split": "split_moa_heldout",
                "baseline": baseline,
                "sample_id": sample_id,
                "rowwise_pearson": c,
                "rmse": e,
                "de_rowwise_pearson": dc,
                "de_rmse": de,
            }
            for sample_id, c, e, dc, de in zip(test_meta["sample_id"], corr, err, de_corr, de_err)
        )
    return pd.DataFrame(rows), pd.DataFrame(per_rows)


def evaluate_ridge(seed: int, splits: pd.DataFrame, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix_meta, y, de_mask = load_matrix(
        ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad",
        ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad",
        "clipped_sign_log10_pval",
        "is_de",
    )
    meta = matrix_meta.merge(
        splits[["sample_id", "drug_name", "cell_context", "condition", "smiles", "split_moa_heldout"]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    treated = meta["condition"].eq("treated").to_numpy()
    meta = meta.loc[treated].reset_index(drop=True)
    y = y[treated]
    de_mask = de_mask[treated]

    train_mask = meta["split_moa_heldout"].eq("train").to_numpy()
    test_mask = meta["split_moa_heldout"].eq("test").to_numpy()
    train_meta = meta.loc[train_mask].reset_index(drop=True)
    test_meta = meta.loc[test_mask].reset_index(drop=True)
    train_y = y[train_mask]
    test_y = y[test_mask]
    test_de_mask = de_mask[test_mask]
    rows = []
    per_rows = []
    for feature_set in ["cell", "drug_fp", "cell_drug_fp"]:
        pred = fit_predict(train_meta, test_meta, train_y, feature_set, 10.0, 1024)
        corr = rowwise_pearson(test_y, pred)
        err = rmse(test_y, pred)
        de_corr = masked_rowwise_pearson(test_y, pred, test_de_mask, 10)
        de_err = masked_rmse(test_y, pred, test_de_mask, 10)
        baseline = f"ridge_{feature_set}"
        rows.append(
            {
                "seed": seed,
                "split": "split_moa_heldout",
                "baseline": baseline,
                "n_train": len(train_meta),
                "n_test": len(test_meta),
                "alpha": 10.0,
                "n_bits": 1024,
                "mean_rowwise_pearson": float(np.nanmean(corr)),
                "mean_rmse": float(np.nanmean(err)),
                "mean_de_rowwise_pearson": float(np.nanmean(de_corr)),
                "mean_de_rmse": float(np.nanmean(de_err)),
                "n_de_evaluable": int(np.isfinite(de_corr).sum()),
            }
        )
        per_rows.extend(
            {
                "seed": seed,
                "split": "split_moa_heldout",
                "baseline": baseline,
                "sample_id": sample_id,
                "rowwise_pearson": c,
                "rmse": e,
                "de_rowwise_pearson": dc,
                "de_rmse": de,
            }
            for sample_id, c, e, dc, de in zip(test_meta["sample_id"], corr, err, de_corr, de_err)
        )
    return pd.DataFrame(rows), pd.DataFrame(per_rows)


def leakage_summary(seed: int, splits: pd.DataFrame) -> dict[str, float | int]:
    treated = splits.loc[splits["condition"].eq("treated")].copy()
    train = treated.loc[treated["split_moa_heldout"].eq("train")]
    test = treated.loc[treated["split_moa_heldout"].eq("test")]
    train_drugs = set(train["drug_name"])
    train_scaffolds = set(train["scaffold"])
    train_moas = set(train["moa"].dropna())
    return {
        "seed": seed,
        "split": "split_moa_heldout",
        "n_train": len(train),
        "n_test": len(test),
        "n_heldout_moa": test["moa"].nunique(),
        "n_test_drugs": test["drug_name"].nunique(),
        "n_test_scaffolds": test["scaffold"].nunique(),
        "same_drug_in_train": float(test["drug_name"].isin(train_drugs).mean()),
        "same_scaffold_in_train": float(test["scaffold"].isin(train_scaffolds).mean()),
        "same_moa_in_train": float(test["moa"].isin(train_moas).mean()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    per_rows = []
    leak_rows = []
    for seed in range(1, 101):
        seed_dir = OUT / f"seed_{seed:03d}"
        (seed_dir / "splits").mkdir(parents=True, exist_ok=True)
        splits = build_seed_split(seed)
        splits.to_csv(seed_dir / "splits/candidate_splits.csv", index=False)
        simple_summary, simple_per = evaluate_simple(seed, splits, seed_dir)
        ridge_summary, ridge_per = evaluate_ridge(seed, splits, seed_dir)
        simple_summary.to_csv(seed_dir / "simple_baseline_summary.csv", index=False)
        ridge_summary.to_csv(seed_dir / "ridge_baseline_summary.csv", index=False)
        pd.concat([simple_per, ridge_per], ignore_index=True).to_csv(seed_dir / "baseline_per_row.csv", index=False)
        summaries.extend([simple_summary, ridge_summary])
        per_rows.extend([simple_per, ridge_per])
        leak_rows.append(leakage_summary(seed, splits))
        if seed % 10 == 0:
            print(f"completed seed {seed}")

    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(OUT / "all_baseline_summary.csv", index=False)
    per_row = pd.concat(per_rows, ignore_index=True)
    per_row.to_csv(OUT / "baseline_per_row_all.csv", index=False)
    leakage = pd.DataFrame(leak_rows)
    leakage.to_csv(OUT / "leakage_summary_by_seed.csv", index=False)

    agg = (
        summary.groupby(["split", "baseline"], as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_train_records=("n_train", "mean"),
            mean_test_records=("n_test", "mean"),
            mean_rowwise_pearson=("mean_rowwise_pearson", "mean"),
            sd_rowwise_pearson=("mean_rowwise_pearson", "std"),
            mean_de_rowwise_pearson=("mean_de_rowwise_pearson", "mean"),
            sd_de_rowwise_pearson=("mean_de_rowwise_pearson", "std"),
        )
    )
    agg.to_csv(OUT / "baseline_summary_aggregate.csv", index=False)
    leak_agg = leakage.agg(
        {
            "n_train": ["mean", "std"],
            "n_test": ["mean", "std"],
            "n_heldout_moa": ["mean", "std"],
            "n_test_drugs": ["mean", "std"],
            "n_test_scaffolds": ["mean", "std"],
            "same_drug_in_train": "mean",
            "same_scaffold_in_train": "mean",
            "same_moa_in_train": "mean",
        }
    )
    leak_agg.to_csv(OUT / "leakage_summary_aggregate_matrix.csv")
    print(agg.to_string(index=False))
    print(leakage.mean(numeric_only=True).to_string())


if __name__ == "__main__":
    main()
