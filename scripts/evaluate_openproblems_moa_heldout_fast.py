#!/usr/bin/env python3
"""Fast MoA-held-out baseline evaluation for OpenProblems."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_ridge_openproblems import (
    fit_predict,
    load_matrix,
    masked_rmse,
    masked_rowwise_pearson,
    rmse,
    rowwise_pearson,
)
from evaluate_simple_baselines_openproblems import grouped_mean, predict_grouped


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/openproblems_moa_heldout_100_fast"


def normalize_name(value: object) -> str:
    return str(value).strip().lower()


def load_base() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    matrix_meta, y, de_mask = load_matrix(
        ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad",
        ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad",
        "clipped_sign_log10_pval",
        "is_de",
    )
    split_meta = pd.read_csv(ROOT / "results/openproblems_multiseed_100_de/seed_001/splits/candidate_splits.csv")
    moa = pd.read_csv(ROOT / "data/raw/openproblems_neurips2023/moa_annotations.csv")
    moa["drug_key"] = moa["sm_name"].map(normalize_name)
    moa = moa.dropna(subset=["moa"]).drop_duplicates("drug_key")[["drug_key", "moa"]]
    split_meta = split_meta.copy()
    split_meta["drug_key"] = split_meta["drug_name"].map(normalize_name)
    split_meta = split_meta.merge(moa, on="drug_key", how="left")
    meta = matrix_meta.merge(
        split_meta[["sample_id", "drug_name", "cell_context", "condition", "smiles", "scaffold", "moa"]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    treated = meta["condition"].eq("treated").to_numpy()
    return meta.loc[treated].reset_index(drop=True), y[treated], de_mask[treated]


def make_assignment(meta: pd.DataFrame, seed: int, holdout_fraction: float = 0.2) -> np.ndarray:
    moas = np.array(sorted(meta["moa"].dropna().unique()))
    rng = np.random.default_rng(seed + 20260624)
    n_holdout = max(1, int(round(holdout_fraction * len(moas))))
    heldout = set(rng.choice(moas, size=n_holdout, replace=False))
    return np.where(meta["moa"].isin(heldout), "test", "train")


def eval_prediction(
    seed: int,
    baseline: str,
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    test_y: np.ndarray,
    pred: np.ndarray,
    test_de_mask: np.ndarray,
) -> dict[str, float | int | str]:
    corr = rowwise_pearson(test_y, pred)
    err = rmse(test_y, pred)
    de_corr = masked_rowwise_pearson(test_y, pred, test_de_mask, 10)
    de_err = masked_rmse(test_y, pred, test_de_mask, 10)
    return {
        "seed": seed,
        "split": "split_moa_heldout",
        "baseline": baseline,
        "n_train": len(train_meta),
        "n_test": len(test_meta),
        "mean_rowwise_pearson": float(np.nanmean(corr)),
        "sd_rowwise_pearson_per_sample": float(np.nanstd(corr, ddof=1)),
        "mean_rmse": float(np.nanmean(err)),
        "mean_de_rowwise_pearson": float(np.nanmean(de_corr)),
        "sd_de_rowwise_pearson_per_sample": float(np.nanstd(de_corr, ddof=1)),
        "mean_de_rmse": float(np.nanmean(de_err)),
        "n_de_evaluable": int(np.isfinite(de_corr).sum()),
    }


def leakage_row(seed: int, train_meta: pd.DataFrame, test_meta: pd.DataFrame) -> dict[str, float | int | str]:
    train_drugs = set(train_meta["drug_name"])
    train_scaffolds = set(train_meta["scaffold"])
    train_moas = set(train_meta["moa"].dropna())
    return {
        "seed": seed,
        "split": "split_moa_heldout",
        "n_train": len(train_meta),
        "n_test": len(test_meta),
        "n_heldout_moa": test_meta["moa"].nunique(),
        "n_test_drugs": test_meta["drug_name"].nunique(),
        "n_test_scaffolds": test_meta["scaffold"].nunique(),
        "same_drug_in_train": float(test_meta["drug_name"].isin(train_drugs).mean()),
        "same_scaffold_in_train": float(test_meta["scaffold"].isin(train_scaffolds).mean()),
        "same_moa_in_train": float(test_meta["moa"].isin(train_moas).mean()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta, y, de_mask = load_base()
    summary_rows = []
    leakage_rows = []
    for seed in range(1, 101):
        assignment = make_assignment(meta, seed)
        train_mask = assignment == "train"
        test_mask = assignment == "test"
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_y = y[train_mask]
        test_y = y[test_mask]
        test_de_mask = de_mask[test_mask]

        global_mean = train_y.mean(axis=0)
        simple_preds = {
            "global_train_mean": np.tile(global_mean, (len(test_meta), 1)),
            "cell_context_mean": predict_grouped(
                test_meta, grouped_mean(train_meta, train_y, "cell_context"), "cell_context", global_mean
            ),
            "drug_mean": predict_grouped(test_meta, grouped_mean(train_meta, train_y, "drug_name"), "drug_name", global_mean),
        }
        for baseline, pred in simple_preds.items():
            summary_rows.append(eval_prediction(seed, baseline, train_meta, test_meta, test_y, pred, test_de_mask))

        for feature_set in ["cell", "drug_fp", "cell_drug_fp"]:
            pred = fit_predict(train_meta, test_meta, train_y, feature_set, 10.0, 1024)
            summary_rows.append(
                eval_prediction(seed, f"ridge_{feature_set}", train_meta, test_meta, test_y, pred, test_de_mask)
            )

        leakage_rows.append(leakage_row(seed, train_meta, test_meta))
        if seed % 10 == 0:
            print(f"completed seed {seed}", flush=True)

    summary = pd.DataFrame(summary_rows)
    leakage = pd.DataFrame(leakage_rows)
    summary.to_csv(OUT / "all_baseline_summary.csv", index=False)
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
    leak_summary = pd.DataFrame(
        [
            {
                "seeds": leakage["seed"].nunique(),
                "mean_train_records": leakage["n_train"].mean(),
                "sd_train_records": leakage["n_train"].std(),
                "mean_test_records": leakage["n_test"].mean(),
                "sd_test_records": leakage["n_test"].std(),
                "mean_heldout_moa": leakage["n_heldout_moa"].mean(),
                "mean_test_drugs": leakage["n_test_drugs"].mean(),
                "mean_test_scaffolds": leakage["n_test_scaffolds"].mean(),
                "mean_same_drug_overlap": leakage["same_drug_in_train"].mean(),
                "mean_same_scaffold_overlap": leakage["same_scaffold_in_train"].mean(),
                "mean_same_moa_overlap": leakage["same_moa_in_train"].mean(),
            }
        ]
    )
    leak_summary.to_csv(OUT / "leakage_summary_aggregate.csv", index=False)
    print(agg.to_string(index=False))
    print(leak_summary.to_string(index=False))


if __name__ == "__main__":
    main()
