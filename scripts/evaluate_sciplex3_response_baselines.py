#!/usr/bin/env python3
"""Evaluate simple and ridge baselines on a Sci-Plex response AnnData."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from evaluate_ridge_openproblems import morgan_fingerprint, rmse, rowwise_pearson


def mean_predict(train_meta: pd.DataFrame, train_y: np.ndarray, test_meta: pd.DataFrame, key: str | None) -> np.ndarray:
    if key is None:
        return np.repeat(train_y.mean(axis=0, keepdims=True), len(test_meta), axis=0)
    means = {k: train_y[idx].mean(axis=0) for k, idx in train_meta.groupby(key, observed=True).indices.items()}
    global_mean = train_y.mean(axis=0)
    return np.vstack([means.get(v, global_mean) for v in test_meta[key]])


def build_features(train_meta: pd.DataFrame, test_meta: pd.DataFrame, feature_set: str, n_bits: int) -> tuple[np.ndarray, np.ndarray]:
    train_blocks = []
    test_blocks = []
    if "cell" in feature_set:
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        train_blocks.append(enc.fit_transform(train_meta[["cell_context"]]).astype(np.float64))
        test_blocks.append(enc.transform(test_meta[["cell_context"]]).astype(np.float64))
    if "dose" in feature_set:
        scaler = StandardScaler()
        train_dose = np.log10(train_meta[["dose_value"]].astype(float).to_numpy())
        test_dose = np.log10(test_meta[["dose_value"]].astype(float).to_numpy())
        train_blocks.append(scaler.fit_transform(train_dose).astype(np.float64))
        test_blocks.append(scaler.transform(test_dose).astype(np.float64))
    if "drug_fp" in feature_set:
        train_blocks.append(np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in train_meta["smiles"]]))
        test_blocks.append(np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in test_meta["smiles"]]))
    return np.hstack(train_blocks), np.hstack(test_blocks)


def ridge_predict(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    train_y: np.ndarray,
    feature_set: str,
    alpha: float,
    n_bits: int,
) -> np.ndarray:
    train_x, test_x = build_features(train_meta, test_meta, feature_set, n_bits)
    train_x_aug = np.column_stack([np.ones(train_x.shape[0]), train_x])
    test_x_aug = np.column_stack([np.ones(test_x.shape[0]), test_x])
    penalty = np.eye(train_x_aug.shape[1], dtype=np.float64) * alpha
    penalty[0, 0] = 0.0
    coef = np.linalg.solve(train_x_aug.T @ train_x_aug + penalty, train_x_aug.T @ train_y)
    return test_x_aug @ coef


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/sciplex3_baselines"))
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    meta = adata.obs.copy().reset_index(names="sample_id")
    y = np.asarray(adata.X, dtype=np.float64)
    splits = pd.read_csv(args.splits)
    split_cols = [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]
    meta = meta.merge(
        splits[["sample_id", "drug_name", "cell_context", "dose_value", "smiles", "scaffold", *split_cols]],
        on="sample_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_split"),
    )
    rows = []
    per_row = []
    mean_specs = [
        ("global_train_mean", None),
        ("cell_context_mean", "cell_context"),
        ("drug_mean", "drug_name"),
        ("cell_dose_mean", "cell_dose_key"),
    ]
    meta["cell_dose_key"] = meta["cell_context"].astype(str) + "||" + meta["dose_value"].astype(str)
    ridge_specs = [
        ("ridge_cell_dose", "cell+dose"),
        ("ridge_drug_fp", "drug_fp"),
        ("ridge_cell_dose_drug_fp", "cell+dose+drug_fp"),
    ]

    for split_col in split_cols:
        train_mask = meta[split_col].eq("train").to_numpy()
        test_mask = meta[split_col].eq("test").to_numpy()
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_y = y[train_mask]
        test_y = y[test_mask]
        if len(train_meta) == 0 or len(test_meta) == 0:
            continue
        predictions = []
        for baseline, key in mean_specs:
            predictions.append((baseline, mean_predict(train_meta, train_y, test_meta, key)))
        for baseline, feature_set in ridge_specs:
            predictions.append((baseline, ridge_predict(train_meta, test_meta, train_y, feature_set, args.alpha, args.n_bits)))
        for baseline, pred in predictions:
            corr = rowwise_pearson(test_y, pred)
            err = rmse(test_y, pred)
            rows.append(
                {
                    "split": split_col,
                    "baseline": baseline,
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "mean_rowwise_pearson": float(np.nanmean(corr)),
                    "median_rowwise_pearson": float(np.nanmedian(corr)),
                    "mean_rmse": float(np.nanmean(err)),
                    "median_rmse": float(np.nanmedian(err)),
                    "alpha": args.alpha if baseline.startswith("ridge") else np.nan,
                    "n_bits": args.n_bits if "drug_fp" in baseline else np.nan,
                }
            )
            per_row.extend(
                {
                    "split": split_col,
                    "baseline": baseline,
                    "sample_id": sample_id,
                    "rowwise_pearson": c,
                    "rmse": e,
                }
                for sample_id, c, e in zip(test_meta["sample_id"], corr, err)
            )
    summary = pd.DataFrame(rows)
    summary.to_csv(args.outdir / "baseline_summary.csv", index=False)
    pd.DataFrame(per_row).to_csv(args.outdir / "baseline_per_row.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
