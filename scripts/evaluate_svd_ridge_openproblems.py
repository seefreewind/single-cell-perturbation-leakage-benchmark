#!/usr/bin/env python3
"""Evaluate low-rank SVD plus ridge baselines for OpenProblems."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent))
from evaluate_ridge_openproblems import (  # noqa: E402
    build_features,
    load_matrix,
    masked_rmse,
    masked_rowwise_pearson,
    rmse,
    rowwise_pearson,
)


def ridge_predict_scores(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    train_scores: np.ndarray,
    feature_set: str,
    alpha: float,
    n_bits: int,
) -> np.ndarray:
    train_x, test_x = build_features(train_meta, test_meta, feature_set, n_bits)
    train_design = np.hstack([np.ones((train_x.shape[0], 1), dtype=np.float64), train_x])
    test_design = np.hstack([np.ones((test_x.shape[0], 1), dtype=np.float64), test_x])
    if not np.isfinite(train_design).all() or not np.isfinite(test_design).all() or not np.isfinite(train_scores).all():
        raise SystemExit("Non-finite values found in SVD-ridge inputs.")
    penalty = np.eye(train_design.shape[1], dtype=np.float64) * alpha
    penalty[0, 0] = 0.0
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*encountered in matmul.*")
        coef = np.linalg.solve(train_design.T @ train_design + penalty, train_design.T @ train_scores)
        pred_scores = test_design @ coef
    if not np.isfinite(pred_scores).all():
        raise SystemExit("Non-finite values produced by SVD-ridge prediction.")
    return pred_scores


def fit_svd_ridge_predict(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    train_y: np.ndarray,
    feature_set: str,
    alpha: float,
    n_bits: int,
    n_components: int,
) -> tuple[np.ndarray, int]:
    max_components = min(n_components, train_y.shape[0] - 1, train_y.shape[1])
    if max_components < 1:
        raise SystemExit("Not enough training rows for SVD-ridge.")

    train_mean = train_y.mean(axis=0, keepdims=True)
    centered = train_y - train_mean
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    components = vt[:max_components]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*encountered in matmul.*")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in matmul")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="overflow encountered in matmul")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in matmul")
        train_scores = centered @ components.T
    if not np.isfinite(train_scores).all():
        raise SystemExit("Non-finite values produced by SVD score projection.")
    pred_scores = ridge_predict_scores(train_meta, test_meta, train_scores, feature_set, alpha, n_bits)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*encountered in matmul.*")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in matmul")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="overflow encountered in matmul")
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in matmul")
        pred = pred_scores @ components + train_mean
    if not np.isfinite(pred).all():
        raise SystemExit("Non-finite values produced after SVD reconstruction.")
    return pred, max_components


def evaluate(
    de_train: Path,
    de_test: Path,
    splits_path: Path,
    outdir: Path,
    layer: str,
    mask_layer: str,
    min_de_genes: int,
    alpha: float,
    n_bits: int,
    n_components: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    matrix_meta, y, de_mask = load_matrix(de_train, de_test, layer, mask_layer)
    splits = pd.read_csv(splits_path)
    split_cols = [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]
    meta = matrix_meta.merge(
        splits[
            [
                "sample_id",
                "drug_name",
                "cell_context",
                "condition",
                "smiles",
                *split_cols,
            ]
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    treated = meta["condition"].eq("treated").to_numpy()
    meta = meta.loc[treated].reset_index(drop=True)
    y = y[treated]
    de_mask = de_mask[treated]

    rows = []
    per_row = []
    for split_col in split_cols:
        train_mask = meta[split_col].eq("train").to_numpy()
        test_mask = meta[split_col].eq("test").to_numpy()
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_y = y[train_mask]
        test_y = y[test_mask]
        test_de_mask = de_mask[test_mask]
        if len(train_meta) == 0 or len(test_meta) == 0:
            continue

        for feature_set in ["cell", "drug_fp", "cell_drug_fp"]:
            pred, used_components = fit_svd_ridge_predict(
                train_meta,
                test_meta,
                train_y,
                feature_set,
                alpha,
                n_bits,
                n_components,
            )
            corr = rowwise_pearson(test_y, pred)
            err = rmse(test_y, pred)
            de_corr = masked_rowwise_pearson(test_y, pred, test_de_mask, min_de_genes)
            de_err = masked_rmse(test_y, pred, test_de_mask, min_de_genes)
            baseline = f"svd{used_components}_ridge_{feature_set}"
            rows.append(
                {
                    "split": split_col,
                    "baseline": baseline,
                    "n_train": len(train_meta),
                    "n_test": len(test_meta),
                    "layer": layer,
                    "mask_layer": mask_layer,
                    "min_de_genes": min_de_genes,
                    "alpha": alpha,
                    "n_bits": n_bits,
                    "requested_components": n_components,
                    "used_components": used_components,
                    "mean_rowwise_pearson": float(np.nanmean(corr)),
                    "median_rowwise_pearson": float(np.nanmedian(corr)),
                    "mean_rmse": float(np.nanmean(err)),
                    "median_rmse": float(np.nanmedian(err)),
                    "mean_de_rowwise_pearson": float(np.nanmean(de_corr)),
                    "median_de_rowwise_pearson": float(np.nanmedian(de_corr)),
                    "mean_de_rmse": float(np.nanmean(de_err)),
                    "median_de_rmse": float(np.nanmedian(de_err)),
                    "n_de_evaluable": int(np.isfinite(de_corr).sum()),
                }
            )
            per_row.extend(
                {
                    "split": split_col,
                    "baseline": baseline,
                    "sample_id": sample_id,
                    "rowwise_pearson": c,
                    "rmse": e,
                    "de_rowwise_pearson": dc,
                    "de_rmse": de,
                }
                for sample_id, c, e, dc, de in zip(test_meta["sample_id"], corr, err, de_corr, de_err)
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "svd_ridge_baseline_summary.csv", index=False)
    pd.DataFrame(per_row).to_csv(outdir / "svd_ridge_baseline_per_row.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/svd_ridge_openproblems"))
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
    parser.add_argument("--n-components", type=int, default=50)
    args = parser.parse_args()
    evaluate(
        args.de_train,
        args.de_test,
        args.splits,
        args.outdir,
        args.layer,
        args.mask_layer,
        args.min_de_genes,
        args.alpha,
        args.n_bits,
        args.n_components,
    )


if __name__ == "__main__":
    main()
