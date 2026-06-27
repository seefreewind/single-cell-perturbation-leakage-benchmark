#!/usr/bin/env python3
"""Evaluate simple mean baselines on OpenProblems DGE matrices and candidate splits."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def dense_x(adata, layer: str):
    if layer:
        if layer not in adata.layers:
            raise SystemExit(f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")
        x = adata.layers[layer]
    else:
        x = adata.X
    if x is None:
        raise SystemExit("AnnData.X is empty. Specify a layer such as --layer clipped_sign_log10_pval.")
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def rowwise_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    y_true = y_true - y_true.mean(axis=1, keepdims=True)
    y_pred = y_pred - y_pred.mean(axis=1, keepdims=True)
    numerator = (y_true * y_pred).sum(axis=1)
    denominator = np.sqrt((y_true**2).sum(axis=1) * (y_pred**2).sum(axis=1))
    out = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator != 0)
    return out


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    return np.sqrt(((y_true - y_pred) ** 2).mean(axis=1))


def masked_rowwise_pearson(y_true: np.ndarray, y_pred: np.ndarray, mask: np.ndarray, min_genes: int) -> np.ndarray:
    out = np.full(y_true.shape[0], np.nan, dtype=np.float64)
    for i in range(y_true.shape[0]):
        keep = mask[i].astype(bool)
        if keep.sum() >= min_genes:
            out[i] = rowwise_pearson(y_true[i : i + 1, keep], y_pred[i : i + 1, keep])[0]
    return out


def masked_rmse(y_true: np.ndarray, y_pred: np.ndarray, mask: np.ndarray, min_genes: int) -> np.ndarray:
    out = np.full(y_true.shape[0], np.nan, dtype=np.float64)
    for i in range(y_true.shape[0]):
        keep = mask[i].astype(bool)
        if keep.sum() >= min_genes:
            out[i] = rmse(y_true[i : i + 1, keep], y_pred[i : i + 1, keep])[0]
    return out


def grouped_mean(train_meta: pd.DataFrame, train_x: np.ndarray, key: str) -> dict[str, np.ndarray]:
    means = {}
    for value, idx in train_meta.groupby(key).indices.items():
        means[value] = train_x[list(idx)].mean(axis=0)
    return means


def predict_grouped(test_meta: pd.DataFrame, means: dict[str, np.ndarray], key: str, fallback: np.ndarray) -> np.ndarray:
    rows = []
    for value in test_meta[key]:
        rows.append(means.get(value, fallback))
    return np.vstack(rows)


def load_combined(de_train: Path, de_test: Path, layer: str, mask_layer: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    import anndata as ad

    parts = []
    xs = []
    masks = []
    for path, source_split in [(de_train, "de_train"), (de_test, "de_test")]:
        adata = ad.read_h5ad(path)
        if mask_layer not in adata.layers:
            raise SystemExit(f"Mask layer '{mask_layer}' not found. Available layers: {list(adata.layers.keys())}")
        obs = adata.obs.copy().reset_index(names="sample_id")
        obs["source_split"] = source_split
        parts.append(obs)
        xs.append(dense_x(adata, layer))
        masks.append(np.asarray(adata.layers[mask_layer], dtype=bool))
    meta = pd.concat(parts, ignore_index=True)
    x = np.vstack(xs)
    mask = np.vstack(masks)
    return meta, x, mask


def evaluate(
    de_train: Path,
    de_test: Path,
    splits_path: Path,
    outdir: Path,
    layer: str,
    mask_layer: str,
    min_de_genes: int,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    matrix_meta, x, de_mask = load_combined(de_train, de_test, layer, mask_layer)
    splits = pd.read_csv(splits_path)

    meta = matrix_meta.merge(
        splits[
            [
                "sample_id",
                "drug_name",
                "cell_context",
                "condition",
                "split_random",
                "split_cell_heldout",
                "split_scaffold_heldout",
                "split_cell_scaffold_heldout",
            ]
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    if meta[["split_random", "split_cell_heldout", "split_scaffold_heldout"]].isna().any().any():
        missing = meta.loc[meta["split_random"].isna(), "sample_id"].head(10).tolist()
        raise SystemExit(f"Split metadata did not match all matrix rows. Examples: {missing}")

    # Keep the main perturbation task focused on treated rows.
    treated_mask = meta["condition"].eq("treated").to_numpy()
    meta = meta.loc[treated_mask].reset_index(drop=True)
    x = x[treated_mask]
    de_mask = de_mask[treated_mask]

    split_cols = [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]

    rows = []
    prediction_rows = []
    for split_col in split_cols:
        train_mask = meta[split_col].eq("train").to_numpy()
        test_mask = meta[split_col].eq("test").to_numpy()
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_x = x[train_mask]
        test_x = x[test_mask]
        test_de_mask = de_mask[test_mask]
        if len(train_meta) == 0 or len(test_meta) == 0:
            continue

        global_mean = train_x.mean(axis=0)
        baselines = {
            "global_train_mean": np.tile(global_mean, (len(test_meta), 1)),
            "cell_context_mean": predict_grouped(
                test_meta, grouped_mean(train_meta, train_x, "cell_context"), "cell_context", global_mean
            ),
            "drug_mean": predict_grouped(test_meta, grouped_mean(train_meta, train_x, "drug_name"), "drug_name", global_mean),
        }

        for baseline, pred in baselines.items():
            corr = rowwise_pearson(test_x, pred)
            err = rmse(test_x, pred)
            de_corr = masked_rowwise_pearson(test_x, pred, test_de_mask, min_de_genes)
            de_err = masked_rmse(test_x, pred, test_de_mask, min_de_genes)
            rows.append(
                {
                    "split": split_col,
                    "baseline": baseline,
                    "n_train": len(train_meta),
                    "n_test": len(test_meta),
                    "layer": layer,
                    "mask_layer": mask_layer,
                    "min_de_genes": min_de_genes,
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
            prediction_rows.extend(
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

    pd.DataFrame(rows).to_csv(outdir / "simple_baseline_summary.csv", index=False)
    pd.DataFrame(prediction_rows).to_csv(outdir / "simple_baseline_per_row.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/baselines_openproblems"))
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    args = parser.parse_args()
    evaluate(args.de_train, args.de_test, args.splits, args.outdir, args.layer, args.mask_layer, args.min_de_genes)


if __name__ == "__main__":
    main()
