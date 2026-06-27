#!/usr/bin/env python3
"""Evaluate ridge baselines with drug fingerprints and cell-context features."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder


def load_matrix(de_train: Path, de_test: Path, layer: str, mask_layer: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    import anndata as ad

    parts = []
    xs = []
    masks = []
    for path, source_split in [(de_train, "de_train"), (de_test, "de_test")]:
        adata = ad.read_h5ad(path)
        if layer not in adata.layers:
            raise SystemExit(f"Layer '{layer}' not found in {path}.")
        if mask_layer not in adata.layers:
            raise SystemExit(f"Mask layer '{mask_layer}' not found in {path}.")
        obs = adata.obs.copy().reset_index(names="sample_id")
        obs["source_split"] = source_split
        parts.append(obs)
        xs.append(np.asarray(adata.layers[layer], dtype=np.float64))
        masks.append(np.asarray(adata.layers[mask_layer], dtype=bool))
    return pd.concat(parts, ignore_index=True), np.vstack(xs), np.vstack(masks)


def morgan_fingerprint(smiles: str, n_bits: int = 1024, radius: int = 2) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    arr = np.zeros((n_bits,), dtype=np.float64)
    if not isinstance(smiles, str) or not smiles:
        return arr
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    on_bits = list(fp.GetOnBits())
    arr[on_bits] = 1.0
    return arr


def rowwise_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    y_true = y_true - y_true.mean(axis=1, keepdims=True)
    y_pred = y_pred - y_pred.mean(axis=1, keepdims=True)
    numerator = (y_true * y_pred).sum(axis=1)
    denominator = np.sqrt((y_true**2).sum(axis=1) * (y_pred**2).sum(axis=1))
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator != 0)


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


def build_features(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    feature_set: str,
    n_bits: int,
) -> tuple[np.ndarray, np.ndarray]:
    train_blocks = []
    test_blocks = []
    if feature_set in {"cell", "cell_drug_fp"}:
        encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        train_cell = encoder.fit_transform(train_meta[["cell_context"]]).astype(np.float64)
        test_cell = encoder.transform(test_meta[["cell_context"]]).astype(np.float64)
        train_blocks.append(train_cell)
        test_blocks.append(test_cell)
    if feature_set in {"drug_fp", "cell_drug_fp"}:
        train_fp = np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in train_meta["smiles"]])
        test_fp = np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in test_meta["smiles"]])
        train_blocks.append(train_fp)
        test_blocks.append(test_fp)
    if not train_blocks:
        raise ValueError(feature_set)
    return np.hstack(train_blocks).astype(np.float64), np.hstack(test_blocks).astype(np.float64)


def fit_predict(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    train_y: np.ndarray,
    feature_set: str,
    alpha: float,
    n_bits: int,
) -> np.ndarray:
    train_x, test_x = build_features(train_meta, test_meta, feature_set, n_bits)
    train_design = np.hstack([np.ones((train_x.shape[0], 1), dtype=np.float64), train_x])
    test_design = np.hstack([np.ones((test_x.shape[0], 1), dtype=np.float64), test_x])
    if not np.isfinite(train_design).all() or not np.isfinite(test_design).all() or not np.isfinite(train_y).all():
        raise SystemExit("Non-finite values found in ridge inputs.")
    penalty = np.eye(train_design.shape[1], dtype=np.float64) * alpha
    penalty[0, 0] = 0.0
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*encountered in matmul.*")
        coef = np.linalg.solve(train_design.T @ train_design + penalty, train_design.T @ train_y)
        pred = test_design @ coef
    if not np.isfinite(pred).all():
        raise SystemExit("Non-finite values produced by ridge prediction.")
    return pred


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
            pred = fit_predict(train_meta, test_meta, train_y, feature_set, alpha, n_bits)
            corr = rowwise_pearson(test_y, pred)
            err = rmse(test_y, pred)
            de_corr = masked_rowwise_pearson(test_y, pred, test_de_mask, min_de_genes)
            de_err = masked_rmse(test_y, pred, test_de_mask, min_de_genes)
            rows.append(
                {
                    "split": split_col,
                    "baseline": f"ridge_{feature_set}",
                    "n_train": len(train_meta),
                    "n_test": len(test_meta),
                    "layer": layer,
                    "mask_layer": mask_layer,
                    "min_de_genes": min_de_genes,
                    "alpha": alpha,
                    "n_bits": n_bits,
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
                    "baseline": f"ridge_{feature_set}",
                    "sample_id": sample_id,
                    "rowwise_pearson": c,
                    "rmse": e,
                    "de_rowwise_pearson": dc,
                    "de_rmse": de,
                }
                for sample_id, c, e, dc, de in zip(test_meta["sample_id"], corr, err, de_corr, de_err)
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "ridge_baseline_summary.csv", index=False)
    pd.DataFrame(per_row).to_csv(outdir / "ridge_baseline_per_row.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/ridge_openproblems"))
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
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
    )


if __name__ == "__main__":
    main()
