#!/usr/bin/env python3
"""Evaluate nearest-drug baselines using Morgan fingerprint Tanimoto similarity."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_matrix(de_train: Path, de_test: Path, layer: str, mask_layer: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    import anndata as ad

    parts = []
    xs = []
    masks = []
    for path, source_split in [(de_train, "de_train"), (de_test, "de_test")]:
        adata = ad.read_h5ad(path)
        obs = adata.obs.copy().reset_index(names="sample_id")
        obs["source_split"] = source_split
        parts.append(obs)
        xs.append(np.asarray(adata.layers[layer], dtype=np.float64))
        masks.append(np.asarray(adata.layers[mask_layer], dtype=bool))
    return pd.concat(parts, ignore_index=True), np.vstack(xs), np.vstack(masks)


def fp_from_smiles(smiles: str):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def tanimoto(fp, fps: list) -> np.ndarray:
    from rdkit import DataStructs

    if fp is None or not fps:
        return np.array([], dtype=float)
    return np.asarray(DataStructs.BulkTanimotoSimilarity(fp, fps), dtype=float)


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


def predict_nearest(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    train_y: np.ndarray,
    mode: str,
) -> tuple[np.ndarray, list[float]]:
    global_mean = train_y.mean(axis=0)
    fp_cache = {drug: fp_from_smiles(smiles) for drug, smiles in train_meta[["drug_name", "smiles"]].drop_duplicates().values}
    test_fp_cache = {drug: fp_from_smiles(smiles) for drug, smiles in test_meta[["drug_name", "smiles"]].drop_duplicates().values}

    preds = []
    max_sims = []
    for _, row in test_meta.iterrows():
        candidates = train_meta.copy()
        if mode == "nearest_drug_same_cell":
            same_cell = candidates["cell_context"].eq(row["cell_context"])
            if same_cell.any():
                candidates = candidates.loc[same_cell].copy()
        fps = [fp_cache.get(d) for d in candidates["drug_name"]]
        valid = [i for i, fp in enumerate(fps) if fp is not None]
        if not valid:
            preds.append(global_mean)
            max_sims.append(np.nan)
            continue
        candidates = candidates.iloc[valid].copy()
        fps = [fps[i] for i in valid]
        sims = tanimoto(test_fp_cache.get(row["drug_name"]), fps)
        if sims.size == 0:
            preds.append(global_mean)
            max_sims.append(np.nan)
            continue
        best = int(np.argmax(sims))
        train_index = candidates.index[best]
        preds.append(train_y[train_meta.index.get_loc(train_index)])
        max_sims.append(float(sims[best]))
    return np.vstack(preds), max_sims


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
    matrix_meta, y, de_mask = load_matrix(de_train, de_test, layer, mask_layer)
    splits = pd.read_csv(splits_path)
    split_cols = [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]
    meta = matrix_meta.merge(
        splits[["sample_id", "drug_name", "cell_context", "condition", "smiles", *split_cols]],
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

        for mode in ["nearest_drug_any_cell", "nearest_drug_same_cell"]:
            pred, sims = predict_nearest(train_meta, test_meta, train_y, mode)
            corr = rowwise_pearson(test_y, pred)
            err = rmse(test_y, pred)
            de_corr = masked_rowwise_pearson(test_y, pred, test_de_mask, min_de_genes)
            de_err = masked_rmse(test_y, pred, test_de_mask, min_de_genes)
            rows.append(
                {
                    "split": split_col,
                    "baseline": mode,
                    "n_train": len(train_meta),
                    "n_test": len(test_meta),
                    "layer": layer,
                    "mask_layer": mask_layer,
                    "min_de_genes": min_de_genes,
                    "mean_nearest_tanimoto": float(np.nanmean(sims)),
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
                    "baseline": mode,
                    "sample_id": sid,
                    "nearest_tanimoto": sim,
                    "rowwise_pearson": c,
                    "rmse": e,
                    "de_rowwise_pearson": dc,
                    "de_rmse": de,
                }
                for sid, sim, c, e, dc, de in zip(test_meta["sample_id"], sims, corr, err, de_corr, de_err)
            )

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "nearest_drug_baseline_summary.csv", index=False)
    pd.DataFrame(per_row).to_csv(outdir / "nearest_drug_baseline_per_row.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/nearest_drug_openproblems"))
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    args = parser.parse_args()
    evaluate(args.de_train, args.de_test, args.splits, args.outdir, args.layer, args.mask_layer, args.min_de_genes)


if __name__ == "__main__":
    main()
