#!/usr/bin/env python3
"""Evaluate top-k differential-response gene overlap and direction agreement."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def load_matrix(de_train: Path, de_test: Path, layer: str) -> tuple[pd.DataFrame, np.ndarray]:
    import anndata as ad

    parts = []
    xs = []
    for path, source_split in [(de_train, "de_train"), (de_test, "de_test")]:
        adata = ad.read_h5ad(path)
        if layer not in adata.layers:
            raise SystemExit(f"Layer '{layer}' not found in {path}.")
        obs = adata.obs.copy().reset_index(names="sample_id")
        obs["source_split"] = source_split
        parts.append(obs)
        xs.append(np.asarray(adata.layers[layer], dtype=np.float64))
    return pd.concat(parts, ignore_index=True), np.vstack(xs)


def morgan_fingerprint(smiles: str, n_bits: int = 1024, radius: int = 2) -> np.ndarray:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    arr = np.zeros((n_bits,), dtype=np.float64)
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr[list(fp.GetOnBits())] = 1.0
    return arr


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


def build_ridge_features(train_meta: pd.DataFrame, test_meta: pd.DataFrame, n_bits: int) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.preprocessing import OneHotEncoder

    encoder = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    train_cell = encoder.fit_transform(train_meta[["cell_context"]]).astype(np.float64)
    test_cell = encoder.transform(test_meta[["cell_context"]]).astype(np.float64)
    train_fp = np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in train_meta["smiles"]])
    test_fp = np.vstack([morgan_fingerprint(s, n_bits=n_bits) for s in test_meta["smiles"]])
    return np.hstack([train_cell, train_fp]), np.hstack([test_cell, test_fp])


def ridge_predict(train_meta: pd.DataFrame, test_meta: pd.DataFrame, train_y: np.ndarray, alpha: float, n_bits: int) -> np.ndarray:
    train_x, test_x = build_ridge_features(train_meta, test_meta, n_bits)
    train_design = np.hstack([np.ones((train_x.shape[0], 1)), train_x])
    test_design = np.hstack([np.ones((test_x.shape[0], 1)), test_x])
    if not np.isfinite(train_design).all() or not np.isfinite(test_design).all() or not np.isfinite(train_y).all():
        raise SystemExit("Non-finite values found in top-k ridge inputs.")
    penalty = np.eye(train_design.shape[1]) * alpha
    penalty[0, 0] = 0.0
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*encountered in matmul.*")
        coef = np.linalg.solve(train_design.T @ train_design + penalty, train_design.T @ train_y)
        pred = test_design @ coef
    if not np.isfinite(pred).all():
        raise SystemExit("Non-finite values produced by top-k ridge prediction.")
    return pred


def nearest_drug_predict(train_meta: pd.DataFrame, test_meta: pd.DataFrame, train_y: np.ndarray) -> np.ndarray:
    global_mean = train_y.mean(axis=0)
    fp_cache = {drug: fp_from_smiles(smiles) for drug, smiles in train_meta[["drug_name", "smiles"]].drop_duplicates().values}
    test_fp_cache = {drug: fp_from_smiles(smiles) for drug, smiles in test_meta[["drug_name", "smiles"]].drop_duplicates().values}
    preds = []
    for _, row in test_meta.iterrows():
        candidates = train_meta.copy()
        fps = [fp_cache.get(d) for d in candidates["drug_name"]]
        valid = [i for i, fp in enumerate(fps) if fp is not None]
        if not valid:
            preds.append(global_mean)
            continue
        candidates = candidates.iloc[valid].copy()
        fps = [fps[i] for i in valid]
        sims = tanimoto(test_fp_cache.get(row["drug_name"]), fps)
        if sims.size == 0:
            preds.append(global_mean)
            continue
        best = int(np.argmax(sims))
        train_index = candidates.index[best]
        preds.append(train_y[train_meta.index.get_loc(train_index)])
    return np.vstack(preds)


def topk_metrics(y_true: np.ndarray, y_pred: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = y_true.shape[0]
    overlap = np.full(n, np.nan)
    recall = np.full(n, np.nan)
    direction = np.full(n, np.nan)
    for i in range(n):
        true_order = np.argsort(np.abs(y_true[i]))[::-1]
        pred_order = np.argsort(np.abs(y_pred[i]))[::-1]
        true_top = set(true_order[:k])
        pred_top = set(pred_order[:k])
        shared = sorted(true_top & pred_top)
        overlap[i] = len(shared) / k
        recall[i] = len(shared) / len(true_top) if true_top else np.nan
        if shared:
            direction[i] = float(np.mean(np.sign(y_true[i, shared]) == np.sign(y_pred[i, shared])))
    return overlap, recall, direction


def evaluate(de_train: Path, de_test: Path, splits_path: Path, outdir: Path, layer: str, topks: list[int]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    matrix_meta, y = load_matrix(de_train, de_test, layer)
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

    rows = []
    per_row = []
    for split_col in split_cols:
        train_mask = meta[split_col].eq("train").to_numpy()
        test_mask = meta[split_col].eq("test").to_numpy()
        train_meta = meta.loc[train_mask].reset_index(drop=True)
        test_meta = meta.loc[test_mask].reset_index(drop=True)
        train_y = y[train_mask]
        test_y = y[test_mask]
        if len(train_meta) == 0 or len(test_meta) == 0:
            continue

        global_mean = np.tile(train_y.mean(axis=0), (len(test_meta), 1))
        ridge = ridge_predict(train_meta, test_meta, train_y, alpha=10.0, n_bits=1024)
        nearest = nearest_drug_predict(train_meta, test_meta, train_y)
        baselines = {
            "global_train_mean": global_mean,
            "ridge_cell_drug_fp": ridge,
            "nearest_drug_any_cell": nearest,
        }

        for baseline, pred in baselines.items():
            for k in topks:
                overlap, recall, direction = topk_metrics(test_y, pred, k)
                rows.append(
                    {
                        "split": split_col,
                        "baseline": baseline,
                        "top_k": k,
                        "n_train": len(train_meta),
                        "n_test": len(test_meta),
                        "mean_topk_overlap": float(np.nanmean(overlap)),
                        "sd_topk_overlap": float(np.nanstd(overlap, ddof=1)),
                        "mean_direction_agreement": float(np.nanmean(direction)),
                        "sd_direction_agreement": float(np.nanstd(direction, ddof=1)),
                        "n_direction_evaluable": int(np.isfinite(direction).sum()),
                    }
                )
                per_row.extend(
                    {
                        "split": split_col,
                        "baseline": baseline,
                        "sample_id": sid,
                        "top_k": k,
                        "topk_overlap": o,
                        "direction_agreement": d,
                    }
                    for sid, o, d in zip(test_meta["sample_id"], overlap, direction)
                )

    pd.DataFrame(rows).to_csv(outdir / "topk_de_overlap_summary.csv", index=False)
    pd.DataFrame(per_row).to_csv(outdir / "topk_de_overlap_per_row.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/topk_de_overlap"))
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--topks", default="50,100,200")
    args = parser.parse_args()
    topks = [int(x.strip()) for x in args.topks.split(",") if x.strip()]
    evaluate(args.de_train, args.de_test, args.splits, args.outdir, args.layer, topks)


if __name__ == "__main__":
    main()
