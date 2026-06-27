#!/usr/bin/env python3
"""Audit whether test-drug proximity to training drugs explains prediction performance."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def fp_from_smiles(smiles: str):
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def max_tanimoto(test_fp, train_fps: list) -> float:
    from rdkit import DataStructs

    if test_fp is None or not train_fps:
        return np.nan
    sims = DataStructs.BulkTanimotoSimilarity(test_fp, train_fps)
    return float(max(sims)) if sims else np.nan


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    data = pd.concat([x, y], axis=1).dropna()
    if data.shape[0] < 3:
        return np.nan
    if data.iloc[:, 0].nunique() < 2 or data.iloc[:, 1].nunique() < 2:
        return np.nan
    return float(data.iloc[:, 0].corr(data.iloc[:, 1]))


def audit(splits_path: Path, per_row_path: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    splits = pd.read_csv(splits_path)
    per_row = pd.read_csv(per_row_path)

    split_cols = [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]

    drug_table = (
        splits[["drug_name", "smiles", "scaffold"]]
        .drop_duplicates("drug_name")
        .assign(fp=lambda d: d["smiles"].map(fp_from_smiles))
    )
    fp_map = dict(zip(drug_table["drug_name"], drug_table["fp"]))

    rows = []
    for split_col in split_cols:
        train_drugs = sorted(splits.loc[splits[split_col].eq("train"), "drug_name"].dropna().unique())
        train_fps = [fp_map[d] for d in train_drugs if fp_map.get(d) is not None]
        test_meta = splits.loc[splits[split_col].eq("test"), ["sample_id", "drug_name", "cell_context", "scaffold"]].copy()
        for _, row in test_meta.iterrows():
            rows.append(
                {
                    "split": split_col,
                    "sample_id": row["sample_id"],
                    "drug_name": row["drug_name"],
                    "cell_context": row["cell_context"],
                    "scaffold": row["scaffold"],
                    "max_train_tanimoto": max_tanimoto(fp_map.get(row["drug_name"]), train_fps),
                    "same_scaffold_in_train": bool(
                        splits.loc[
                            splits[split_col].eq("train") & splits["scaffold"].eq(row["scaffold"]),
                            "drug_name",
                        ].nunique()
                    ),
                }
            )

    similarity = pd.DataFrame(rows)
    merged = per_row.merge(similarity, on=["split", "sample_id"], how="left")
    merged.to_csv(outdir / "drug_similarity_vs_error_per_row.csv", index=False)

    summary = (
        merged.groupby(["split", "baseline"])
        .agg(
            n=("sample_id", "size"),
            mean_max_train_tanimoto=("max_train_tanimoto", "mean"),
            corr_tanimoto_pearson=("max_train_tanimoto", lambda s: safe_corr(s, merged.loc[s.index, "rowwise_pearson"])),
            corr_tanimoto_rmse=("max_train_tanimoto", lambda s: safe_corr(s, merged.loc[s.index, "rmse"])),
        )
        .reset_index()
    )
    summary.to_csv(outdir / "drug_similarity_vs_error_summary.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--per-row", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/drug_similarity_audit"))
    args = parser.parse_args()
    audit(args.splits, args.per_row, args.outdir)


if __name__ == "__main__":
    main()
