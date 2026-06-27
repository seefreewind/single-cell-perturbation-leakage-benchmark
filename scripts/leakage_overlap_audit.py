#!/usr/bin/env python3
"""Audit train-test overlap leakage for split tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SPLIT_COLS = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]


def audit_one(splits_path: Path) -> pd.DataFrame:
    df = pd.read_csv(splits_path)
    rows = []
    required = {"drug_name", "cell_context", "scaffold", *SPLIT_COLS}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns in {splits_path}: {sorted(missing)}")
    df["drug_cell_pair"] = df["drug_name"].astype(str) + "||" + df["cell_context"].astype(str)

    for split_col in SPLIT_COLS:
        train = df[df[split_col].eq("train")].copy()
        test = df[df[split_col].eq("test")].copy()
        if test.empty:
            continue
        train_drugs = set(train["drug_name"].dropna())
        train_scaffolds = set(train["scaffold"].dropna())
        train_cells = set(train["cell_context"].dropna())
        train_pairs = set(train["drug_cell_pair"].dropna())
        test_drugs = set(test["drug_name"].dropna())
        test_scaffolds = set(test["scaffold"].dropna())
        test_cells = set(test["cell_context"].dropna())
        test_pairs = set(test["drug_cell_pair"].dropna())

        rows.append(
            {
                "split": split_col,
                "n_train_records": len(train),
                "n_test_records": len(test),
                "n_excluded_records": int(df[split_col].eq("excluded").sum()),
                "test_drug_same_in_train_frac": float(test["drug_name"].isin(train_drugs).mean()),
                "test_scaffold_same_in_train_frac": float(test["scaffold"].isin(train_scaffolds).mean()),
                "test_cell_same_in_train_frac": float(test["cell_context"].isin(train_cells).mean()),
                "test_drug_cell_pair_same_in_train_frac": float(test["drug_cell_pair"].isin(train_pairs).mean()),
                "unique_test_drug_overlap_frac": len(test_drugs & train_drugs) / len(test_drugs) if test_drugs else np.nan,
                "unique_test_scaffold_overlap_frac": len(test_scaffolds & train_scaffolds) / len(test_scaffolds) if test_scaffolds else np.nan,
                "unique_test_cell_overlap_frac": len(test_cells & train_cells) / len(test_cells) if test_cells else np.nan,
                "unique_test_pair_overlap_frac": len(test_pairs & train_pairs) / len(test_pairs) if test_pairs else np.nan,
                "n_unique_test_drugs": len(test_drugs),
                "n_unique_test_scaffolds": len(test_scaffolds),
                "n_unique_test_cells": len(test_cells),
                "n_unique_test_pairs": len(test_pairs),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    parser.add_argument("--outdir", type=Path, default=Path("results/leakage_overlap_audit"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]

    frames = []
    for seed in seeds:
        split_path = args.split_root / f"seed_{seed:03d}" / "splits" / "candidate_splits.csv"
        frame = audit_one(split_path)
        frame["seed"] = seed
        frames.append(frame)

    all_rows = pd.concat(frames, ignore_index=True)
    all_rows.to_csv(args.outdir / "all_leakage_overlap_by_seed.csv", index=False)
    metric_cols = [c for c in all_rows.columns if c.endswith("_frac") or c.startswith("n_")]
    grouped = all_rows.groupby("split")[metric_cols].agg(["mean", "std", "min", "max"])
    grouped.columns = ["_".join(c).strip("_") for c in grouped.columns]
    grouped = grouped.reset_index()
    grouped.to_csv(args.outdir / "leakage_overlap_summary.csv", index=False)
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
