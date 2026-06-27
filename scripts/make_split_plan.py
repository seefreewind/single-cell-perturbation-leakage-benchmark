#!/usr/bin/env python3
"""Create candidate split tables for random, cell-held-out, scaffold-held-out, and joint held-out tasks."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def assign_by_values(values: pd.Series, test_fraction: float, rng: np.random.Generator) -> pd.Series:
    unique_values = pd.Series(values.dropna().unique())
    n_test = max(1, int(round(len(unique_values) * test_fraction)))
    test_values = set(rng.choice(unique_values.to_numpy(), size=n_test, replace=False))
    return values.isin(test_values).map({True: "test", False: "train"})


def make_splits(metadata_path: Path, outdir: Path, test_fraction: float, seed: int) -> None:
    df = pd.read_csv(metadata_path)
    required = {"drug_id", "cell_context", "scaffold"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    split_df = df.copy()
    split_df["split_random"] = np.where(rng.random(len(split_df)) < test_fraction, "test", "train")
    split_df["split_cell_heldout"] = assign_by_values(split_df["cell_context"], test_fraction, rng)
    split_df["split_scaffold_heldout"] = assign_by_values(split_df["scaffold"], test_fraction, rng)

    heldout_cells = set(split_df.loc[split_df["split_cell_heldout"] == "test", "cell_context"])
    heldout_scaffolds = set(split_df.loc[split_df["split_scaffold_heldout"] == "test", "scaffold"])
    joint_test = split_df["cell_context"].isin(heldout_cells) & split_df["scaffold"].isin(heldout_scaffolds)
    joint_train = ~split_df["cell_context"].isin(heldout_cells) & ~split_df["scaffold"].isin(heldout_scaffolds)
    split_df["split_cell_scaffold_heldout"] = "excluded"
    split_df.loc[joint_train, "split_cell_scaffold_heldout"] = "train"
    split_df.loc[joint_test, "split_cell_scaffold_heldout"] = "test"

    split_df.to_csv(outdir / "candidate_splits.csv", index=False)

    rows = []
    for col in [
        "split_random",
        "split_cell_heldout",
        "split_scaffold_heldout",
        "split_cell_scaffold_heldout",
    ]:
        counts = split_df[col].value_counts().to_dict()
        rows.append(
            {
                "split": col,
                "train_records": counts.get("train", 0),
                "test_records": counts.get("test", 0),
                "excluded_records": counts.get("excluded", 0),
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "candidate_split_summary.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("metadata/splits"))
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    make_splits(args.metadata, args.outdir, args.test_fraction, args.seed)


if __name__ == "__main__":
    main()

