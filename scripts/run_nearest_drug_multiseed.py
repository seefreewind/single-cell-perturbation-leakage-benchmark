#!/usr/bin/env python3
"""Run nearest-drug baseline over existing split seeds."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default="env/.venv/bin/python")
    parser.add_argument("--split-root", type=Path, default=Path("results/openproblems_multiseed_10"))
    parser.add_argument("--de-train", type=Path, default=Path("data/raw/openproblems_neurips2023/de_train.h5ad"))
    parser.add_argument("--de-test", type=Path, default=Path("data/raw/openproblems_neurips2023/de_test.h5ad"))
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--outdir", type=Path, default=Path("results/nearest_drug_multiseed_10"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    frames = []
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    for seed in seeds:
        seed_dir = args.outdir / f"seed_{seed:03d}"
        split_path = args.split_root / f"seed_{seed:03d}" / "splits" / "candidate_splits.csv"
        cmd = [
            args.python,
            "scripts/evaluate_nearest_drug_openproblems.py",
            "--de-train",
            str(args.de_train),
            "--de-test",
            str(args.de_test),
            "--splits",
            str(split_path),
            "--outdir",
            str(seed_dir),
            "--layer",
            args.layer,
        ]
        print(">>", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        summary = pd.read_csv(seed_dir / "nearest_drug_baseline_summary.csv")
        summary["seed"] = seed
        frames.append(summary)

    all_summary = pd.concat(frames, ignore_index=True)
    all_summary.to_csv(args.outdir / "all_nearest_drug_baseline_summary.csv", index=False)
    grouped = (
        all_summary.groupby(["split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            mean_nearest_tanimoto=("mean_nearest_tanimoto", "mean"),
            sd_nearest_tanimoto=("mean_nearest_tanimoto", "std"),
            mean_rowwise_pearson=("mean_rowwise_pearson", "mean"),
            sd_rowwise_pearson=("mean_rowwise_pearson", "std"),
            mean_rmse=("mean_rmse", "mean"),
            sd_rmse=("mean_rmse", "std"),
            mean_de_rowwise_pearson=("mean_de_rowwise_pearson", "mean"),
            sd_de_rowwise_pearson=("mean_de_rowwise_pearson", "std"),
            mean_de_rmse=("mean_de_rmse", "mean"),
            sd_de_rmse=("mean_de_rmse", "std"),
            mean_n_de_evaluable=("n_de_evaluable", "mean"),
        )
        .reset_index()
    )
    grouped.to_csv(args.outdir / "nearest_drug_summary_by_seed_mean.csv", index=False)
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
