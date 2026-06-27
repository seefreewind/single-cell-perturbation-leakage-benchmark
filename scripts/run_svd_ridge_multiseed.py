#!/usr/bin/env python3
"""Run SVD-ridge baselines over existing OpenProblems split seeds."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def run(cmd: list[str], quiet: bool) -> None:
    print(">>", " ".join(cmd))
    if quiet:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default="env/.venv/bin/python")
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--de-train", type=Path, default=Path("data/raw/openproblems_neurips2023/de_train.h5ad"))
    parser.add_argument("--de-test", type=Path, default=Path("data/raw/openproblems_neurips2023/de_test.h5ad"))
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--outdir", type=Path, default=Path("results/svd_ridge_multiseed_10"))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    frames = []

    for seed in seeds:
        seed_dir = args.outdir / f"seed_{seed:03d}"
        split_path = args.split_root / f"seed_{seed:03d}" / "splits" / "candidate_splits.csv"
        run(
            [
                args.python,
                "scripts/evaluate_svd_ridge_openproblems.py",
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
                "--mask-layer",
                args.mask_layer,
                "--min-de-genes",
                str(args.min_de_genes),
                "--alpha",
                str(args.alpha),
                "--n-bits",
                str(args.n_bits),
                "--n-components",
                str(args.n_components),
            ],
            quiet=args.quiet,
        )
        summary = pd.read_csv(seed_dir / "svd_ridge_baseline_summary.csv")
        summary["seed"] = seed
        frames.append(summary)

    all_summary = pd.concat(frames, ignore_index=True)
    all_summary.to_csv(args.outdir / "all_svd_ridge_baseline_summary.csv", index=False)
    grouped = (
        all_summary.groupby(["split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_train=("n_train", "mean"),
            mean_n_test=("n_test", "mean"),
            mean_used_components=("used_components", "mean"),
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
    grouped.to_csv(args.outdir / "svd_ridge_summary_by_seed_mean.csv", index=False)
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
