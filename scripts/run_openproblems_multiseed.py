#!/usr/bin/env python3
"""Run repeated split generation, baseline evaluation, and similarity audit for OpenProblems."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def run(cmd: list[str], quiet: bool = False) -> None:
    print(">>", " ".join(cmd))
    if quiet:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default="env/.venv/bin/python")
    parser.add_argument("--metadata", type=Path, default=Path("results/feasibility_openproblems_full/perturbation_metadata_with_scaffolds.csv"))
    parser.add_argument("--de-train", type=Path, default=Path("data/raw/openproblems_neurips2023/de_train.h5ad"))
    parser.add_argument("--de-test", type=Path, default=Path("data/raw/openproblems_neurips2023/de_test.h5ad"))
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--outdir", type=Path, default=Path("results/openproblems_multiseed"))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]

    mean_frames = []
    ridge_frames = []
    sim_frames = []
    split_frames = []

    for seed in seeds:
        seed_dir = args.outdir / f"seed_{seed:03d}"
        split_dir = seed_dir / "splits"
        mean_dir = seed_dir / "mean_baselines"
        ridge_dir = seed_dir / "ridge_baselines"
        sim_dir = seed_dir / "drug_similarity"

        run(
            [
                args.python,
                "scripts/make_split_plan.py",
                "--metadata",
                str(args.metadata),
                "--outdir",
                str(split_dir),
                "--test-fraction",
                str(args.test_fraction),
                "--seed",
                str(seed),
            ],
            quiet=args.quiet,
        )
        run(
            [
                args.python,
                "scripts/evaluate_simple_baselines_openproblems.py",
                "--de-train",
                str(args.de_train),
                "--de-test",
                str(args.de_test),
                "--splits",
                str(split_dir / "candidate_splits.csv"),
                "--outdir",
                str(mean_dir),
                "--layer",
                args.layer,
            ],
            quiet=args.quiet,
        )
        run(
            [
                args.python,
                "scripts/evaluate_ridge_openproblems.py",
                "--de-train",
                str(args.de_train),
                "--de-test",
                str(args.de_test),
                "--splits",
                str(split_dir / "candidate_splits.csv"),
                "--outdir",
                str(ridge_dir),
                "--layer",
                args.layer,
                "--alpha",
                "10",
                "--n-bits",
                "1024",
            ],
            quiet=args.quiet,
        )
        run(
            [
                args.python,
                "scripts/drug_similarity_audit.py",
                "--splits",
                str(split_dir / "candidate_splits.csv"),
                "--per-row",
                str(ridge_dir / "ridge_baseline_per_row.csv"),
                "--outdir",
                str(sim_dir),
            ],
            quiet=args.quiet,
        )

        mean = pd.read_csv(mean_dir / "simple_baseline_summary.csv")
        mean["seed"] = seed
        mean_frames.append(mean)

        ridge = pd.read_csv(ridge_dir / "ridge_baseline_summary.csv")
        ridge["seed"] = seed
        ridge_frames.append(ridge)

        sim = pd.read_csv(sim_dir / "drug_similarity_vs_error_summary.csv")
        sim["seed"] = seed
        sim_frames.append(sim)

        split_summary = pd.read_csv(split_dir / "candidate_split_summary.csv")
        split_summary["seed"] = seed
        split_frames.append(split_summary)

    all_mean = pd.concat(mean_frames, ignore_index=True)
    all_ridge = pd.concat(ridge_frames, ignore_index=True)
    all_sim = pd.concat(sim_frames, ignore_index=True)
    all_splits = pd.concat(split_frames, ignore_index=True)

    all_mean.to_csv(args.outdir / "all_simple_baseline_summary.csv", index=False)
    all_ridge.to_csv(args.outdir / "all_ridge_baseline_summary.csv", index=False)
    all_sim.to_csv(args.outdir / "all_drug_similarity_summary.csv", index=False)
    all_splits.to_csv(args.outdir / "all_split_summary.csv", index=False)

    combined = pd.concat([all_mean, all_ridge], ignore_index=True, sort=False)
    combined.to_csv(args.outdir / "all_baseline_summary.csv", index=False)

    grouped = (
        combined.groupby(["split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_train=("n_train", "mean"),
            mean_n_test=("n_test", "mean"),
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
    grouped.to_csv(args.outdir / "baseline_summary_by_seed_mean.csv", index=False)
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
