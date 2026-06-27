#!/usr/bin/env python3
"""Run ridge hyperparameter sensitivity on existing OpenProblems split seeds."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


def parse_csv_numbers(text: str, cast):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def run(cmd: list[str], quiet: bool) -> None:
    print(">>", " ".join(cmd))
    if quiet:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default="env/.venv/bin/python")
    parser.add_argument("--de-train", type=Path, default=Path("data/raw/openproblems_neurips2023/de_train.h5ad"))
    parser.add_argument("--de-test", type=Path, default=Path("data/raw/openproblems_neurips2023/de_test.h5ad"))
    parser.add_argument("--split-root", type=Path, default=Path("results/openproblems_multiseed_100_de"))
    parser.add_argument("--outdir", type=Path, default=Path("results/ridge_sensitivity_30"))
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30")
    parser.add_argument("--alphas", default="0.1,1,10,100,1000")
    parser.add_argument("--n-bits", default="1024")
    parser.add_argument("--extra-n-bits", default="512,2048")
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    seeds = parse_csv_numbers(args.seeds, int)
    alphas = parse_csv_numbers(args.alphas, float)
    n_bits_main = parse_csv_numbers(args.n_bits, int)
    extra_bits = parse_csv_numbers(args.extra_n_bits, int)

    configs = [(alpha, n_bits) for alpha in alphas for n_bits in n_bits_main]
    configs.extend((10.0, n_bits) for n_bits in extra_bits)
    configs = sorted(set(configs), key=lambda x: (x[1], x[0]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    frames = []
    for alpha, n_bits in configs:
        config_dir = args.outdir / f"alpha_{alpha:g}_bits_{n_bits}"
        for seed in seeds:
            split_path = args.split_root / f"seed_{seed:03d}" / "splits" / "candidate_splits.csv"
            if not split_path.exists():
                raise SystemExit(f"Missing split file: {split_path}")
            seed_dir = config_dir / f"seed_{seed:03d}"
            summary_path = seed_dir / "ridge_baseline_summary.csv"
            if not summary_path.exists():
                run(
                    [
                        args.python,
                        "scripts/evaluate_ridge_openproblems.py",
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
                        "--alpha",
                        str(alpha),
                        "--n-bits",
                        str(n_bits),
                    ],
                    args.quiet,
                )
            frame = pd.read_csv(summary_path)
            frame["seed"] = seed
            frame["sensitivity_alpha"] = alpha
            frame["sensitivity_n_bits"] = n_bits
            frames.append(frame)

    all_rows = pd.concat(frames, ignore_index=True)
    all_rows.to_csv(args.outdir / "all_ridge_sensitivity_summary.csv", index=False)
    grouped = (
        all_rows.groupby(["sensitivity_alpha", "sensitivity_n_bits", "split", "baseline"])
        .agg(
            n_seeds=("seed", "nunique"),
            mean_n_test=("n_test", "mean"),
            mean_rowwise_pearson=("mean_rowwise_pearson", "mean"),
            sd_rowwise_pearson=("mean_rowwise_pearson", "std"),
            mean_de_rowwise_pearson=("mean_de_rowwise_pearson", "mean"),
            sd_de_rowwise_pearson=("mean_de_rowwise_pearson", "std"),
            mean_de_rmse=("mean_de_rmse", "mean"),
            sd_de_rmse=("mean_de_rmse", "std"),
        )
        .reset_index()
    )
    grouped.to_csv(args.outdir / "ridge_sensitivity_by_config.csv", index=False)

    primary = grouped[grouped["baseline"].eq("ridge_cell_drug_fp")].copy()
    primary.to_csv(args.outdir / "ridge_cell_drug_fp_sensitivity.csv", index=False)
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
