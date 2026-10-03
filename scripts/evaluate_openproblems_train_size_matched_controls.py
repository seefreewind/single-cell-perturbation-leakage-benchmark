#!/usr/bin/env python3
"""Evaluate train/test-size-matched random controls for the OP3 primary ridge.

Each control is derived solely from the random train/test pools for the same
fixed seed. The joint manifest supplies target counts and, for composition
matching, test-cell and DE-count strata only. Joint test/excluded records never
enter the random-control training set.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from evaluate_ridge_openproblems import (
    fit_predict,
    load_matrix,
    masked_rowwise_pearson,
    rmse,
    rowwise_pearson,
)


def de_counts(de_train: Path, de_test: Path, mask_layer: str) -> pd.Series:
    import anndata as ad

    items = []
    for path in (de_train, de_test):
        data = ad.read_h5ad(path, backed="r")
        mask = data.layers[mask_layer]
        count = np.asarray(mask.toarray() if hasattr(mask, "toarray") else mask, dtype=bool).sum(axis=1)
        items.append(pd.Series(count, index=data.obs_names.astype(str)))
        data.file.close()
    return pd.concat(items)


def draw_composition_matched(
    random_test: pd.DataFrame,
    joint_test: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[str]]:
    """Draw one random-test ID per joint-test row, with documented fallback."""
    combined = pd.concat([random_test[["n_de_genes"]], joint_test[["n_de_genes"]]], ignore_index=True)
    try:
        bins = pd.qcut(combined["n_de_genes"], q=4, labels=False, duplicates="drop")
        random = random_test.copy()
        target = joint_test.copy()
        random["de_quartile"] = bins.iloc[: len(random)].to_numpy()
        target["de_quartile"] = bins.iloc[len(random) :].to_numpy()
    except ValueError:
        random = random_test.assign(de_quartile=0)
        target = joint_test.assign(de_quartile=0)
    groups = random.groupby(["cell_context", "de_quartile"], dropna=False).groups
    selected, fallbacks = [], []
    for row in target.itertuples(index=False):
        candidates = np.asarray(groups.get((row.cell_context, row.de_quartile), []))
        source = "cell_and_de_quartile"
        if not len(candidates):
            candidates = random.index[random["cell_context"].eq(row.cell_context)].to_numpy()
            source = "cell_only"
        if not len(candidates):
            candidates = random.index.to_numpy()
            source = "random_pool"
        selected.append(int(rng.choice(candidates)))
        fallbacks.append(source)
    return np.asarray(selected), fallbacks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, default=ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad")
    parser.add_argument("--de-test", type=Path, default=ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad")
    parser.add_argument("--manifest-root", type=Path, default=ROOT / "results/openproblems_multiseed_100_de")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results/openproblems_train_size_matched_100")
    parser.add_argument("--seeds", default="1-100")
    parser.add_argument("--layer", default="clipped_sign_log10_pval")
    parser.add_argument("--mask-layer", default="is_de")
    parser.add_argument("--min-de-genes", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=10.0)
    parser.add_argument("--n-bits", type=int, default=1024)
    args = parser.parse_args()

    if args.seeds == "1-100":
        seeds = list(range(1, 101))
    else:
        seeds = [int(x) for x in args.seeds.split(",")]
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrix_meta, y, mask = load_matrix(args.de_train, args.de_test, args.layer, args.mask_layer)
    count_by_id = de_counts(args.de_train, args.de_test, args.mask_layer)
    rows, draws, failures = [], [], []
    for seed in seeds:
        try:
            manifest = pd.read_csv(args.manifest_root / f"seed_{seed:03d}" / "splits/candidate_splits.csv")
            meta = matrix_meta.merge(
                manifest[["sample_id", "drug_name", "cell_context", "condition", "smiles", "split_random", "split_cell_scaffold_heldout"]],
                on="sample_id", how="left", validate="one_to_one",
            )
            treated = meta["condition"].eq("treated").to_numpy()
            meta = meta.loc[treated].reset_index(drop=True)
            target_y, target_mask = y[treated], mask[treated]
            meta["n_de_genes"] = meta["sample_id"].astype(str).map(count_by_id).astype(int)
            random_train = meta[meta.split_random.eq("train")].copy()
            random_test = meta[meta.split_random.eq("test")].copy()
            joint_train = meta[meta.split_cell_scaffold_heldout.eq("train")].copy()
            joint_test = meta[meta.split_cell_scaffold_heldout.eq("test")].copy()
            if not len(random_train) or not len(random_test) or not len(joint_train) or not len(joint_test):
                raise ValueError("empty random or joint partition")
            rng = np.random.default_rng(81000 + seed)
            test_size_idx = rng.choice(random_test.index.to_numpy(), size=len(joint_test), replace=False)
            train_size_idx = rng.choice(random_train.index.to_numpy(), size=len(joint_train), replace=False)
            composition_test_idx, fallback = draw_composition_matched(random_test, joint_test, np.random.default_rng(91000 + seed))

            controls = {
                "test_size_matched_random": {
                    "train": random_train,
                    "test": random_test.loc[test_size_idx],
                    "matched_train": False,
                    "matched_composition": False,
                },
                "train_test_size_matched_random": {
                    "train": random_train.loc[train_size_idx],
                    "test": random_test.loc[test_size_idx],
                    "matched_train": True,
                    "matched_composition": False,
                },
                "composition_train_test_size_matched_random": {
                    "train": random_train.loc[train_size_idx],
                    "test": random_test.loc[composition_test_idx],
                    "matched_train": True,
                    "matched_composition": True,
                },
            }
            for position, sample_id in enumerate(random_test.loc[composition_test_idx, "sample_id"].astype(str)):
                draws.append({"seed": seed, "control": "composition_train_test_size_matched_random", "sample_id": sample_id, "fallback": fallback[position]})
            for control, parts in controls.items():
                train = parts["train"].reset_index()
                test = parts["test"].reset_index()
                train_y = target_y[train["index"].to_numpy()]
                test_y = target_y[test["index"].to_numpy()]
                test_mask = target_mask[test["index"].to_numpy()]
                pred = fit_predict(train, test, train_y, "cell_drug_fp", args.alpha, args.n_bits)
                corr = rowwise_pearson(test_y, pred)
                de_corr = masked_rowwise_pearson(test_y, pred, test_mask, args.min_de_genes)
                rows.append({
                    "seed": seed,
                    "control": control,
                    "model": "ridge_cell_drug_fp",
                    "metric": "rowwise_pearson",
                    "value": float(np.nanmean(corr)),
                    "n_train": len(train),
                    "n_test": len(test),
                    "mean_rmse": float(np.nanmean(rmse(test_y, pred))),
                    "mean_de_gene_pearson": float(np.nanmean(de_corr)),
                    "n_de_evaluable": int(np.isfinite(de_corr).sum()),
                    "matched_train_size": parts["matched_train"],
                    "matched_test_size": True,
                    "matched_cell_and_de_quartile": parts["matched_composition"],
                    "source_random_train_only": True,
                    "source_random_test_only": True,
                })
        except Exception as exc:
            failures.append({"seed": seed, "status": "failure", "error": f"{type(exc).__name__}: {exc}"})
            print(f"seed {seed}: {failures[-1]['error']}", flush=True)
        else:
            print(f"seed {seed}: completed", flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "train_size_matched_random_metrics_by_seed.csv", index=False)
    pd.DataFrame(draws).to_csv(args.outdir / "train_size_matched_random_draws.csv", index=False)
    pd.DataFrame(failures).to_csv(args.outdir / "train_size_matched_random_failures.csv", index=False)
    if not result.empty:
        summary = result.groupby(["control", "model", "metric"], as_index=False).agg(
            n_seeds=("seed", "nunique"),
            mean=("value", "mean"),
            sd=("value", "std"),
            median=("value", "median"),
            iqr=("value", lambda x: x.quantile(0.75) - x.quantile(0.25)),
            mean_n_train=("n_train", "mean"),
            mean_n_test=("n_test", "mean"),
        )
        summary.to_csv(args.outdir / "train_size_matched_random_summary.csv", index=False)


if __name__ == "__main__":
    main()
