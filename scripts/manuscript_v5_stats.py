#!/usr/bin/env python3
"""Additional manuscript v5 robustness and leakage audits."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/manuscript_v5_stats"


SPLIT_COLS = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]


def ci(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    return tuple(np.quantile(values, [0.025, 0.975]))


def openproblems_de_counts() -> pd.DataFrame:
    import anndata as ad

    rows = []
    for path in [
        ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad",
        ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad",
    ]:
        adata = ad.read_h5ad(path)
        obs = adata.obs.copy().reset_index(names="sample_id")
        obs["n_de_genes"] = np.asarray(adata.layers["is_de"], dtype=bool).sum(axis=1)
        rows.append(obs[["sample_id", "n_de_genes"]])
    return pd.concat(rows, ignore_index=True)


def matched_random_controls(rng_seed: int = 1729, n_resamples: int = 50) -> None:
    rng = np.random.default_rng(rng_seed)
    de_counts = openproblems_de_counts()
    records = []
    for seed_dir in sorted((ROOT / "results/openproblems_multiseed_100_de").glob("seed_*")):
        seed = int(seed_dir.name.split("_")[1])
        splits = pd.read_csv(seed_dir / "splits/candidate_splits.csv")
        splits = splits.loc[splits["condition"].eq("treated")].merge(de_counts, on="sample_id", how="left")
        splits["de_bin"] = pd.qcut(
            splits["n_de_genes"].rank(method="first"),
            q=4,
            labels=False,
            duplicates="drop",
        )
        joint_meta = splits.loc[splits["split_cell_scaffold_heldout"].eq("test")]
        random_meta = splits.loc[splits["split_random"].eq("test")]
        if joint_meta.empty or random_meta.empty:
            continue
        joint_n = len(joint_meta)
        joint_keys = list(zip(joint_meta["cell_context"], joint_meta["de_bin"]))

        per_rows = []
        for subdir in ["ridge_baselines", "mean_baselines"]:
            per_path = seed_dir / subdir / ("ridge_baseline_per_row.csv" if subdir == "ridge_baselines" else "simple_baseline_per_row.csv")
            df = pd.read_csv(per_path)
            per_rows.append(df)
        per = pd.concat(per_rows, ignore_index=True)
        per = per.merge(splits[["sample_id", "cell_context", "n_de_genes", "de_bin"]], on="sample_id", how="left")

        for baseline in ["ridge_cell_drug_fp", "global_train_mean", "drug_mean"]:
            random_per = per.loc[per["split"].eq("split_random") & per["baseline"].eq(baseline)].copy()
            joint_per = per.loc[per["split"].eq("split_cell_scaffold_heldout") & per["baseline"].eq(baseline)].copy()
            if random_per.empty or joint_per.empty:
                continue

            records.append(
                {
                    "seed": seed,
                    "baseline": baseline,
                    "control": "actual_joint",
                    "n_test": len(joint_per),
                    "mean_rowwise_pearson": joint_per["rowwise_pearson"].mean(),
                    "mean_de_rowwise_pearson": joint_per["de_rowwise_pearson"].mean(),
                }
            )

            size_vals = []
            comp_vals = []
            grouped = {}
            by_cell = {}
            for key, group in random_per.groupby(["cell_context", "de_bin"], dropna=False):
                grouped[key] = group
            for key, group in random_per.groupby("cell_context", dropna=False):
                by_cell[key] = group
            for _ in range(n_resamples):
                size_sample = random_per.sample(n=joint_n, replace=joint_n > len(random_per), random_state=int(rng.integers(0, 2**31 - 1)))
                size_vals.append((size_sample["rowwise_pearson"].mean(), size_sample["de_rowwise_pearson"].mean()))

                chosen = []
                for cell, de_bin in joint_keys:
                    pool = grouped.get((cell, de_bin))
                    if pool is None or pool.empty:
                        pool = by_cell.get(cell)
                    if pool.empty:
                        pool = random_per
                    pick = pool.sample(n=1, random_state=int(rng.integers(0, 2**31 - 1)))
                    chosen.append(pick)
                comp_sample = pd.concat(chosen, ignore_index=False)
                comp_vals.append((comp_sample["rowwise_pearson"].mean(), comp_sample["de_rowwise_pearson"].mean()))

            for label, vals in [("size_matched_random", size_vals), ("composition_matched_random", comp_vals)]:
                arr = np.asarray(vals, dtype=float)
                records.append(
                    {
                        "seed": seed,
                        "baseline": baseline,
                        "control": label,
                        "n_test": joint_n,
                        "mean_rowwise_pearson": float(np.nanmean(arr[:, 0])),
                        "mean_de_rowwise_pearson": float(np.nanmean(arr[:, 1])),
                    }
                )

    by_seed = pd.DataFrame(records)
    by_seed.to_csv(OUT / "openproblems_matched_random_controls_by_seed.csv", index=False)
    summary_rows = []
    for (baseline, control), g in by_seed.groupby(["baseline", "control"]):
        lo, hi = ci(g["mean_rowwise_pearson"])
        delo, dehi = ci(g["mean_de_rowwise_pearson"])
        summary_rows.append(
            {
                "baseline": baseline,
                "control": control,
                "seeds": g["seed"].nunique(),
                "mean_rowwise_pearson": g["mean_rowwise_pearson"].mean(),
                "rowwise_pearson_ci_low": lo,
                "rowwise_pearson_ci_high": hi,
                "mean_de_rowwise_pearson": g["mean_de_rowwise_pearson"].mean(),
                "de_rowwise_pearson_ci_low": delo,
                "de_rowwise_pearson_ci_high": dehi,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "openproblems_matched_random_controls_summary.csv", index=False)

    contrasts = []
    pivot = by_seed.pivot_table(
        index=["seed", "baseline"],
        columns="control",
        values=["mean_rowwise_pearson", "mean_de_rowwise_pearson"],
    )
    for baseline in by_seed["baseline"].unique():
        p = pivot.xs(baseline, level="baseline")
        for control in ["size_matched_random", "composition_matched_random"]:
            diff = p[("mean_rowwise_pearson", control)] - p[("mean_rowwise_pearson", "actual_joint")]
            de_diff = p[("mean_de_rowwise_pearson", control)] - p[("mean_de_rowwise_pearson", "actual_joint")]
            lo, hi = ci(diff)
            delo, dehi = ci(de_diff)
            contrasts.append(
                {
                    "baseline": baseline,
                    "contrast": f"{control} - actual_joint",
                    "mean_rowwise_pearson_difference": diff.mean(),
                    "rowwise_difference_ci_low": lo,
                    "rowwise_difference_ci_high": hi,
                    "mean_de_rowwise_pearson_difference": de_diff.mean(),
                    "de_difference_ci_low": delo,
                    "de_difference_ci_high": dehi,
                }
            )
    pd.DataFrame(contrasts).to_csv(OUT / "openproblems_matched_random_controls_contrasts.csv", index=False)


def normalize_name(x: object) -> str:
    return str(x).strip().lower()


def moa_overlap() -> None:
    moa = pd.read_csv(ROOT / "data/raw/openproblems_neurips2023/moa_annotations.csv")
    moa["drug_key"] = moa["sm_name"].map(normalize_name)
    moa = moa.dropna(subset=["moa"]).drop_duplicates("drug_key")
    rows = []
    for seed_dir in sorted((ROOT / "results/openproblems_multiseed_100_de").glob("seed_*")):
        seed = int(seed_dir.name.split("_")[1])
        splits = pd.read_csv(seed_dir / "splits/candidate_splits.csv")
        splits = splits.loc[splits["condition"].eq("treated")].copy()
        splits["drug_key"] = splits["drug_name"].map(normalize_name)
        splits = splits.merge(moa[["drug_key", "moa"]], on="drug_key", how="left")
        for split_col in SPLIT_COLS:
            train = splits.loc[splits[split_col].eq("train")]
            test = splits.loc[splits[split_col].eq("test")]
            train_moas = set(train["moa"].dropna())
            annotated = test["moa"].notna()
            same_moa = test["moa"].isin(train_moas) & annotated
            rows.append(
                {
                    "seed": seed,
                    "split": split_col,
                    "n_test": len(test),
                    "n_test_moa_annotated": int(annotated.sum()),
                    "fraction_test_moa_annotated": float(annotated.mean()) if len(test) else np.nan,
                    "same_moa_in_train_among_all": float(same_moa.mean()) if len(test) else np.nan,
                    "same_moa_in_train_among_annotated": float(same_moa.sum() / annotated.sum()) if annotated.sum() else np.nan,
                }
            )
    by_seed = pd.DataFrame(rows)
    by_seed.to_csv(OUT / "openproblems_moa_overlap_by_seed.csv", index=False)
    summary = (
        by_seed.groupby("split", as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_test_records=("n_test", "mean"),
            mean_annotated_fraction=("fraction_test_moa_annotated", "mean"),
            mean_same_moa_all=("same_moa_in_train_among_all", "mean"),
            sd_same_moa_all=("same_moa_in_train_among_all", "std"),
            mean_same_moa_annotated=("same_moa_in_train_among_annotated", "mean"),
            sd_same_moa_annotated=("same_moa_in_train_among_annotated", "std"),
        )
    )
    summary.to_csv(OUT / "openproblems_moa_overlap_summary.csv", index=False)


def sciplex_dose_leakage() -> None:
    rows = []
    for seed_dir in sorted((ROOT / "results/sciplex3_24h_top2000_multiseed_30").glob("seed_*")):
        seed = int(seed_dir.name.split("_")[1])
        splits = pd.read_csv(seed_dir / "splits/candidate_splits.csv")
        for split_col in SPLIT_COLS:
            train = splits.loc[splits[split_col].eq("train")]
            test = splits.loc[splits[split_col].eq("test")]
            train_drugs = set(train["drug_id"])
            train_drug_doses = set(zip(train["drug_id"], train["dose_value"]))
            train_drug_cell_doses = set(zip(train["drug_id"], train["cell_context"], train["dose_value"]))
            train_doses = set(train["dose_value"])
            same_drug = test["drug_id"].isin(train_drugs)
            same_drug_dose = pd.Series(
                [(d, dose) in train_drug_doses for d, dose in zip(test["drug_id"], test["dose_value"])],
                index=test.index,
            )
            same_drug_cell_dose = pd.Series(
                [
                    (d, c, dose) in train_drug_cell_doses
                    for d, c, dose in zip(test["drug_id"], test["cell_context"], test["dose_value"])
                ],
                index=test.index,
            )
            same_dose = test["dose_value"].isin(train_doses)
            rows.append(
                {
                    "seed": seed,
                    "split": split_col,
                    "n_test": len(test),
                    "same_drug_in_train": float(same_drug.mean()) if len(test) else np.nan,
                    "same_drug_dose_in_train": float(same_drug_dose.mean()) if len(test) else np.nan,
                    "same_drug_cell_dose_in_train": float(same_drug_cell_dose.mean()) if len(test) else np.nan,
                    "same_numeric_dose_in_train": float(same_dose.mean()) if len(test) else np.nan,
                    "unique_test_doses": ",".join(map(str, sorted(test["dose_value"].dropna().unique()))),
                }
            )
    by_seed = pd.DataFrame(rows)
    by_seed.to_csv(OUT / "sciplex3_dose_leakage_by_seed.csv", index=False)
    summary = (
        by_seed.groupby("split", as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_test_records=("n_test", "mean"),
            mean_same_drug=("same_drug_in_train", "mean"),
            mean_same_drug_dose=("same_drug_dose_in_train", "mean"),
            mean_same_drug_cell_dose=("same_drug_cell_dose_in_train", "mean"),
            mean_same_numeric_dose=("same_numeric_dose_in_train", "mean"),
        )
    )
    summary.to_csv(OUT / "sciplex3_dose_leakage_summary.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matched_random_controls()
    moa_overlap()
    sciplex_dose_leakage()
    print(f"Wrote v5 statistics to {OUT}")


if __name__ == "__main__":
    main()
