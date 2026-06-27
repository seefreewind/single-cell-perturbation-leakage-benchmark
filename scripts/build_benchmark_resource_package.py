#!/usr/bin/env python3
"""Build reusable benchmark-resource artifacts for the manuscript package."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESOURCE = ROOT / "benchmark_resource"
RESULTS = ROOT / "results/benchmark_resource"

SPLIT_COLUMNS = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]


def seed_from_dir(path: Path) -> int:
    return int(path.name.split("_")[1])


def dataset_label(path: Path) -> str:
    if "sciplex3" in str(path):
        return "sciplex3_24h_top2000"
    return "openproblems_neurips2023"


def build_split_manifest(seed_dirs: list[Path], out_path: Path) -> pd.DataFrame:
    rows = []
    keep = [
        "sample_id",
        "drug_id",
        "drug_name",
        "cell_context",
        "condition",
        "smiles",
        "scaffold",
        "dose",
        "dose_value",
        "time",
    ]
    for seed_dir in seed_dirs:
        seed = seed_from_dir(seed_dir)
        dataset = dataset_label(seed_dir)
        splits = pd.read_csv(seed_dir / "splits/candidate_splits.csv")
        present_keep = [c for c in keep if c in splits.columns]
        for split_col in SPLIT_COLUMNS:
            base = splits[present_keep].copy()
            base["dataset"] = dataset
            base["seed"] = seed
            base["split_name"] = split_col.replace("split_", "")
            base["assignment"] = splits[split_col].values
            rows.append(base)
    manifest = pd.concat(rows, ignore_index=True)
    ordered = ["dataset", "seed", "split_name", "sample_id", "assignment"]
    rest = [c for c in manifest.columns if c not in ordered]
    manifest = manifest[ordered + rest]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out_path, index=False)
    return manifest


def summarize_manifest(manifest: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    summary = (
        manifest.groupby(["dataset", "seed", "split_name", "assignment"], as_index=False)
        .size()
        .rename(columns={"size": "records"})
    )
    summary.to_csv(out_path, index=False)
    return summary


def combine_leakage_tables(out_path: Path) -> pd.DataFrame:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    parts = []
    open_overlap = pd.read_csv(ROOT / "results/leakage_overlap_100/all_leakage_overlap_by_seed.csv")
    open_overlap["dataset"] = "openproblems_neurips2023"
    parts.append(open_overlap)

    sciplex_overlap = pd.read_csv(ROOT / "results/sciplex3_24h_top2000_multiseed_30/leakage_overlap_by_seed.csv")
    sciplex_overlap["dataset"] = "sciplex3_24h_top2000"
    parts.append(sciplex_overlap)

    overlap = pd.concat(parts, ignore_index=True, sort=False)
    moa = pd.read_csv(ROOT / "results/manuscript_v5_stats/openproblems_moa_overlap_by_seed.csv")
    moa["dataset"] = "openproblems_neurips2023"
    dose = pd.read_csv(ROOT / "results/manuscript_v5_stats/sciplex3_dose_leakage_by_seed.csv")
    dose["dataset"] = "sciplex3_24h_top2000"

    leakage = overlap.merge(
        moa[["dataset", "seed", "split", "same_moa_in_train_among_all"]],
        on=["dataset", "seed", "split"],
        how="left",
    ).merge(
        dose[
            [
                "dataset",
                "seed",
                "split",
                "same_drug_dose_in_train",
                "same_drug_cell_dose_in_train",
                "same_numeric_dose_in_train",
            ]
        ],
        on=["dataset", "seed", "split"],
        how="left",
    )
    leakage.to_csv(out_path, index=False)
    return leakage


def moa_heldout_feasibility(outdir: Path) -> None:
    moa = pd.read_csv(ROOT / "data/raw/openproblems_neurips2023/moa_annotations.csv")
    moa["drug_key"] = moa["sm_name"].str.strip().str.lower()
    moa = moa.dropna(subset=["moa"]).drop_duplicates("drug_key")

    rows = []
    rng = np.random.default_rng(20260623)
    for seed_dir in sorted((ROOT / "results/openproblems_multiseed_100_de").glob("seed_*")):
        seed = seed_from_dir(seed_dir)
        splits = pd.read_csv(seed_dir / "splits/candidate_splits.csv")
        splits = splits.loc[splits["condition"].eq("treated")].copy()
        splits["drug_key"] = splits["drug_name"].str.strip().str.lower()
        splits = splits.merge(moa[["drug_key", "moa"]], on="drug_key", how="left")
        moas = np.array(sorted(splits["moa"].dropna().unique()))
        if len(moas) == 0:
            continue
        n_holdout = max(1, int(round(0.2 * len(moas))))
        heldout = set(rng.choice(moas, size=n_holdout, replace=False))
        test = splits.loc[splits["moa"].isin(heldout)]
        train = splits.loc[~splits["moa"].isin(heldout)]
        train_drugs = set(train["drug_name"])
        train_scaffolds = set(train["scaffold"])
        train_moas = set(train["moa"].dropna())
        rows.append(
            {
                "seed": seed,
                "n_unique_moa": len(moas),
                "n_heldout_moa": len(heldout),
                "train_records": len(train),
                "test_records": len(test),
                "test_drugs": test["drug_name"].nunique(),
                "test_scaffolds": test["scaffold"].nunique(),
                "same_drug_overlap": float(test["drug_name"].isin(train_drugs).mean()) if len(test) else np.nan,
                "same_scaffold_overlap": float(test["scaffold"].isin(train_scaffolds).mean()) if len(test) else np.nan,
                "same_moa_overlap": float(test["moa"].isin(train_moas).mean()) if len(test) else np.nan,
            }
        )
    by_seed = pd.DataFrame(rows)
    by_seed.to_csv(outdir / "openproblems_moa_heldout_feasibility_by_seed.csv", index=False)
    summary = by_seed.agg(
        {
            "n_unique_moa": "mean",
            "n_heldout_moa": "mean",
            "train_records": ["mean", "std"],
            "test_records": ["mean", "std"],
            "test_drugs": ["mean", "std"],
            "test_scaffolds": ["mean", "std"],
            "same_drug_overlap": "mean",
            "same_scaffold_overlap": "mean",
            "same_moa_overlap": "mean",
        }
    )
    summary.to_csv(outdir / "openproblems_moa_heldout_feasibility_summary_matrix.csv")
    flat = {
        "seeds": by_seed["seed"].nunique(),
        "mean_unique_moa": by_seed["n_unique_moa"].mean(),
        "mean_heldout_moa": by_seed["n_heldout_moa"].mean(),
        "mean_train_records": by_seed["train_records"].mean(),
        "sd_train_records": by_seed["train_records"].std(),
        "mean_test_records": by_seed["test_records"].mean(),
        "sd_test_records": by_seed["test_records"].std(),
        "mean_test_drugs": by_seed["test_drugs"].mean(),
        "mean_test_scaffolds": by_seed["test_scaffolds"].mean(),
        "mean_same_drug_overlap": by_seed["same_drug_overlap"].mean(),
        "mean_same_scaffold_overlap": by_seed["same_scaffold_overlap"].mean(),
        "mean_same_moa_overlap": by_seed["same_moa_overlap"].mean(),
    }
    pd.DataFrame([flat]).to_csv(outdir / "openproblems_moa_heldout_feasibility_summary.csv", index=False)


def write_schema() -> None:
    schema = {
        "schema_name": "single_cell_drug_perturbation_leakage_audit",
        "version": "0.1.0",
        "required_fields": {
            "dataset": "Dataset identifier.",
            "seed": "Integer split seed.",
            "split_name": "random, cell_heldout, scaffold_heldout, cell_scaffold_heldout, or another declared split.",
            "n_train": "Number of training records.",
            "n_test": "Number of testing records.",
            "same_drug_overlap": "Fraction of test records with same drug in train.",
            "same_scaffold_overlap": "Fraction of test records with same scaffold in train.",
            "same_cell_context_overlap": "Fraction of test records with same cell context in train.",
            "nearest_train_drug_tanimoto": "Mean nearest training-drug Tanimoto similarity.",
            "model_name": "Model or baseline name.",
            "mean_rowwise_pearson": "Mean row-wise Pearson correlation.",
            "mean_de_rowwise_pearson": "Mean DE-gene row-wise Pearson correlation where available.",
        },
        "recommended_fields": {
            "same_moa_overlap": "Fraction of test records with same annotated MoA in train.",
            "same_drug_dose_overlap": "Fraction of test drug-dose combinations observed in train.",
            "random_to_strict_delta": "Paired random-to-strict performance difference.",
            "rank_correlation_with_random": "Model-rank Spearman correlation with random split.",
            "external_replication_dataset": "Independent dataset used for replication.",
        },
    }
    (RESOURCE / "schemas").mkdir(parents=True, exist_ok=True)
    with (RESOURCE / "schemas/benchmark_report_schema.json").open("w", encoding="utf-8") as fh:
        json.dump(schema, fh, indent=2, ensure_ascii=False)

    yaml_lines = [
        "schema_name: single_cell_drug_perturbation_leakage_audit",
        "version: 0.1.0",
        "required_fields:",
    ]
    for key, value in schema["required_fields"].items():
        yaml_lines.append(f"  {key}: {value}")
    yaml_lines.append("recommended_fields:")
    for key, value in schema["recommended_fields"].items():
        yaml_lines.append(f"  {key}: {value}")
    (RESOURCE / "schemas/benchmark_report_schema.yaml").write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")


def write_templates() -> None:
    (RESOURCE / "templates").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        columns=[
            "dataset",
            "seed",
            "split_name",
            "n_train",
            "n_test",
            "same_drug_overlap",
            "same_scaffold_overlap",
            "same_moa_overlap",
            "same_cell_context_overlap",
            "same_drug_dose_overlap",
            "nearest_train_drug_tanimoto",
            "model_name",
            "mean_rowwise_pearson",
            "mean_de_rowwise_pearson",
            "rank_correlation_with_random",
        ]
    ).to_csv(RESOURCE / "templates/leakage_audit_template.csv", index=False)


def write_readme() -> None:
    text = """# Scaffold-strict single-cell perturbation benchmark resource

This local resource package accompanies the manuscript draft and provides reusable split manifests, leakage-audit tables, and a report schema for single-cell drug perturbation prediction benchmarks.

## Contents

- `splits/split_manifest_openproblems.csv`: 100-seed long-format train/test/excluded assignments for OpenProblems.
- `splits/split_manifest_sciplex3.csv`: 30-seed long-format assignments for the Sci-Plex 3 24 h pseudobulk benchmark.
- `audit/leakage_audit_table.csv`: same-drug, same-scaffold, same-cell, same-MoA, and dose-neighbor audit fields where available.
- `schemas/benchmark_report_schema.{json,yaml}`: minimum reporting schema for future benchmark reports.
- `templates/leakage_audit_template.csv`: empty template for applying the audit to new datasets.
- `results/openproblems_moa_heldout_feasibility_summary.csv`: feasibility summary for a secondary MoA-held-out split.

## Reproducibility entry points

Run `python scripts/build_benchmark_resource_package.py` from the project root to regenerate this package.

For formal submission, upload this folder together with the analysis scripts and environment files to a public GitHub repository and archive the release on Zenodo. The manuscript should then replace the local-package note with the final repository URL and DOI.
"""
    (RESOURCE / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    RESOURCE.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    open_dirs = sorted((ROOT / "results/openproblems_multiseed_100_de").glob("seed_*"))
    sciplex_dirs = sorted((ROOT / "results/sciplex3_24h_top2000_multiseed_30").glob("seed_*"))
    build_split_manifest(open_dirs, RESOURCE / "splits/split_manifest_openproblems.csv")
    build_split_manifest(sciplex_dirs, RESOURCE / "splits/split_manifest_sciplex3.csv")
    summarize_manifest(pd.read_csv(RESOURCE / "splits/split_manifest_openproblems.csv"), RESOURCE / "splits/split_manifest_openproblems_summary.csv")
    summarize_manifest(pd.read_csv(RESOURCE / "splits/split_manifest_sciplex3.csv"), RESOURCE / "splits/split_manifest_sciplex3_summary.csv")
    combine_leakage_tables(RESOURCE / "audit/leakage_audit_table.csv")
    moa_heldout_feasibility(RESULTS)
    write_schema()
    write_templates()
    write_readme()
    print(f"Wrote benchmark resource package to {RESOURCE}")
    print(f"Wrote resource analysis results to {RESULTS}")


if __name__ == "__main__":
    main()
