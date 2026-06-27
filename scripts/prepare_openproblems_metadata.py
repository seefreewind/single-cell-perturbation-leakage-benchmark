#!/usr/bin/env python3
"""Prepare project metadata from OpenProblems NeurIPS 2023 DGE h5ad files."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_obs(path: Path, source_split: str) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs.copy().reset_index(names="obs_id")
    obs["source_file"] = path.name
    obs["source_split"] = source_split
    return obs


def prepare(de_train: Path, de_test: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    obs = pd.concat([read_obs(de_train, "de_train"), read_obs(de_test, "de_test")], ignore_index=True)

    required = {"sm_name", "SMILES", "cell_type", "split", "control"}
    missing = required - set(obs.columns)
    if missing:
        raise SystemExit(f"Missing OpenProblems columns: {sorted(missing)}")

    obs["drug_id"] = obs["sm_name"].astype(str)
    obs["drug_name"] = obs["sm_name"].astype(str)
    obs["cell_context"] = obs["cell_type"].astype(str)
    obs["dataset"] = "openproblems_neurips2023"
    obs["cell_line"] = "PBMC"
    obs["tissue"] = "blood"
    obs["batch"] = obs["source_split"].astype(str) + ":" + obs["split"].astype(str)
    obs["dose"] = "not_available_in_de_h5ad"
    obs["time"] = "24h"
    obs["condition"] = obs["control"].map({True: "control", False: "treated"}).fillna("unknown")
    obs["sample_id"] = obs["obs_id"].astype(str)
    obs["control_sample_id"] = ""

    perturbation = obs[
        [
            "sample_id",
            "dataset",
            "drug_id",
            "drug_name",
            "cell_context",
            "cell_type",
            "cell_line",
            "tissue",
            "batch",
            "dose",
            "time",
            "condition",
            "control_sample_id",
            "source_file",
            "source_split",
            "split",
            "control",
        ]
    ].copy()
    perturbation.to_csv(outdir / "openproblems_perturbation_metadata_full.csv", index=False)

    drug = (
        obs[["drug_id", "drug_name", "SMILES"]]
        .rename(columns={"SMILES": "smiles"})
        .drop_duplicates("drug_id")
        .sort_values("drug_id")
    )
    drug["inchikey"] = ""
    drug["pubchem_cid"] = ""
    drug["chembl_id"] = ""
    drug["drugbank_id"] = ""
    drug["target"] = ""
    drug["moa"] = ""
    drug["atc_class"] = ""
    drug.to_csv(outdir / "openproblems_drug_metadata_full.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_rows": len(obs),
                "n_drugs": obs["drug_id"].nunique(),
                "n_cell_contexts": obs["cell_context"].nunique(),
                "n_treated_rows": int((obs["condition"] == "treated").sum()),
                "n_control_rows": int((obs["condition"] == "control").sum()),
                "n_train_rows": int((obs["source_split"] == "de_train").sum()),
                "n_test_rows": int((obs["source_split"] == "de_test").sum()),
            }
        ]
    )
    summary.to_csv(outdir / "openproblems_metadata_summary.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-train", type=Path, required=True)
    parser.add_argument("--de-test", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("metadata/openproblems"))
    args = parser.parse_args()
    prepare(args.de_train, args.de_test, args.outdir)


if __name__ == "__main__":
    main()

