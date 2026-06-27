#!/usr/bin/env python3
"""Feasibility audit for leakage-resistant drug perturbation benchmarks.

Expected inputs:
  --perturbations metadata/perturbation_metadata.csv
  --drugs metadata/drug_metadata.csv

Required perturbation columns:
  drug_id, cell_context

Required drug columns:
  drug_id, smiles

Optional columns are retained when present.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def compute_scaffold(smiles: str) -> str | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:
        raise SystemExit(
            "RDKit is required for scaffold calculation. Install dependencies from env/requirements.txt."
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaffold, canonical=True)


def audit(perturbations_path: Path, drugs_path: Path, outdir: Path) -> None:
    perturbations = pd.read_csv(perturbations_path)
    drugs = pd.read_csv(drugs_path)

    required_perturbation_cols = {"drug_id", "cell_context"}
    required_drug_cols = {"drug_id", "smiles"}
    missing_perturbation = required_perturbation_cols - set(perturbations.columns)
    missing_drug = required_drug_cols - set(drugs.columns)
    if missing_perturbation:
        raise SystemExit(f"Missing perturbation columns: {sorted(missing_perturbation)}")
    if missing_drug:
        raise SystemExit(f"Missing drug columns: {sorted(missing_drug)}")

    outdir.mkdir(parents=True, exist_ok=True)

    drugs = drugs.copy()
    drugs["scaffold"] = drugs["smiles"].map(compute_scaffold)
    drugs.to_csv(outdir / "drug_metadata_with_scaffolds.csv", index=False)

    merged = perturbations.merge(
        drugs[["drug_id", "smiles", "scaffold"]].drop_duplicates("drug_id"),
        on="drug_id",
        how="left",
    )
    merged.to_csv(outdir / "perturbation_metadata_with_scaffolds.csv", index=False)

    scaffold_context = (
        merged.dropna(subset=["scaffold"])
        .groupby("scaffold")
        .agg(
            n_drugs=("drug_id", "nunique"),
            n_cell_contexts=("cell_context", "nunique"),
            n_records=("drug_id", "size"),
        )
        .reset_index()
        .sort_values(["n_cell_contexts", "n_drugs", "n_records"], ascending=False)
    )
    scaffold_context.to_csv(outdir / "scaffold_cell_context_coverage.csv", index=False)

    perturbation_drugs = set(perturbations["drug_id"].dropna())
    overlapping_drugs = drugs[drugs["drug_id"].isin(perturbation_drugs)]

    summary = {
        "n_records": len(perturbations),
        "n_drugs_in_perturbation_metadata": perturbations["drug_id"].nunique(),
        "n_drugs_with_smiles": overlapping_drugs.loc[overlapping_drugs["smiles"].notna(), "drug_id"].nunique(),
        "n_valid_scaffolds": merged["scaffold"].nunique(dropna=True),
        "n_cell_contexts": perturbations["cell_context"].nunique(),
        "n_drug_cell_pairs": merged[["drug_id", "cell_context"]].drop_duplicates().shape[0],
        "n_scaffolds_ge2_cell_contexts": int((scaffold_context["n_cell_contexts"] >= 2).sum()),
        "n_scaffolds_ge3_cell_contexts": int((scaffold_context["n_cell_contexts"] >= 3).sum()),
        "n_scaffolds_ge5_cell_contexts": int((scaffold_context["n_cell_contexts"] >= 5).sum()),
    }

    if "batch" in perturbations.columns:
        summary["n_batches"] = perturbations["batch"].nunique()
    if "dataset" in perturbations.columns:
        summary["n_datasets"] = perturbations["dataset"].nunique()

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(outdir / "feasibility_summary.csv", index=False)

    print(summary_df.to_string(index=False))
    ge2 = summary["n_scaffolds_ge2_cell_contexts"]
    if ge2 >= 50:
        print("Decision: proceed with full Level 4-5 benchmark.")
    elif ge2 >= 30:
        print("Decision: proceed as a pilot benchmark or methods-focused short article.")
    else:
        print("Decision: emphasize leakage-risk auditing rather than strict unseen-drug generalization.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--perturbations", type=Path, required=True)
    parser.add_argument("--drugs", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("results/feasibility"))
    args = parser.parse_args()
    audit(args.perturbations, args.drugs, args.outdir)


if __name__ == "__main__":
    main()
