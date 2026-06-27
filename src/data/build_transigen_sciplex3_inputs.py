"""Build Sci-Plex 3 pseudobulk inputs for the TranSiGen adaptation."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/processed/transigen/sciplex3"
H5AD = ROOT / "data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad"
SPLITS = ROOT / "benchmark_resource/split_manifests/split_manifest_sciplex3.csv"


def morgan_fp(smiles: str, n_bits: int = 1024) -> np.ndarray:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return np.zeros(n_bits, dtype=np.float32)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H5AD)
    if "control_mean_log1p_cp10k" not in adata.layers:
        raise ValueError("Missing control_mean_log1p_cp10k layer.")
    if "treated_mean_log1p_cp10k" not in adata.layers:
        raise ValueError("Missing treated_mean_log1p_cp10k layer.")

    record_ids = adata.obs_names.astype(str).to_numpy()
    basal = np.asarray(adata.layers["control_mean_log1p_cp10k"], dtype=np.float32)
    treated = np.asarray(adata.layers["treated_mean_log1p_cp10k"], dtype=np.float32)
    response = np.asarray(adata.X, dtype=np.float32)
    if response.shape != basal.shape:
        response = treated - basal
    genes = adata.var_names.astype(str).to_list()
    np.savez_compressed(
        OUT / "matrices.npz",
        basal=basal,
        treated=treated,
        response=response,
        record_id=record_ids,
        gene=np.asarray(genes),
    )
    (OUT / "gene_names.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

    meta = adata.obs.reset_index(names="record_id").rename(
        columns={
            "perturbation": "drug_name",
            "canonical_smiles": "smiles",
            "cell_line": "cell_line",
        }
    )
    meta["sample_id"] = meta["record_id"]
    meta["dose_value"] = meta["dose_value"].astype(float)
    meta["time"] = meta["time"].astype(float)
    meta["cell_line_24h_control_profile"] = meta["cell_line"].astype(str) + "_24h_control_profile"
    meta.to_csv(OUT / "record_metadata.csv", index=False)

    compounds = meta[["drug_name", "smiles", "scaffold", "pubchem_cid"]].drop_duplicates("drug_name").copy()
    fps = np.vstack([morgan_fp(s) for s in compounds["smiles"]]).astype(np.float32)
    compounds["fingerprint_row"] = np.arange(len(compounds))
    compounds.to_csv(OUT / "compound_metadata.csv", index=False)
    np.savez_compressed(OUT / "compound_fingerprints.npz", drug_name=compounds["drug_name"].to_numpy(), morgan1024=fps)

    splits = pd.read_csv(SPLITS)
    keep_cols = [
        "dataset",
        "seed",
        "split_name",
        "sample_id",
        "assignment",
        "drug_name",
        "cell_context",
        "smiles",
        "scaffold",
        "dose_value",
        "time",
    ]
    split_long = splits[keep_cols].copy()
    split_long = split_long.rename(columns={"sample_id": "record_id", "cell_context": "cell_line"})
    missing = set(split_long["record_id"]) - set(record_ids)
    if missing:
        raise ValueError(f"Split manifest contains {len(missing)} record_ids absent from response matrix.")
    split_long.to_csv(OUT / "split_assignments_long.csv", index=False)

    readme = f"""# TranSiGen Sci-Plex 3 Input Package

Source response matrix: `{H5AD}`

Records: {response.shape[0]}
Genes: {response.shape[1]}

Matrices in `matrices.npz`:
- `basal`: matched cell-line 24 h control mean log1p(CP10K), shape {basal.shape}
- `treated`: treated mean log1p(CP10K), shape {treated.shape}
- `response`: treated minus matched control response target, shape {response.shape}

Compound fingerprints:
- `compound_fingerprints.npz`, Morgan radius 2, 1024 bits.

Split assignments:
- `split_assignments_long.csv` is copied from the existing Sci-Plex 3 split manifest.
- Excluded records must not enter training or testing.
- Test records must not be used for normalization, early stopping, or hyperparameter selection.

Adaptation note:
This package supports `transigen_adapted_sciplex3`, a pseudobulk adaptation that uses basal/control expression, compound fingerprints, dose, and cell-line features to predict treated-control response vectors.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
