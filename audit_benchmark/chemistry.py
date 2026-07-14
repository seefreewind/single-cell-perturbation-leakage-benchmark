"""Chemistry feature helpers for benchmark audits."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ChemistryStatus:
    rdkit_available: bool
    scaffold_method: str
    fingerprint_method: str
    radius: int
    n_bits: int
    notes: str = ""


def _load_rdkit() -> tuple[Any, Any, Any, Any]:
    try:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import AllChem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except Exception as exc:  # pragma: no cover - depends on optional runtime
        raise RuntimeError(
            "RDKit is required for scaffold and Tanimoto audits. "
            "Install rdkit-pypi or conda-forge rdkit."
        ) from exc
    RDLogger.DisableLog("rdApp.*")
    return Chem, DataStructs, AllChem, MurckoScaffold


def rdkit_status(radius: int = 2, n_bits: int = 2048) -> ChemistryStatus:
    try:
        _load_rdkit()
    except RuntimeError as exc:
        return ChemistryStatus(
            rdkit_available=False,
            scaffold_method="unavailable",
            fingerprint_method="unavailable",
            radius=radius,
            n_bits=n_bits,
            notes=str(exc),
        )
    return ChemistryStatus(
        rdkit_available=True,
        scaffold_method="Bemis-Murcko scaffold from RDKit MurckoScaffoldSmiles",
        fingerprint_method="Morgan bit vector with Tanimoto similarity",
        radius=radius,
        n_bits=n_bits,
    )


def scaffold_from_smiles(smiles: object) -> str | float:
    if not isinstance(smiles, str) or not smiles.strip():
        return np.nan
    Chem, _, _, MurckoScaffold = _load_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    return scaffold if scaffold else np.nan


def fingerprint_from_smiles(smiles: object, radius: int, n_bits: int):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    Chem, _, AllChem, _ = _load_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def max_tanimoto(test_fp, train_fps: list) -> float:
    if test_fp is None or not train_fps:
        return float("nan")
    _, DataStructs, _, _ = _load_rdkit()
    sims = DataStructs.BulkTanimotoSimilarity(test_fp, train_fps)
    return float(max(sims)) if sims else float("nan")
