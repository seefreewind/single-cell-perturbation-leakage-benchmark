#!/usr/bin/env python3
"""Map scPerturb perturbation names to PubChem CIDs and SMILES."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd


def candidate_names(name: str) -> list[str]:
    name = str(name).strip()
    candidates = [name]
    no_paren = re.sub(r"\s*\([^)]*\)", "", name).strip()
    if no_paren and no_paren not in candidates:
        candidates.append(no_paren)
    for inside in re.findall(r"\(([^)]*)\)", name):
        cleaned_inside = inside.strip()
        if cleaned_inside and cleaned_inside not in candidates:
            candidates.append(cleaned_inside)
    before_comma = name.split(",")[0].strip()
    if before_comma and before_comma not in candidates:
        candidates.append(before_comma)
    replacements = [
        (" 2HCl", ""),
        (" HCl", ""),
        (" hcl", ""),
        (" hydrochloride", ""),
        (" phosphate", ""),
        (" mesylate", ""),
        (" diphosphate", ""),
        (" sodium", ""),
        (" Sodium", ""),
        (" trihydrate", ""),
    ]
    for old, new in replacements:
        cleaned = name.replace(old, new).strip()
        if cleaned and cleaned not in candidates:
            candidates.append(cleaned)
    return candidates


def fetch_json_with_curl(url: str, timeout: int) -> dict | None:
    try:
        proc = subprocess.run(
            ["curl", "-L", "--silent", "--show-error", "--max-time", str(timeout), url],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError:
        return None
    if proc.returncode != 0 or not proc.stdout.strip():
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None


def pubchem_lookup(name: str, timeout: int = 30) -> dict[str, str | int | None]:
    props = "CanonicalSMILES,IsomericSMILES,IUPACName,MolecularFormula"
    for candidate in candidate_names(name):
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{quote(candidate)}/property/{props}/JSON"
        )
        data = fetch_json_with_curl(url, timeout)
        if data is None:
            time.sleep(0.15)
            continue
        items = data.get("PropertyTable", {}).get("Properties", [])
        if not items:
            time.sleep(0.15)
            continue
        item = items[0]
        canonical = item.get("CanonicalSMILES") or item.get("SMILES") or item.get("ConnectivitySMILES")
        isomeric = item.get("IsomericSMILES") or item.get("SMILES")
        return {
            "query_name": name,
            "matched_name": candidate,
            "pubchem_cid": item.get("CID"),
            "canonical_smiles": canonical,
            "isomeric_smiles": isomeric,
            "iupac_name": item.get("IUPACName"),
            "molecular_formula": item.get("MolecularFormula"),
            "status": "matched",
        }
    return {
        "query_name": name,
        "matched_name": None,
        "pubchem_cid": None,
        "canonical_smiles": None,
        "isomeric_smiles": None,
        "iupac_name": None,
        "molecular_formula": None,
        "status": "unmatched",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=Path("data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad"))
    parser.add_argument("--out", type=Path, default=Path("metadata/sciplex3_pubchem_smiles.csv"))
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs[["perturbation"]].copy()
    adata.file.close()

    names = sorted(
        n
        for n in obs["perturbation"].dropna().astype(str).unique()
        if n.lower() != "control"
    )
    existing = pd.DataFrame()
    done = set()
    if args.out.exists() and not args.force:
        existing = pd.read_csv(args.out)
        done = set(existing.loc[existing["status"].eq("matched"), "query_name"].astype(str))

    rows = [] if existing.empty else existing.to_dict("records")
    for i, name in enumerate(names, start=1):
        if name in done:
            continue
        result = pubchem_lookup(name)
        rows.append(result)
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"[{i}/{len(names)}] {name}: {result['status']}")
        time.sleep(args.sleep)

    result_df = pd.DataFrame(rows)
    result_df.to_csv(args.out, index=False)
    print(result_df["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
