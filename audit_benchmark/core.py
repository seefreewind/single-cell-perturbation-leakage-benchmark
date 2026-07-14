"""Core leakage-aware neighborhood audit implementation."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .chemistry import fingerprint_from_smiles, max_tanimoto, rdkit_status, scaffold_from_smiles


@dataclass(frozen=True)
class AuditColumns:
    record_id: str
    drug: str
    smiles: str | None
    cell: str | None
    dose: str | None
    moa: str | None
    scaffold: str | None


def _norm_missing(value: str | None) -> str | None:
    if value is None or value == "" or str(value).lower() in {"none", "nan", "na"}:
        return None
    return value


def _require(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [c for c in columns if c and c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def load_records_and_manifest(
    records_path: Path,
    split_manifest_path: Path,
    columns: AuditColumns,
    split_column: str,
    assignment_column: str,
) -> pd.DataFrame:
    records = pd.read_csv(records_path)
    manifest = pd.read_csv(split_manifest_path)
    _require(records, [columns.record_id, columns.drug], "records")
    _require(manifest, [columns.record_id], "split manifest")

    if assignment_column not in manifest.columns:
        if split_column in manifest.columns:
            assignment_column = split_column
        else:
            raise ValueError(
                f"split manifest must contain assignment column '{assignment_column}' "
                f"or split column '{split_column}'"
            )

    optional_cols = [columns.smiles, columns.cell, columns.dose, columns.moa, columns.scaffold]
    meta_cols = [columns.record_id, columns.drug] + [c for c in optional_cols if c and c in records.columns]
    manifest_cols = [columns.record_id, assignment_column]
    if split_column in manifest.columns and split_column != assignment_column:
        manifest_cols.append(split_column)
    merged = records[meta_cols].merge(manifest[manifest_cols], on=columns.record_id, how="inner")
    if merged.empty:
        raise ValueError("No records matched between records and split manifest.")
    merged = merged.rename(columns={assignment_column: "assignment"})
    merged["assignment"] = merged["assignment"].astype(str).str.lower()
    allowed = {"train", "test", "excluded"}
    bad = sorted(set(merged["assignment"]) - allowed)
    if bad:
        raise ValueError(f"Unsupported assignment values: {bad}; expected train/test/excluded.")
    return merged


def prepare_chemistry(df: pd.DataFrame, columns: AuditColumns, radius: int, n_bits: int) -> tuple[pd.DataFrame, dict]:
    out = df.copy()
    status = rdkit_status(radius=radius, n_bits=n_bits)
    if columns.scaffold and columns.scaffold in out.columns:
        out["_audit_scaffold"] = out[columns.scaffold]
        scaffold_source = "provided"
    elif columns.smiles and columns.smiles in out.columns:
        out["_audit_scaffold"] = out[columns.smiles].map(scaffold_from_smiles)
        scaffold_source = "computed_from_smiles"
    else:
        out["_audit_scaffold"] = np.nan
        scaffold_source = "missing"

    fp_map = {}
    if columns.smiles and columns.smiles in out.columns:
        for drug, smiles in out[[columns.drug, columns.smiles]].drop_duplicates(columns.drug).itertuples(index=False):
            fp_map[drug] = fingerprint_from_smiles(smiles, radius=radius, n_bits=n_bits)

    meta = asdict(status)
    meta.update(
        {
            "scaffold_source": scaffold_source,
            "unique_drugs": int(out[columns.drug].nunique(dropna=True)),
            "fingerprints_available": int(sum(v is not None for v in fp_map.values())),
            "smiles_column_present": bool(columns.smiles and columns.smiles in out.columns),
        }
    )
    return out, {"fingerprints": fp_map, "metadata": meta}


def audit_dataframe(df: pd.DataFrame, columns: AuditColumns, radius: int = 2, n_bits: int = 2048) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    df, chem = prepare_chemistry(df, columns, radius=radius, n_bits=n_bits)
    train = df[df["assignment"].eq("train")].copy()
    test = df[df["assignment"].eq("test")].copy()
    excluded = df[df["assignment"].eq("excluded")].copy()
    if train.empty or test.empty:
        raise ValueError("Both train and test assignments are required for audit.")

    train_drugs = set(train[columns.drug].dropna())
    train_scaffolds = set(train["_audit_scaffold"].dropna())
    train_cells = set(train[columns.cell].dropna()) if columns.cell and columns.cell in train.columns else set()
    train_doses = set(train[columns.dose].dropna()) if columns.dose and columns.dose in train.columns else set()
    train_moa = set(train[columns.moa].dropna()) if columns.moa and columns.moa in train.columns else set()
    train_drug_dose = set()
    if columns.dose and columns.dose in train.columns:
        train_drug_dose = set(zip(train[columns.drug], train[columns.dose]))

    fps = chem["fingerprints"]
    train_fps = [fps[d] for d in train_drugs if fps.get(d) is not None]

    rows = []
    for _, row in test.iterrows():
        drug = row[columns.drug]
        scaffold = row["_audit_scaffold"]
        cell = row[columns.cell] if columns.cell and columns.cell in row.index else np.nan
        dose = row[columns.dose] if columns.dose and columns.dose in row.index else np.nan
        moa = row[columns.moa] if columns.moa and columns.moa in row.index else np.nan
        moa_missing = pd.isna(moa) or moa == ""
        rows.append(
            {
                "record_id": row[columns.record_id],
                "assignment": "test",
                "drug": drug,
                "cell": cell,
                "dose": dose,
                "moa": moa,
                "scaffold": scaffold,
                "same_drug_in_train": bool(drug in train_drugs),
                "same_scaffold_in_train": bool(pd.notna(scaffold) and scaffold in train_scaffolds),
                "same_cell_in_train": bool(pd.notna(cell) and cell in train_cells),
                "same_dose_in_train": bool(pd.notna(dose) and dose in train_doses),
                "same_drug_dose_in_train": bool((drug, dose) in train_drug_dose),
                "moa_missing": bool(moa_missing),
                "same_moa_in_train": bool((not moa_missing) and moa in train_moa),
                "nearest_train_drug_tanimoto": max_tanimoto(fps.get(drug), train_fps),
            }
        )
    record = pd.DataFrame(rows)
    summary = pd.DataFrame(
        [
            {
                "n_records": len(df),
                "n_train": len(train),
                "n_test": len(test),
                "n_excluded": len(excluded),
                "same_drug_overlap": record["same_drug_in_train"].mean(),
                "same_scaffold_overlap": record["same_scaffold_in_train"].mean(),
                "same_cell_overlap": record["same_cell_in_train"].mean(),
                "same_dose_overlap": record["same_dose_in_train"].mean(),
                "same_drug_dose_overlap": record["same_drug_dose_in_train"].mean(),
                "same_moa_overlap_among_all": record["same_moa_in_train"].mean(),
                "moa_missing_rate": record["moa_missing"].mean(),
                "mean_nearest_train_drug_tanimoto": record["nearest_train_drug_tanimoto"].mean(),
                "median_nearest_train_drug_tanimoto": record["nearest_train_drug_tanimoto"].median(),
            }
        ]
    )
    integrity = {
        "passed": True,
        "n_records": int(len(df)),
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "n_excluded": int(len(excluded)),
        "test_records_in_train": int(set(test[columns.record_id]) & set(train[columns.record_id]) != set()),
        "excluded_records_in_train": int(set(excluded[columns.record_id]) & set(train[columns.record_id]) != set()),
        "excluded_records_in_test": int(set(excluded[columns.record_id]) & set(test[columns.record_id]) != set()),
    }
    integrity["passed"] = not any(
        integrity[k] for k in ["test_records_in_train", "excluded_records_in_train", "excluded_records_in_test"]
    )
    return record, summary, chem["metadata"], integrity


def make_size_matched_random_manifest(
    df: pd.DataFrame,
    columns: AuditColumns,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_train = int(df["assignment"].eq("train").sum())
    n_test = int(df["assignment"].eq("test").sum())
    n_excluded = int(df["assignment"].eq("excluded").sum())
    ids = df[columns.record_id].to_numpy()
    if n_train + n_test + n_excluded != len(ids):
        raise ValueError("Assignment counts do not sum to record count.")
    shuffled = ids.copy()
    rng.shuffle(shuffled)
    assignments = {}
    for rid in shuffled[:n_test]:
        assignments[rid] = "test"
    for rid in shuffled[n_test : n_test + n_excluded]:
        assignments[rid] = "excluded"
    for rid in shuffled[n_test + n_excluded :]:
        assignments[rid] = "train"
    return pd.DataFrame({columns.record_id: ids, "assignment": [assignments[rid] for rid in ids]})


def write_outputs(
    outdir: Path,
    record: pd.DataFrame,
    summary: pd.DataFrame,
    chemistry_metadata: dict,
    integrity: dict,
    matched_random: pd.DataFrame | None = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    record.to_csv(outdir / "record_level_audit.csv", index=False)
    summary.to_csv(outdir / "split_level_summary.csv", index=False)
    (outdir / "similarity_matrix_metadata.json").write_text(json.dumps(chemistry_metadata, indent=2))
    (outdir / "integrity_report.json").write_text(json.dumps(integrity, indent=2))
    if matched_random is not None:
        matched_random.to_csv(outdir / "matched_random_manifest.csv", index=False)
    html = [
        "<html><head><title>Benchmark Audit Report</title></head><body>",
        "<h1>Leakage-aware benchmark audit report</h1>",
        "<h2>Split-level summary</h2>",
        summary.to_html(index=False),
        "<h2>Integrity report</h2>",
        f"<pre>{json.dumps(integrity, indent=2)}</pre>",
        "<h2>Chemistry metadata</h2>",
        f"<pre>{json.dumps(chemistry_metadata, indent=2)}</pre>",
        "</body></html>",
    ]
    (outdir / "audit_report.html").write_text("\n".join(html))
