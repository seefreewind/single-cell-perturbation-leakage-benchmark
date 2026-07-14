"""Command-line entry point for leakage-aware benchmark auditing."""

from __future__ import annotations

import argparse
from pathlib import Path

from .core import AuditColumns, audit_dataframe, load_records_and_manifest, make_size_matched_random_manifest, write_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m audit_benchmark",
        description="Audit train-test biological and chemical neighborhoods for a benchmark split.",
    )
    parser.add_argument("--records", type=Path, required=True, help="CSV containing benchmark records and metadata.")
    parser.add_argument("--split-manifest", type=Path, required=True, help="CSV containing record assignments.")
    parser.add_argument("--record-id-column", default="sample_id")
    parser.add_argument("--drug-column", default="drug_name")
    parser.add_argument("--smiles-column", default="smiles")
    parser.add_argument("--cell-column", default="cell_context")
    parser.add_argument("--dose-column", default="dose")
    parser.add_argument("--moa-column", default="moa")
    parser.add_argument("--scaffold-column", default="scaffold")
    parser.add_argument("--split-column", default="assignment")
    parser.add_argument("--assignment-column", default="assignment")
    parser.add_argument("--fingerprint-radius", type=int, default=2)
    parser.add_argument("--fingerprint-bits", type=int, default=2048)
    parser.add_argument("--matched-random-seed", type=int, default=1)
    parser.add_argument("--no-matched-random", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    columns = AuditColumns(
        record_id=args.record_id_column,
        drug=args.drug_column,
        smiles=args.smiles_column,
        cell=args.cell_column,
        dose=args.dose_column,
        moa=args.moa_column,
        scaffold=args.scaffold_column,
    )
    df = load_records_and_manifest(
        args.records,
        args.split_manifest,
        columns=columns,
        split_column=args.split_column,
        assignment_column=args.assignment_column,
    )
    record, summary, chemistry_metadata, integrity = audit_dataframe(
        df,
        columns=columns,
        radius=args.fingerprint_radius,
        n_bits=args.fingerprint_bits,
    )
    matched = None if args.no_matched_random else make_size_matched_random_manifest(df, columns, seed=args.matched_random_seed)
    write_outputs(args.output_dir, record, summary, chemistry_metadata, integrity, matched_random=matched)
    print(f"Wrote audit outputs to {args.output_dir}")
