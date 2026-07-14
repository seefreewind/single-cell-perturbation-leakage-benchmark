# Leakage-aware benchmark audit CLI

This project now includes a minimal reusable command-line audit tool:

```bash
python -m audit_benchmark \
  --records metadata/examples/splits/candidate_splits.csv \
  --split-manifest metadata/examples/splits/candidate_splits.csv \
  --split-column split_cell_scaffold_heldout \
  --assignment-column split_cell_scaffold_heldout \
  --output-dir results/example_audit_cli
```

The command writes:

- `record_level_audit.csv`
- `split_level_summary.csv`
- `similarity_matrix_metadata.json`
- `integrity_report.json`
- `matched_random_manifest.csv`
- `audit_report.html`

## Input requirements

The records table must include a record ID column, drug identifier/name, and enough metadata to define the requested audit dimensions. The default columns follow the project metadata examples:

| Argument | Default |
|---|---|
| `--record-id-column` | `sample_id` |
| `--drug-column` | `drug_name` |
| `--smiles-column` | `smiles` |
| `--cell-column` | `cell_context` |
| `--dose-column` | `dose` |
| `--moa-column` | `moa` |
| `--scaffold-column` | `scaffold` |

The split manifest must contain the record ID column and one assignment column with values `train`, `test`, or `excluded`.

## Chemistry implementation

The CLI uses RDKit to compute Bemis-Murcko scaffolds when a scaffold column is not supplied, Morgan fingerprints, and nearest-training-drug Tanimoto similarity. The default fingerprint parameters are radius 2 and 2048 bits.

## Current scope

This is a minimal reusable audit CLI for tabular benchmark records and fixed split manifests. It does not train models, run third-party adapters, or replace the project-specific scripts used to reproduce all manuscript figures.
