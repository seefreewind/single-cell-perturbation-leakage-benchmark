# Data Source Plan

## Primary benchmark candidates

### 1. Sci-Plex / scPerturb

Recommended role: primary scaffold-strict benchmark candidate.

Rationale:

- Chemical perturbation dataset with multiple cancer cell lines.
- The sci-Plex study profiled A549, K562, and MCF7 cells exposed to 188 compounds or DMSO.
- scPerturb provides harmonized `.h5ad` files through Zenodo.
- `SrivatsanTrapnell2020_sciplex3.h5ad` is the strongest main-data candidate but is large, around 2.35 GB.
- `SrivatsanTrapnell2020_sciplex2.h5ad` and `SrivatsanTrapnell2020_sciplex4.h5ad` are smaller smoke-test candidates.

Immediate action:

1. Download a small Sci-Plex file first when network speed permits.
2. Inspect `.obs` columns to identify drug, dose, cell line, control, and batch fields.
3. Extract harmonized perturbation and drug metadata.
4. Compute scaffolds and test Level 3-5 split feasibility.

Known download issue:

- A Zenodo download attempt on `SrivatsanTrapnell2020_sciplex2.h5ad` was extremely slow, so large downloads should be resumed rather than restarted.

### 2. OpenProblems / NeurIPS 2023 perturbation prediction

Recommended role: secondary benchmark or external validation-style benchmark.

Rationale:

- Predicts how small molecules change gene expression in different cell types.
- Includes PBMCs, multiple cell types, three healthy donors, 144 compounds, SMILES, and LINCS identifiers.
- Provides a ready benchmark structure with `sc_counts.h5ad`, `de_train.h5ad`, `de_test.h5ad`, and `id_map.csv`.

Immediate action:

1. Clone or download the OpenProblems task repository.
2. Inspect resource download scripts and file URLs.
3. Determine whether the resource can be downloaded without Kaggle credentials.
4. If accessible, use it as a clean secondary benchmark with known cell-type and drug metadata.

### 3. McFarlandTsherniak2020 / scPerturb

Recommended role: later cancer-cell-line robustness analysis.

Rationale:

- Contains many cancer cell lines and a smaller number of drugs.
- Useful for cell-context generalization and drug-response transfer, but less ideal for scaffold-strict drug generalization because the number of drugs is limited.

Immediate action:

1. Consider after Sci-Plex metadata is inspected.
2. Use mainly for cell-held-out and drug-sensitivity-linked analyses.

## Data-source decision rules

- Use Sci-Plex as primary if it provides enough scaffold diversity and at least two cell contexts per scaffold.
- Use OpenProblems/NeurIPS as a secondary benchmark if it can be downloaded without manual credentials.
- Use McFarland as a cell-context generalization analysis if drug diversity is too limited for Level 4-5 scaffold testing.

## Minimum metadata required

For every included dataset, construct:

- `metadata/perturbation_metadata.csv`
- `metadata/drug_metadata.csv`
- `results/feasibility/feasibility_summary.csv`
- `results/feasibility/scaffold_cell_context_coverage.csv`

The benchmark should not proceed to model training until these files exist for the selected primary dataset.

