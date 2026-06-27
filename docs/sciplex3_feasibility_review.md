# Sci-Plex 3 Feasibility Review

## Bottom Line

Sci-Plex 3 from scPerturb is a strong second-dataset candidate for the high-score revision. It is independent of the OpenProblems/OP3 PBMC benchmark, has multiple cancer cell lines, many small molecules, multiple doses, and enough valid scaffolds for scaffold-strict and joint cell-plus-scaffold-held-out validation.

## Local Data

Downloaded file:

- `data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad`
- Size: 2,526,631,614 bytes
- Shape: 799,317 cells x 110,983 genes/features

SMILES mapping:

- `metadata/sciplex3_pubchem_smiles.csv`
- `metadata/sciplex3_pubchem_smiles_with_scaffolds.csv`
- PubChem mapping success: 188/188 drug perturbations
- RDKit-valid molecules: 188/188
- Unique Bemis-Murcko scaffolds: 165

Candidate record table:

- `metadata/sciplex3_candidate_records.csv`
- Unit: perturbation x cell_line x dose_value x time

## Metadata Summary

Cell-line composition:

| Cell line | Cells |
|---|---:|
| MCF7 | 344,862 |
| A549 | 244,281 |
| K562 | 173,652 |
| missing cell_line | 36,522 |

Perturbation structure:

| Field | Count |
|---|---:|
| Drug perturbations, excluding control | 188 |
| Control cells | 17,578 |
| Dose levels | 10, 100, 1000, 10000 nM |
| Time points | 24 h and 72 h |
| Replicates | rep1 and rep2 |

Controls by cell line and time:

| Cell line | Time | Control cells |
|---|---:|---:|
| A549 | 24 h | 3,773 |
| A549 | 72 h | 2,084 |
| K562 | 24 h | 3,935 |
| MCF7 | 24 h | 7,786 |

## Candidate Records

Candidate records are defined as drug x cell_line x dose x time groups. The table below shows how many records remain after filtering by minimum cells per group.

| Minimum cells per group | Records | Drugs | Cell lines | Scaffolds |
|---:|---:|---:|---:|---:|
| 10 | 2,444 | 188 | 3 | 165 |
| 20 | 2,440 | 188 | 3 | 165 |
| 50 | 2,423 | 188 | 3 | 165 |
| 100 | 2,377 | 188 | 3 | 165 |

With a 100-cell threshold:

| Cell line | Records |
|---|---:|
| A549 | 914 |
| K562 | 719 |
| MCF7 | 744 |

Time-point structure:

| Time | Records at min 100 cells |
|---:|---:|
| 24 h | 2,200 |
| 72 h | 177 |

## Recommended Benchmark Scope

Primary Sci-Plex 3 benchmark:

- Use 24 h records only.
- Include MCF7, A549, and K562.
- Use drug x cell_line x dose as the prediction record.
- Require at least 100 cells per treated group.
- Compare treated pseudobulk expression to matched control pseudobulk within cell line and time.
- Treat dose as either:
  - a separate condition feature, or
  - restrict to one standard dose, likely 1000 nM, for the cleanest first replication.

Recommended first pass:

> 24 h, all three cell lines, all four doses, minimum 100 cells per treated group.

This keeps 2,200 candidate records and preserves enough repeated drug-cell-dose structure to test random, cell-line-held-out, scaffold-held-out, and joint cell-line-plus-scaffold-held-out validation.

## Why This Solves The Reviewer Risk

Sci-Plex 3 addresses the largest current critique: the OpenProblems/OP3 analysis is one dataset. Compared with the current PBMC benchmark, Sci-Plex 3 adds:

- a different biological system: cancer cell lines rather than PBMC cell types;
- a larger perturbation grid: 188 drugs and 2,377 high-cell-count records;
- a large scaffold set: 165 unique Bemis-Murcko scaffolds;
- dose structure, enabling either dose-aware modeling or dose-restricted robustness checks.

## Remaining Work

## Completed Smoke-test Benchmark

Processed response matrix:

- `data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad`
- Shape: 2,200 records x 2,000 genes
- Response definition: treated mean log1p(CP10K) minus matched cell-line control mean log1p(CP10K)
- Scope: 24 h, A549/K562/MCF7, minimum 100 cells per treated group, top 2,000 genes by total count

Split output:

- `results/sciplex3_24h_top2000_seed_001/splits/candidate_splits.csv`

Seed 1 split sizes:

| Split | Train records | Test records | Excluded records |
|---|---:|---:|---:|
| Random | 1,767 | 433 | 0 |
| Cell-line held-out | 1,463 | 737 | 0 |
| Scaffold held-out | 1,766 | 434 | 0 |
| Cell-line + scaffold held-out | 1,173 | 144 | 883 |

Smoke-test baseline result for `ridge_cell_dose_drug_fp`:

| Split | Mean row-wise Pearson |
|---|---:|
| Random | 0.278 |
| Cell-line held-out | 0.149 |
| Scaffold held-out | 0.167 |
| Cell-line + scaffold held-out | 0.061 |

Chemical-neighbor audit for the same baseline:

| Split | Mean nearest training-drug Tanimoto |
|---|---:|
| Random | 1.000 |
| Cell-line held-out | 1.000 |
| Scaffold held-out | 0.283 |
| Cell-line + scaffold held-out | 0.284 |

Interpretation:

The first Sci-Plex 3 smoke test reproduces the OpenProblems/OP3 pattern in an independent dataset. Random validation keeps identical drugs in the training set and gives the highest drug-fingerprint ridge performance. Scaffold-strict and joint validation remove same-scaffold chemical neighbors and substantially reduce performance.

## Remaining Work

1. Replace the smoke-test ridge solver with a numerically more stable implementation.
2. Run 30-100 repeated Sci-Plex 3 split seeds.
3. Add top-k or pathway-focused response metrics if needed.
4. Decide whether the main Sci-Plex analysis should use all doses or a standard-dose subset.
5. Add a two-dataset comparison figure/table to the manuscript.
