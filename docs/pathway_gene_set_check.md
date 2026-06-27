# Pathway Gene Set Check

Local GMT file: `/Users/zy/Documents/New project 4/resources/gene_sets/hallmark_symbols.gmt`
Gene sets parsed: 50
Minimum overlap threshold: 10 genes

## Dataset Overlap

### openproblems_neurips2023
- Gene sets with >= 10 overlapping genes: 47/50
- Median overlap: 58.0
- Minimum overlap: 4
- Maximum overlap: 174

### sciplex3_24h_top2000
- Gene sets with >= 10 overlapping genes: 39/50
- Median overlap: 20.5
- Minimum overlap: 2
- Maximum overlap: 71

## Model Inclusion Notes

- Included models are limited to runs with full predicted response vectors or reconstructable train-only baseline predictions.
- Nearest-drug and SVD-ridge OpenProblems outputs were not projected because the available local files contain per-row metrics but not full prediction vectors.

## Run Notes

- OpenProblems nearest-drug and SVD-ridge were not projected because only per-row metrics, not full prediction vectors, were available locally.
