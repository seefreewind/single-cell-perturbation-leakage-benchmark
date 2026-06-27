# TranSiGen Sci-Plex 3 Input Package

Source response matrix: `/Users/zy/Documents/New project 4/data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad`

Records: 2200
Genes: 2000

Matrices in `matrices.npz`:
- `basal`: matched cell-line 24 h control mean log1p(CP10K), shape (2200, 2000)
- `treated`: treated mean log1p(CP10K), shape (2200, 2000)
- `response`: treated minus matched control response target, shape (2200, 2000)

Compound fingerprints:
- `compound_fingerprints.npz`, Morgan radius 2, 1024 bits.

Split assignments:
- `split_assignments_long.csv` is copied from the existing Sci-Plex 3 split manifest.
- Excluded records must not enter training or testing.
- Test records must not be used for normalization, early stopping, or hyperparameter selection.

Adaptation note:
This package supports `transigen_adapted_sciplex3`, a pseudobulk adaptation that uses basal/control expression, compound fingerprints, dose, and cell-line features to predict treated-control response vectors.
