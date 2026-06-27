# Preparing gene sets for pathway-level perturbation recovery

This project does not currently contain a local MSigDB Hallmark `.gmt` file or another approved gene-set file. Pathway-level evaluation was therefore not run in the JGG-targeted v13 revision.

## Required input

Preferred file:

- MSigDB Hallmark gene sets in GMT format, for example `h.all.v202x.x.Hs.symbols.gmt`.

Acceptable alternatives:

- Curated pathway GMT files using HGNC gene symbols.
- A project-specific CSV with columns `gene_set_name` and `gene_symbol`.

## Suggested location

Place the file under:

```text
resources/gene_sets/
```

Recommended example:

```text
resources/gene_sets/h.all.Hs.symbols.gmt
```

## Planned analysis once gene sets are available

1. Load Sci-Plex 3 target and prediction matrices for global mean, ridge cell+dose+drug FP, TranSiGen-style pseudobulk, and PRnet-style pseudobulk outputs.
2. Match gene-set symbols to the 2,000 Sci-Plex 3 genes.
3. Retain gene sets with at least 5 matched genes.
4. For each response vector and gene set, compute a signed mean response score or standardized mean response score.
5. Compute pathway-level Pearson between true and predicted pathway score vectors for each seed, split, model, and gene set.
6. Summarize random, scaffold held-out, and joint held-out performance.
7. Report random-to-strict drops and representative Hallmark groups such as stress response, apoptosis, cell cycle, interferon response, and DNA damage if present.

No gene-set file should be hard-downloaded during manuscript revision unless the user explicitly provides or approves the exact source.
