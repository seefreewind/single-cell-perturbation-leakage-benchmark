# Pathway-Level Analysis Feasibility Report

## Status

Superseded by the completed v14/v15 pathway-level analysis.

## Current state

The earlier v13 feasibility limitation was resolved after the user provided a local MSigDB Hallmark GMT file. The pathway-level recovery analysis was completed without downloading gene sets online.

Current pathway-analysis outputs are:

- `docs/pathway_gene_set_check.md`
- `results/pathway_level/pathway_gene_set_overlap.csv`
- `results/pathway_level/pathway_metrics_long.csv`
- `results/pathway_level/pathway_seed_summary.csv`
- `results/pathway_level/pathway_random_to_strict_contrasts.csv`
- `results/pathway_level/pathway_per_gene_set_errors.csv`
- `results/pathway_level/pathway_bootstrap_ci.csv`
- `figures/figure5_pathway_level_recovery.png`
- `figures/figure5_pathway_level_recovery.pdf`
- `figures/source_data/figure5_pathway_level_recovery.csv`

## Public release note

The MSigDB Hallmark GMT file itself must not be redistributed in the public benchmark resource package. Users should obtain Hallmark gene sets directly from MSigDB according to its terms of use and place the file at `resources/gene_sets/hallmark_symbols.gmt`.
