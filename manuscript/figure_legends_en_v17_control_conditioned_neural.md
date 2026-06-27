# Figure legends for v17 control-conditioned neural manuscript

## Figure 1. Leakage-aware split framework for single-cell perturbation transcriptomics.

Benchmark schematic and split coverage for random, cell held-out, scaffold held-out, and joint cell-plus-scaffold held-out validation. The figure summarizes the escalation from local interpolation to stricter chemical and cellular extrapolation tasks.

## Figure 2. Random validation overestimates transcriptomic perturbation prediction.

Completed OpenProblems baseline performance across split types. Pearson correlations are shown for all-gene and DE-gene response targets where available. Error bars summarize variability across split seeds.

## Figure 3. Chemical-neighbor leakage audit across validation splits.

Nearest training-drug Tanimoto similarity and paired split contrasts. Random splits retain close chemical neighbors, whereas scaffold held-out and joint held-out splits remove same-scaffold shortcuts. Same-drug, same-scaffold, and same-MoA audits are reported in the corresponding source tables.

## Figure 4. Sci-Plex 3 external validation and neural pseudobulk adapters.

Repeated-seed Sci-Plex 3 model panel for global mean, ridge drug-fingerprint, ridge cell-plus-dose-plus-drug-fingerprint, the control-conditioned response regressor, and the PRnet-interface pseudobulk adapter. Bars summarize all-gene row-wise Pearson under random, scaffold held-out, and joint cell-line-plus-scaffold held-out validation. The neural arms are benchmark-compatible pseudobulk adaptations rather than full original workflow reproductions.

## Figure 5. Pathway-level recovery also deteriorates under strict extrapolation.

Hallmark gene-set-level response recovery for completed vector-level predictors. Panels summarize pathway-level Pearson in OpenProblems and Sci-Plex 3, paired random-to-joint pathway recovery drops, and representative Hallmark gene-set absolute errors under joint cell-line-plus-scaffold validation. Gene sets were retained only when at least 10 member genes overlapped the evaluated response matrix.

## Figure 6. Random validation does not reliably preserve model ranking under strict extrapolation.

Sci-Plex 3 seed-level rank-transfer analysis after adding the 30-seed PRnet-interface pseudobulk adapter. Ranks are computed within each seed and split for the five-model repeated-seed panel. The right panel summarizes random-to-strict Spearman rank transfer across seeds.
