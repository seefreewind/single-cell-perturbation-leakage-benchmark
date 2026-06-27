# PRnet Main10 Run Report

Date: 2026-06-24

## Status

`prnet_adapted_sciplex3` main10 completed successfully for all required seeds and splits.

- Seeds: 1-10.
- Splits: random, scaffold held-out, joint cell-line-plus-scaffold held-out.
- Successful seed/split combinations: 10/10 for every required split.
- Failed combinations: 0.
- Model seed matched the fixed split-manifest seed in every run.

## Performance Summary

| Split | Seeds | Mean Pearson | SD Pearson | Mean RMSE |
|---|---:|---:|---:|---:|
| Random | 10 | 0.0165 | 0.0095 | 0.2325 |
| Scaffold held-out | 10 | 0.0144 | 0.0089 | 0.2259 |
| Joint held-out | 10 | 0.0013 | 0.0488 | 0.5822 |

Random minus scaffold Pearson difference: 0.0020.

Random minus joint Pearson difference: 0.0152.

## Paired Comparisons

Paired PRnet-minus-comparison differences are available in `results/deep_model_panel/sciplex3/prnet_main10_paired_differences.csv`. The comparisons use matched seeds and split types where available.

## Rank Transfer

Seed-level rank transfer after adding PRnet is available in `results/deep_model_panel/model_rank_stability.csv`. The random-to-joint Spearman mean for the five-model Sci-Plex 3 panel is 0.1100.

## Interpretation

This is a completed repeated-seed pseudobulk adaptation, not a full official PRnet reproduction. The run supports inclusion as a completed Sci-Plex 3 model-panel result, but the low Pearson values indicate that this adapter is not competitive with the ridge or TranSiGen-style repeated-seed results in the current configuration.
