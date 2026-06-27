# PRnet Main30 Run Report

Date: 2026-06-24

## Status

`prnet_adapted_sciplex3` main30 completed successfully for all required seeds and splits.

- Seeds: 1-30.
- Splits: random, cell-line held-out, scaffold held-out, joint cell-line-plus-scaffold held-out.
- Successful seed/split combinations: 120/120 total; 30/30 for every required split.
- Failed combinations: 0.
- Model seed matched the fixed split-manifest seed in every run.

## Performance Summary

| Split | Seeds | Mean Pearson | SD Pearson | Median | IQR | Mean RMSE |
|---|---:|---:|---:|---:|---:|---:|
| Random | 30 | 0.0175 | 0.0120 | 0.0153 | 0.0173 | 0.2324 |
| Cell-line held-out | 30 | -0.0016 | 0.0339 | -0.0028 | 0.0731 | 0.5718 |
| Scaffold held-out | 30 | 0.0151 | 0.0082 | 0.0146 | 0.0109 | 0.2327 |
| Joint held-out | 30 | 0.0006 | 0.0388 | -0.0020 | 0.0479 | 0.5769 |

Random minus scaffold Pearson difference: 0.0024 (bootstrap 95% CI -0.0016 to 0.0063).

Random minus joint Pearson difference: 0.0168 (bootstrap 95% CI 0.0043 to 0.0292).

## Paired Comparisons

Paired PRnet-minus-comparison differences are available in `results/deep_model_panel/sciplex3/prnet_main30_paired_differences.csv`. The comparisons use matched seeds and split types where available.

## Rank Transfer

Seed-level rank transfer after adding PRnet is available in `results/deep_model_panel/model_rank_stability.csv`. The random-to-joint Spearman mean for the five-model Sci-Plex 3 panel is 0.1433.

## Interpretation

This is a completed repeated-seed pseudobulk adaptation, not a full official PRnet reproduction. The run supports inclusion as a completed Sci-Plex 3 model-panel result, but the low Pearson values indicate that this adapter is not competitive with the ridge or control-conditioned repeated-seed results in the current configuration.
