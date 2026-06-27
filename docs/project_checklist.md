# Project Checklist

Last updated for manuscript v7 draft: 2026-06-24.

## Phase 1: feasibility audit

- [x] Identify one primary perturbation dataset for the minimum viable benchmark.
- [x] Collect or construct drug metadata with `drug_id`, `drug_name`, and `smiles`.
- [x] Collect perturbation metadata with `sample_id`, `drug_id`, `cell_context`, `batch`, `dataset`, `dose`, and `time`.
- [x] Compute Bemis-Murcko scaffolds from SMILES.
- [x] Count drugs, scaffolds, cell contexts, batches, and drug-cell combinations.
- [x] Count scaffolds measured in at least two cell contexts.
- [ ] Estimate overlap with PRISM, GDSC, CTRP, or DepMap drug-response resources.
- [x] Generate random, cell-held-out, scaffold-held-out, and cell-plus-scaffold-held-out split candidates.
- [x] Run one simple baseline before running foundation models.
- [ ] Run a small embedding extraction smoke test for one foundation model.

## Phase 2: minimum viable benchmark

- [x] Fit mean baselines.
- [ ] Fit PCA + ridge regression.
- [x] Fit Morgan fingerprint + ridge regression.
- [ ] Fit basal expression + drug fingerprint fusion baseline.
- [ ] Extract embeddings from one or two single-cell foundation models.
- [x] Evaluate Level 1, Level 3, Level 4, and Level 5 splits.
- [x] Report DE-gene correlation and top-k DE-gene overlap.

## Phase 3: leakage audit

- [ ] Train classifiers to predict cell identity from embeddings.
- [ ] Train classifiers to predict batch and dataset source from embeddings.
- [x] Compute test-drug nearest-neighbor similarity to training drugs.
- [x] Relate structural similarity to prediction error.
- [x] Compare ordinary drug split against scaffold-strict split.
- [x] Add same-drug, same-scaffold, and same-cell overlap audit.
- [x] Add 100-seed OpenProblems statistical hardening.
- [x] Add ridge alpha and fingerprint-length sensitivity audit.
- [x] Add Sci-Plex 3 independent-dataset smoke test.
- [x] Upgrade Sci-Plex 3 to 30-seed external validation.
- [x] Add cross-dataset model-ranking instability figure.
- [x] Add minimum leakage-audit reporting checklist.
- [x] Add paired bootstrap confidence intervals for OpenProblems and Sci-Plex 3 random-to-strict performance deltas.
- [x] Report seed-level uncertainty and negative-seed fraction for Sci-Plex 3 model-rank instability.
- [x] Report row-wise Tanimoto-performance correlations or explicitly mark near-constant settings as unstable.
- [x] Add size-matched and composition-matched random controls for OpenProblems joint held-out.
- [x] Add OpenProblems same-MoA/mechanism-neighbor overlap audit.
- [x] Add Sci-Plex 3 same-drug-dose and dose-grid leakage audit.
- [x] Add OpenProblems MoA-held-out secondary split performance audit.
- [x] Generate reusable split manifests for OpenProblems and Sci-Plex 3.
- [x] Generate local leakage audit table, template, and benchmark report schema.
- [x] Add MoA-held-out performance tables to the benchmark resource package.

## Phase 4: manuscript build

- [ ] Freeze the benchmark design before large-scale model runs.
- [x] Save exact split definitions.
- [x] Save software versions and parameters.
- [x] Draft Results in workflow order, not in order of model preference.
- [x] Keep Results descriptive and reserve interpretation for Discussion.
- [x] Include negative findings and failed validations.
- [x] Correct Figure 5/Figure 6 first-citation order.
- [x] Mark mixed seed counts in Table 3 and Figure 2 legend.
- [x] Align Figure 4 legend with the actual two-panel figure.
- [x] Add discussion of OpenProblems versus Sci-Plex 3 cell-held-out/scaffold-held-out difficulty reversal.
- [x] Add limitation for total-expression-based top-2,000 gene selection in Sci-Plex 3.
- [x] Reframe title and abstract around chemical-neighbor shortcuts and model mis-ranking.
- [x] Reorder Table 6-8 to match first citation order.
- [x] Reframe manuscript as a benchmark resource package rather than only an audit analysis.
- [x] Create v8 English deep-model-panel manuscript draft.
- [x] Add unified model adapter interface and model registry.
- [x] Add explicit feasibility/status reporting for PRnet, TranSiGen, chemCPA/CPA, and frozen foundation-model arms.
- [x] Add standardized deep-model-panel long-format output tables.
- [x] Add v8 resource-package README, data manifest, model feasibility report, and reproducibility guide.
- [x] Clone and document the TranSiGen repository without modifying third-party source.
- [x] Build Sci-Plex 3 TranSiGen pseudobulk input matrices with basal/control expression, response targets, compound metadata, and fixed split assignments.
- [x] Run TranSiGen-style Sci-Plex 3 smoke, pilot, 10-seed, and 30-seed stress tests.
- [x] Integrate TranSiGen metrics, predictions manifest, leakage audit, contrasts, rank transfer, figures, and v9 manuscript draft.

## Remaining before submission-grade manuscript

- [x] Run 30 repeated Sci-Plex 3 split seeds.
- [x] Add at least one true perturbation/foundation-model baseline or embedding-based baseline.
- [ ] Add Tanimoto-threshold or cluster-held-out split as a supplemental sensitivity analysis.
- [x] Add main figure for model-ranking instability across OpenProblems and Sci-Plex 3.
- [x] Convert Table 7/8 content into a reusable local benchmark resource with split manifests, leakage-audit template, and report schema.
- [ ] Upload resource package to public GitHub and archive a release DOI.
- [ ] Add forest plot summary of random-to-strict effect sizes across datasets and metrics.
- [ ] Update earlier main figures to fully reflect the current expanded 100-seed and Sci-Plex 3 results.
- [ ] Replace placeholder Acknowledgements, Authors' contributions, Funding, repository URL, and data DOI fields.
- [ ] Decide target journal and final manuscript language.
