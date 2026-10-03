# Fixed-cohort revision reproducibility files

This revision accompanies **Auditing chemical and cellular train-test overlap in perturbation-response prediction**. It supplies an evaluation protocol and cohort-specific scripts, not a general plug-in perturbation-model framework. Published-model adapters are not full reproductions of their original workflows. PRnet results are not primary evidence.

## Frozen outputs

`results/bmc_master/` contains the fixed-manifest master tables. This directory adds paired baseline-adjusted scores, common-test-record comparisons, chemical annotation coverage, and training-only top-2,000-gene sensitivity outputs and actual run logs. No analysis was rerun during public-release preparation.

Intervals are 95% split-seed bootstrap intervals, conditional on the observed cohorts. They are not independent biological-sample uncertainty estimates. OpenProblems uses fixed seeds 1-100; Sci-Plex 3 uses fixed seeds 1-30. The common-test-record comparison excludes empty intersections.

## Reproduction

Install the audit CLI from the repository root with `python -m pip install .`, then run `perturbation-audit --help`. The utility supports record-level audits and test-size-matched manifests. Train-and-test-size controls and composition controls are separate cohort-specific scripts.

```sh
python -m pytest tests/test_audit_benchmark_cli.py
python BMC_revision_20261003/analyse_revision.py
python BMC_revision_20261003/train_only_gene_sensitivity.py
python scripts/plot_bmc_revision_figures.py
python scripts/plot_bmc_final_framework.py
```

The second command regenerates summaries from immutable per-record outputs. The third reruns the existing training-only sensitivity, requiring the original full-gene Sci-Plex H5AD at `data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad`. Obtain it from the original scPerturb source; it is not redistributed here. Existing processed targets, metadata and fixed manifests are included. Do not substitute newly constructed split assignments.

The original cell-line control profiles remain part of Sci-Plex response construction. The nonlinear probe also uses the held-out line's measured basal profile; this is not entirely unseen-cellular zero-shot prediction.

## Provenance and terms

`release_sha256.json` records exact packaged bytes. `analysis_provenance.json` describes the actual revision postprocessing environment; `envs/transigen_env.yml` describes the archived neural training environment. They represent different stages, not one reconstructed training environment. `supplementary_BMC_final_polished.docx` contains the consolidated supplementary methods and Tables S1-S4.

MoA source: https://openproblems-bio.s3.amazonaws.com/public/neurips-2023-competition/moa_annotations.csv

Expected SHA256: `2ab652b6ec33239b5e4172176fb8a68980e7066b0fae9fd54a12efc8f920d30e`.

No raw source H5AD files, MSigDB Hallmark gene-set redistribution, private reviews, author correspondence, credentials, or local runtime installations are included. Upstream data remain governed by their original access and reuse terms. This revision is archived as [Zenodo version 1](https://doi.org/10.5281/zenodo.23116924), corresponding to GitHub release `1` at commit `c644e76c6139b27313542434bf44ed22c8f53919`. The checksums describe the frozen packaged files at that commit; subsequent documentation-only DOI updates do not alter the archived results. The earlier DOI `10.5281/zenodo.21348707` remains a separate historical archive.
