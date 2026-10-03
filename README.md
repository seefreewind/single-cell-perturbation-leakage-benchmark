# Auditing chemical and cellular train-test overlap in perturbation-response prediction

## Current fixed-cohort revision

The current manuscript is **Auditing chemical and cellular train-test overlap in perturbation-response prediction**. Current revision code, frozen outputs, consolidated Supplementary Information, logs, and exact file checksums are described in [BMC_revision_20261003/README.md](BMC_revision_20261003/README.md). This revision evaluates an audit protocol with transparent probes; it does not rank untested model families or reproduce original published neural workflows. The revision resources are archived as [Zenodo version 1](https://doi.org/10.5281/zenodo.23116924), corresponding to GitHub release `1` at commit `c644e76c6139b27313542434bf44ed22c8f53919`. Later documentation-only commits do not change that immutable archive. The sections below document the preserved historical v17 release, separately archived under DOI `10.5281/zenodo.21348707`.

Current clean submission documents and Figure 1-4 are indexed in [BMC_Bioinformatics_submission/README.md](BMC_Bioinformatics_submission/README.md). Figure-specific numerical source data and the immutable analysis release remain available under the revision archive above.

## Historical release / archived previous version

The following sections describe the preserved v17 release. Its manuscript, figures, and supplementary material are historical and should not be interpreted as the current submission.

The historical manuscript draft is:

**Leakage-aware benchmarking reveals chemical-neighbor leakage and model-ranking instability in single-cell chemical perturbation transcriptomics**

## What is included

- Fixed split manifests for OpenProblems NeurIPS 2023 and Sci-Plex 3 pseudobulk benchmark records.
- Leakage-audit outputs for same-drug, same-scaffold, same-cell-context, mechanism-neighbor, dose-neighbor, and nearest-chemical-neighbor checks.
- Long-format model metrics and random-to-strict contrasts for completed baselines and benchmark-compatible neural adapters.
- Source data for main and supplementary figures.
- Scripts used to build the benchmark outputs and manuscript figures.
- Historical manuscript draft, tables, figure legends, rendered DOCX, and rendered PDF.
- QC and provenance reports for the v17 revision.

## What is not included

Large raw input datasets are not redistributed in this GitHub repository. Users should obtain OpenProblems, scPerturb/Sci-Plex, and MSigDB resources from their original providers and according to their terms of use. The repository includes the current processed Sci-Plex 3 pseudobulk matrices and fixed split manifests used for this draft, but not the raw `.h5ad` source files.

The full per-record pathway metric table `pathway_metrics_long.csv`, neural prediction-vector NPZ files, and large TIFF exports are excluded from GitHub because of size. Summary-level pathway outputs, figure source data, processed Sci-Plex 3 pseudobulk matrices, and rendered PNG/PDF figures are included. A fuller archive bundle can be uploaded to Zenodo from the prepared `zenodo_upload/` directory in the working project.

## Key directories

```text
benchmark_resource/          Reusable benchmark resource package
configs/                     Model-panel and adapter configuration
data/processed/              Processed Sci-Plex 3 pseudobulk matrices and split manifest
docs/                        QC, provenance, submission, and run reports
env/                         Requirements file
figures/                     Main/supplementary figures and source data
manuscript/                  Historical v17 manuscript markdown, tables, and legends
metadata/                    Drug, perturbation, and split metadata
results/                     Model metrics, leakage audits, pathway summaries, and source tables
scripts/                     Analysis, QC, plotting, and manuscript-build scripts
src/                         Reusable model and metric code
submission_package_v17_control_conditioned_neural/  Historical DOCX/PDF manuscript draft
tests/                       Split, prediction-alignment, and leakage tests
```

## Historical manuscript artifacts

- `manuscript/manuscript_full_en_v17_control_conditioned_neural.md`
- `manuscript/tables_draft_en_v17_control_conditioned_neural.md`
- `manuscript/figure_legends_en_v17_control_conditioned_neural.md`
- `submission_package_v17_control_conditioned_neural/manuscript/manuscript_full_en_v17_control_conditioned_neural.docx`
- `submission_package_v17_control_conditioned_neural/manuscript/manuscript_full_en_v17_control_conditioned_neural.pdf`

## Reproducibility notes

The main QC checks used in the historical v17 revision were:

```bash
python scripts/qc_deep_model_panel.py
python tests/test_transigen_input_shapes.py
python tests/test_transigen_no_test_leakage.py
python tests/test_transigen_prediction_alignment.py
python tests/test_transigen_split_integrity.py
python tests/test_transigen_input_ablation_shapes.py
python tests/test_prnet_main30_prediction_alignment.py
python tests/test_prnet_main30_split_integrity.py
```

The most important provenance files are:

- `docs/manuscript_number_provenance_v17.md`
- `docs/qc_self_audit_report_v17.md`
- `docs/resource_file_tree_v17.md`
- `summary_v17_control_conditioned_neural_revision.md`

## Historical citation

The historical release is separately archived at:

https://doi.org/10.5281/zenodo.21348707

## License and data-use terms

See `LICENSE.md`. Code and generated summaries are provided for research reuse, but upstream datasets and MSigDB gene sets remain governed by their original licenses and terms of use.
