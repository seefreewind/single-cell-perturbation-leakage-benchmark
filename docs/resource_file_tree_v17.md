# v17 benchmark resource file tree and mapping

## Local resource tree

```text
benchmark_resource/
  README.md
  data_manifest.tsv
  configs/
    deep_model_panel.yaml
    model_registry.yaml
  split_manifests/
    split_manifest_openproblems.csv
    split_manifest_openproblems_summary.csv
    split_manifest_sciplex3.csv
    split_manifest_sciplex3_summary.csv
  leakage_audits/
    leakage_audit_table.csv
  model_results/
    all_model_metrics_long.csv
    all_model_leakage_audit_long.csv
    bootstrap_ci.csv
    model_rank_stability.csv
    random_to_strict_contrasts.csv
    transigen_metrics_long.csv
    transigen_predictions_manifest.csv
    transigen_leakage_audit.csv
    transigen_random_to_strict_contrasts.csv
  pathway_level/
    pathway_seed_summary.csv
    pathway_bootstrap_ci.csv
    pathway_random_to_strict_contrasts.csv
    pathway_per_gene_set_errors.csv
    pathway_gene_set_overlap.csv
  report_schema/
    benchmark_report_schema.json
    benchmark_report_schema.yaml
  source_data/
    figure5_pathway_level_recovery.csv
    figure6_sciplex3_prnet30_transigen30_model_panel.csv
    figure7_rank_transfer_prnet30.csv
    figure7_rank_transfer_prnet30_summary.csv
```

## Manuscript outputs

| Artifact | Path |
|---|---|
| v17 manuscript Markdown | `manuscript/manuscript_full_en_v17_control_conditioned_neural.md` |
| v17 tables Markdown | `manuscript/tables_draft_en_v17_control_conditioned_neural.md` |
| v17 figure legends Markdown | `manuscript/figure_legends_en_v17_control_conditioned_neural.md` |
| v17 DOCX builder | `scripts/build_manuscript_docx_v17_control_conditioned_neural.py` |
| v17 DOCX | `submission_package_v17_control_conditioned_neural/manuscript/manuscript_full_en_v17_control_conditioned_neural.docx` |

## Figure and source-data mapping

| Manuscript figure | Figure file | Source data |
|---|---|---|
| Fig. 1 | `figures/manuscript_main/figure1_benchmark_design.png` | Split manifests and dataset summaries |
| Fig. 2 | `figures/manuscript_main/figure2_baseline_decay.png` | OpenProblems baseline summaries |
| Fig. 3 | `figures/manuscript_main/figure4_chemical_similarity_audit.png` | Leakage and Tanimoto audit tables |
| Fig. 4 | `figures/figure6_sciplex3_prnet30_transigen30_model_panel.png` and `.pdf` | `results/deep_model_panel/source_data/figure6_sciplex3_prnet30_transigen30_model_panel.csv` |
| Fig. 5 | `figures/figure5_pathway_level_recovery.png` and `.pdf` | `results/deep_model_panel/source_data/figure5_pathway_level_recovery.csv` |
| Fig. 6 | `figures/figure7_rank_transfer_prnet30.png` and `.pdf` | `results/deep_model_panel/source_data/figure7_rank_transfer_prnet30.csv`; `results/deep_model_panel/source_data/figure7_rank_transfer_prnet30_summary.csv` |
| Supplementary diagnostic | `figures/figure6_transigen_ablation_sciplex3.png` and `.pdf` | `results/deep_model_panel/source_data/figure6_transigen_ablation_sciplex3.csv` |
| Supplementary diagnostic | `figures/supp_transigen_random_to_strict_ci.png` and `.pdf` | `results/deep_model_panel/source_data/supp_transigen_random_to_strict_ci.csv` |
| Supplementary diagnostic | `figures/supp_transigen_vs_ridge_paired.png` and `.pdf` | `results/deep_model_panel/source_data/supp_transigen_vs_ridge_paired.csv` |

## Public-release blockers

- Replace `[GitHub URL]` and `[Zenodo DOI]` only after a real public or anonymous repository and archive exist.
- Add final commit identifier only after repository creation and commit.
- Decide whether to include local prediction NPZ files in the repository or archive them separately because they are larger than CSV summaries.
- Add fingerprint-preprocessing metadata, especially Morgan radius, in the next resource release.
