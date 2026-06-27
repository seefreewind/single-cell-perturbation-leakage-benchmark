# JGG Submission Resource Checklist

## Manuscript files

| Item | Status | Path or action |
|---|---|---|
| JGG-targeted Word manuscript | Present after build | `submission_package_v13_JGG_targeted_revision/manuscript/manuscript_full_en_v13_JGG_targeted_revision.docx` |
| JGG-targeted Markdown source | Present | `manuscript/manuscript_full_en_v13_JGG_targeted_revision.md` |
| Tables source | Present | `manuscript/tables_draft_en_v13_JGG_targeted_revision.md` |
| Figure legends source | Present | `manuscript/figure_legends_en_v13_JGG_targeted_revision.md` |
| Revision summary | Present | `docs/JGG_revision_summary.md` |

## Figure files and source data

| Figure | Status | Main file | Source data status |
|---|---|---|---|
| Figure 1 benchmark design | Present | `figures/manuscript_main/figure1_benchmark_design.png` | Existing source data should be included from manuscript/source-data folders |
| Figure 2 OpenProblems performance decay | Present | `figures/manuscript_main/figure2_baseline_decay.png` | Existing source data should be included from manuscript/source-data folders |
| Figure 3 leakage audit | Present | `figures/manuscript_main/figure4_chemical_similarity_audit.png` | Existing similarity/leakage audit tables should be included |
| Figure 4 Sci-Plex 3 model panel | Present | `figures/figure6_sciplex3_prnet30_transigen30_model_panel.png` | `results/deep_model_panel/source_data/figure6_sciplex3_prnet30_transigen30_model_panel.csv` |
| Figure 5 rank transfer | Present | `figures/figure7_rank_transfer_prnet30.png` | `results/deep_model_panel/source_data/figure7_rank_transfer_prnet30.csv`; `results/deep_model_panel/source_data/figure7_rank_transfer_prnet30_summary.csv` |
| Supplementary model-status figure | Present | `results/deep_model_panel/figures/figure8_foundation_embedding_feasibility.png` | Include corresponding model feasibility tables |
| Pathway-level figure | Not generated | Not applicable | Gene-set file missing; see `docs/pathway_analysis_feasibility_report.md` |

## Split manifests and leakage audits

| Item | Status | Path or action |
|---|---|---|
| OpenProblems split manifests | Present locally | Include relevant `results/openproblems_*/*/splits/` outputs |
| Sci-Plex 3 split manifest | Present | `data/processed/transigen/sciplex3/split_assignments_long.csv` |
| OpenProblems leakage audits | Present locally | Include leakage audit outputs used for figures and tables |
| Sci-Plex 3 PRnet leakage audit | Present | `results/deep_model_panel/sciplex3/prnet_main30_leakage_audit.csv` |
| Sci-Plex 3 TranSiGen leakage audit | Present | `results/deep_model_panel/sciplex3/transigen_main30_leakage_audit.csv` |

## Model outputs

| Item | Status | Path or action |
|---|---|---|
| Integrated long-format model metrics | Present | `results/deep_model_panel/all_model_metrics_long.csv` |
| Random-to-strict contrasts | Present | `results/deep_model_panel/random_to_strict_contrasts.csv` |
| Model-rank stability | Present | `results/deep_model_panel/model_rank_stability.csv` |
| Sci-Plex 3 completed model summary | Present | `results/deep_model_panel/sciplex3/sciplex3_completed_model_summary.csv` |
| PRnet main30 metrics | Present | `results/deep_model_panel/sciplex3/prnet_main30_metrics_long.csv` |
| PRnet main30 prediction manifest | Present | `results/deep_model_panel/sciplex3/prnet_main30_predictions_manifest.csv` |
| TranSiGen main30 metrics | Present | `results/deep_model_panel/sciplex3/transigen_metrics_long.csv` |
| TranSiGen main30 prediction manifest | Present | `results/deep_model_panel/sciplex3/transigen_main30_predictions_manifest.csv` |

## Feasibility reports and QC logs

| Item | Status | Path or action |
|---|---|---|
| Model feasibility reports | Present locally | Include `results/deep_model_panel/logs/` and related feasibility outputs |
| PRnet main30 run report | Present | `docs/prnet_main30_run_report.md` |
| PRnet main30 stdout log | Present | `logs/prnet_main30_stdout.log` |
| PRnet main30 per-run logs | Present | `logs/prnet_main30/` |
| Pathway feasibility report | Present | `docs/pathway_analysis_feasibility_report.md` |

## Environment and reproducibility

| Item | Status | Path or action |
|---|---|---|
| Existing project scripts | Present | `scripts/` |
| Existing Python environment | Present locally | `envs/transigen_venv/` |
| Environment export | Needs manual packaging | Export a reproducible `requirements.txt` or `environment.yml` before submission |
| GitHub repository | Needs manual upload | No public URL is claimed in the manuscript |
| Zenodo DOI | Needs manual archive | No DOI is claimed in the manuscript |
| Release checklist | Needs final repository structure | Create once GitHub repository is fixed |

## Human information still required

- Final author list and affiliations.
- Confirmed CRediT author-contribution roles.
- Funding statement or declaration of no funding.
- Acknowledgements.
- Competing-interest confirmation.
- Public GitHub URL.
- Zenodo DOI.
- Cover letter.
