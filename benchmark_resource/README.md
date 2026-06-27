# Leakage-Aware Single-Cell Chemical Perturbation Benchmark Resource

This resource package accompanies the manuscript draft:

**Scaffold-strict benchmarking reveals chemical-neighbor shortcuts and model-ranking instability across modern single-cell drug perturbation models**

The package provides reusable split manifests, leakage-audit tables, model-result summaries, report schemas, and feasibility documentation for evaluating single-cell drug perturbation prediction models under leakage-aware validation settings.

## Data Sources

- OpenProblems NeurIPS 2023 single-cell perturbation prediction response matrices.
- Sci-Plex 3 24 h drug-cell-line-dose pseudobulk response matrix derived from the local scPerturb/Sci-Plex 3 AnnData file.
- Project-curated drug metadata, SMILES mappings, Bemis-Murcko scaffolds, MoA annotations, and dose metadata.

See `data_manifest.tsv` for local paths and source notes.

## Split Manifests

The package contains fixed split manifests for:

- OpenProblems: random, cell held-out, scaffold held-out, joint cell+scaffold held-out, and MoA-held-out secondary split.
- Sci-Plex 3: random, cell-line held-out, scaffold held-out, and joint cell-line+scaffold held-out.

All model adapters must reuse these manifests. Test assignments must not be changed for model convenience.

## Leakage Audits

Leakage-audit outputs report:

- same-drug overlap;
- same-scaffold overlap;
- same-cell or same-cell-line overlap;
- same-MoA overlap where available;
- nearest training-drug Tanimoto similarity;
- Sci-Plex 3 dose leakage, including same drug-dose and same drug-cell-line-dose overlap.

## Model Results

`model_results/` contains standardized long-format outputs for completed existing baselines and explicit status rows for planned modern/foundation models that were not runnable in the current local environment.

Core files:

```text
all_model_metrics_long.csv
all_model_leakage_audit_long.csv
model_rank_stability.csv
random_to_strict_contrasts.csv
bootstrap_ci.csv
```

## Running the Current Resource Assembly

From the project root:

```bash
.venv/bin/python scripts/run_deep_model_panel_smoke.py
.venv/bin/python scripts/qc_deep_model_panel.py
.venv/bin/python scripts/plot_deep_model_panel_figures.py
```

## Running the Deep Model Panel

The adapter interface is defined in `src/models/base.py`. Planned adapters are located under:

```text
src/models/deep/
src/models/foundation/
```

Deep perturbation and foundation-model arms require external model repositories, compatible `torch` runtimes, and model-specific input builders. Failures must be recorded in `model_feasibility_report.md` and in long-format status rows.

## Report Schema

`report_schema/` contains JSON and YAML schemas for reporting benchmark metadata, split leakage, metric outputs, and model feasibility status.

## Citation and DOI

Code and benchmark resources will be deposited before submission at `[GitHub URL]` and archived with Zenodo at `[Zenodo DOI]`. Do not replace these placeholders until public release records exist.
