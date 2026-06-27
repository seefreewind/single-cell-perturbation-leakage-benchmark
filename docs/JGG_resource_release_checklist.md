# JGG Resource Release Checklist

| Required item | Status | Notes |
|---|---|---|
| `README.md` | present | `/Users/zy/Documents/New project 4/benchmark_resource/README.md` |
| `data_manifest.tsv` | present | `/Users/zy/Documents/New project 4/benchmark_resource/data_manifest.tsv` |
| `split_manifests` | present | `/Users/zy/Documents/New project 4/benchmark_resource/split_manifests` |
| `leakage_audits` | present | `/Users/zy/Documents/New project 4/benchmark_resource/leakage_audits` |
| `model_results` | present | `/Users/zy/Documents/New project 4/benchmark_resource/model_results` |
| `pathway_level` | present | `/Users/zy/Documents/New project 4/benchmark_resource/pathway_level` |
| `source_data` | present | `/Users/zy/Documents/New project 4/benchmark_resource/source_data` |
| `report_schema` | present | `/Users/zy/Documents/New project 4/benchmark_resource/report_schema` |
| `configs` | present | `/Users/zy/Documents/New project 4/benchmark_resource/configs` |
| `model_feasibility_report.md` | present | `/Users/zy/Documents/New project 4/benchmark_resource/model_feasibility_report.md` |
| `QC_logs` | present | `/Users/zy/Documents/New project 4/benchmark_resource/QC_logs` |
| `environment.yml` or `requirements.txt` | needs manual addition | Add before public release if not present. |

## GitHub release checklist

- Add final code, split manifests, leakage audits, result summaries, source data, schemas, configs, and QC reports.
- Do not include `resources/gene_sets/hallmark_symbols.gmt` or any redistributed MSigDB GMT file.
- Add a license only after the authors confirm the intended license.
- Tag a release and record the release version in the manuscript before submission.

## Zenodo archive checklist

- Archive the GitHub release after GitHub contents are finalized.
- Record the Zenodo DOI in the manuscript only after it exists.
- Exclude MSigDB GMT files from the archive.
- Include README, manifest, source data, and schema files in the archive.
