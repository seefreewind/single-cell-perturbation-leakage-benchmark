# v17 control-conditioned neural revision summary

## Completed

- Created the v17 English manuscript draft:
  - `manuscript/manuscript_full_en_v17_control_conditioned_neural.md`
  - `manuscript/tables_draft_en_v17_control_conditioned_neural.md`
  - `manuscript/figure_legends_en_v17_control_conditioned_neural.md`
  - `submission_package_v17_control_conditioned_neural/manuscript/manuscript_full_en_v17_control_conditioned_neural.docx`
- Reframed the completed Sci-Plex 3 neural arm as a **control-conditioned response regressor** rather than a TranSiGen reproduction.
- Added a Methods-ready architecture/specification table as Table 7.
- Updated figure labels and source data for the Sci-Plex 3 model panel, rank-transfer panel, and pathway-level panel.
- Added a diagnostic 10-seed input/fingerprint ablation paragraph, explicitly marked as diagnostic rather than a 30-seed main result.
- Promoted model-ranking instability in the Abstract and Results while preserving the leakage-audit/resource framing.
- Created audit and resource documentation:
  - `docs/manuscript_number_provenance_v17.md`
  - `docs/qc_self_audit_report_v17.md`
  - `docs/resource_file_tree_v17.md`
  - `docs/v16_to_v17_manuscript.diff`

## Key neural results now reported

| Result | Value |
|---|---|
| Control-conditioned regressor random Pearson | 0.233 |
| Control-conditioned regressor scaffold Pearson | 0.212 |
| Control-conditioned regressor joint Pearson | 0.017 |
| Random minus scaffold paired difference | 0.021, 95% CI 0.013 to 0.029 |
| Random minus joint paired difference | 0.216, 95% CI 0.208 to 0.224 |
| Neural minus ridge, random | -0.059, 95% CI -0.062 to -0.056 |
| Neural minus ridge, scaffold | +0.023, 95% CI 0.019 to 0.027 |
| Neural minus ridge, joint | -0.043, 95% CI -0.054 to -0.032 |
| Drug-only diagnostic ablation random-minus-scaffold | 0.067, 95% CI 0.051 to 0.084 |
| Full diagnostic ablation random-minus-scaffold | 0.019, 95% CI 0.005 to 0.033 |

## QC status

- `scripts/qc_deep_model_panel.py`: passed.
- Control-conditioned adapter tests: input shape, no-test-leakage, prediction alignment, split integrity, and input-ablation shape tests passed.
- PRnet-interface main30 tests: prediction alignment and split integrity tests passed.
- Scaffold held-out same-drug and same-scaffold overlap were 0 for the completed Sci-Plex 3 neural outputs.
- Joint held-out same-drug, same-scaffold, and same-cell-line overlap were 0.
- Prediction NPZ files checked from the completed manifests had 2,000-gene prediction vectors and matched target shapes.
- Reference order check passed: references first appear in order 1-40, 40 references are listed, and no `[REF needed]` marker remains.
- DOCX render QA completed: 16 rendered pages, with Figure 5 and Figure 6 pages visually checked after label updates.

## Important non-fabrication notes

- No GitHub URL, commit ID, or Zenodo DOI was added because the current workspace is not a git repository and no real archive was created.
- Sci-Plex 3 DE-gene-only neural metrics were not added because the current pseudobulk artifact does not include aligned per-record DE masks.
- Training and validation loss curves were not added because they were not persisted during the completed 30-seed neural runs.
- Morgan fingerprint length is recorded as 1,024 bits, but the radius is not stored in the current output metadata.

## Submission readiness

The manuscript is stronger and more internally consistent after this revision, especially around the neural arm and source-data labels. Remaining submission blockers are administrative/resource-release items rather than manuscript-logic blockers:

- Public or anonymous GitHub repository.
- Zenodo DOI or equivalent archive.
- Final supplementary file bundle.
- Final author contribution confirmation.
- Funding statement confirmation.
- Acknowledgements.
- Cover letter and journal-specific formatting.
