# JGG-Targeted Revision Summary

## Manuscript repositioning

The manuscript was revised from a computational drug-model benchmark toward a Journal of Genetics and Genomics-oriented single-cell functional genomics and perturbation transcriptomics benchmark resource.

Major positioning changes:

- Title changed to: "Leakage-aware benchmarking exposes shortcut learning in single-cell chemical perturbation transcriptomics".
- Abstract now emphasizes single-cell perturbation transcriptomics, functional-genomics benchmark validity, leakage-aware split design, and resource outputs.
- Model-name stacking was reduced in the abstract.
- OpenProblems, Sci-Plex 3, MoA-held-out validation, TranSiGen-style pseudobulk adaptation, and PRnet-style pseudobulk adaptation remain in the manuscript.
- chemCPA/CPA, scGPT, and scFoundation remain feasibility/status rows and are not described as completed predictors.
- The Results and Discussion now foreground benchmark validity, mechanism-aware validation, functional-genomics interpretation, and reporting recommendations.

## Author and correspondence information

The manuscript title page was updated with:

- First author: Da Lin.
- Corresponding author: Yu Zhang.
- Corresponding email: zhangyu1@wzhealth.com.
- Affiliation: The Second Affiliated Hospital of Wenzhou Medical University, Wenzhou, Zhejiang, China.

## Results restructuring

The Results section was reorganized into the requested JGG-facing structure:

1. Leakage-aware split construction for single-cell perturbation transcriptomics.
2. Random validation overestimates transcriptomic perturbation prediction.
3. Chemical-neighbor shortcuts account for random-split advantage.
4. Mechanism-held-out validation separates scaffold independence from MoA-level extrapolation.
5. Sci-Plex 3 external validation reproduces random-to-strict performance decay.
6. Neural pseudobulk adaptations remain vulnerable under joint extrapolation.
7. Random validation does not preserve model ranking.

The MoA-held-out analysis was promoted from an auxiliary point to a main Results subsection.

## Discussion restructuring

The Discussion was rewritten with the requested subheadings:

1. Principal findings.
2. Random validation as local interpolation rather than extrapolation.
3. Scaffold independence is not mechanism independence.
4. Why model feasibility is part of benchmark validity.
5. Implications for single-cell foundation model evaluation.
6. Practical reporting recommendations for perturbation transcriptomics benchmarks.
7. Limitations.

The revised Discussion explicitly distinguishes same-drug, same-scaffold, and same-cell-context shortcuts; frames MoA-held-out validation as a mechanism-level secondary stress test; and states that pseudobulk adaptations do not represent the final capability of the original TranSiGen or PRnet workflows.

## Pathway-level biological evaluation

Pathway-level analysis was not added because no local MSigDB Hallmark `.gmt` file or other approved gene-set file was present in the workspace. The instructions prohibited hard-downloading gene sets during this revision.

Generated instead:

- `scripts/prepare_gene_sets.md`
- `docs/pathway_analysis_feasibility_report.md`

The v13 manuscript does not claim pathway-level results. The limitation section states why pathway-level recovery was not included.

## Model-result integrity

No original model result values were changed. Existing OpenProblems, Sci-Plex 3, TranSiGen-style, and PRnet-style values were preserved.

The manuscript continues to state:

- TranSiGen-style and PRnet-style results are pseudobulk adaptations.
- They are not full reproductions of the original model workflows.
- The PRnet-style pseudobulk adaptation was executable but did not improve strict extrapolation under this benchmark-compatible adaptation.
- This should not be interpreted as evidence that the original PRnet workflow fails.

## Figure and table changes

Main figures were reorganized as:

1. Benchmark design and leakage-aware split framework.
2. OpenProblems random-to-strict performance decay.
3. Leakage audit and chemical-neighbor similarity.
4. Sci-Plex 3 external validation and neural pseudobulk adaptations.
5. Model ranking instability.

The previous model-status figure was moved out of the main figure flow and represented in the text through Table 4 and the resource checklist.

## New supporting documents

Generated for this revision:

- `docs/JGG_revision_summary.md`
- `docs/JGG_submission_resource_checklist.md`
- `docs/references_to_add_JGG.md`
- `docs/pathway_analysis_feasibility_report.md`
- `scripts/prepare_gene_sets.md`

## Manuscript output

Generated:

- `submission_package_v13_JGG_targeted_revision/manuscript/manuscript_full_en_v13_JGG_targeted_revision.docx`
- `submission_package_v13_JGG_targeted_revision/rendered_docx/manuscript_full_en_v13_JGG_targeted_revision.pdf`

The DOCX was rendered to page images and visually checked for the title page, main figures, Table 4, Figure 5, and declaration sections.

## Remaining human tasks before submission

- Replace `[REF needed]` markers with verified references.
- Complete and verify the expanded reference list.
- Confirm final author list and affiliations.
- Fill exact CRediT author contributions.
- Fill acknowledgements.
- Fill funding statement or declare no funding.
- Confirm competing-interest statement.
- Create or provide public GitHub repository URL.
- Create or provide Zenodo DOI.
- Package supplementary files and source data.
- Prepare cover letter for JGG.
