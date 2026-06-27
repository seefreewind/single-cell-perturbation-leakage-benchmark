# Leakage-aware benchmarking reveals chemical-neighbor leakage and model-ranking instability in single-cell chemical perturbation transcriptomics

Da Lin; Yu Zhang*

Affiliation: The Second Affiliated Hospital of Wenzhou Medical University, Wenzhou, Zhejiang, China

*Correspondence: Yu Zhang, zhangyu1@wzhealth.com

## Abstract

### Background

Single-cell perturbation transcriptomics provides a functional-genomics readout of how cells respond to drugs, genetic state, and cellular context [1-7]. These data are increasingly used to benchmark response-prediction models, but standard random validation can leave the same drugs, related chemical scaffolds, mechanisms of action, or cellular contexts in both training and test records. Such overlap can turn an intended extrapolation benchmark into local interpolation and can distort model selection.

### Methods

We developed a leakage-aware benchmark audit for single-cell chemical perturbation transcriptomics using OpenProblems NeurIPS 2023 perturbation response records as the primary benchmark and Sci-Plex 3 24 h drug-cell-line-dose pseudobulk responses as an external validation cohort. We compared random, cell or cell-line held-out, scaffold held-out, and joint cell-plus-scaffold held-out splits, and we added an OpenProblems mechanism-of-action (MoA)-held-out stress test. The audit quantified chemical-neighbor leakage, MoA-neighbor leakage, model-rank transfer, Hallmark gene-set-level response recovery, and model feasibility. Completed predictors included mean, ridge, nearest-drug, low-rank ridge, MoA-held-out ridge, a control-conditioned neural response regressor, and a PRnet-interface pseudobulk adapter in Sci-Plex 3. chemCPA/CPA, scGPT, and scFoundation were retained as feasibility/status rows rather than completed predictors.

### Results

In OpenProblems, ridge cell-plus-drug fingerprint performance decreased from 0.362 under random validation to 0.234 under scaffold held-out validation and 0.118 under joint cell-plus-scaffold held-out validation. DE-gene Pearson values decreased from 0.568 to 0.392 and 0.209. Random splits retained same-drug and same-scaffold training neighbors for 98.9% and 97.5% of test records, whereas scaffold-strict and joint-strict splits removed these chemical-identity shortcuts. Model-rank transfer was weak: random-to-joint Spearman correlation was 0.060 in OpenProblems and 0.143 in the five-model Sci-Plex 3 repeated-seed panel. Sci-Plex 3 reproduced the random-to-strict decay for ridge cell-plus-dose-plus-drug-fingerprint models, with Pearson values of 0.292, 0.189, and 0.060 under random, scaffold held-out, and joint held-out validation. The control-conditioned neural response regressor showed random and scaffold signal but not joint extrapolation, with all-gene Pearson of 0.233, 0.212, and 0.017. Hallmark pathway-level recovery showed the same pattern for vector-level predictors with clear random-split signal: Sci-Plex 3 ridge cell-plus-dose-plus-drug-fingerprint pathway Pearson decreased from 0.313 under random validation to 0.070 under joint held-out validation.

### Conclusions

Random validation can overestimate transcriptomic perturbation prediction and can mis-rank models relative to scaffold-strict and joint cellular-chemical extrapolation. The resource provides fixed split manifests, leakage audits, model-output schemas, feasibility/status rows, and quality-control checks for perturbation transcriptomics benchmarks. These results support routine reporting of chemical-neighbor leakage, mechanism-neighbor leakage, model adaptation level, and random-to-strict rank transfer before random-split gains are interpreted as evidence of biological or chemical generalization.

## Keywords

single-cell perturbation transcriptomics; functional genomics; chemical perturbation; scaffold split; benchmark leakage; mechanism of action; model ranking; Sci-Plex; OpenProblems

## Introduction

Single-cell perturbation assays make it possible to connect chemical or genetic interventions to transcriptome-wide cellular responses [1-4]. In drug perturbation studies, the measured response can be interpreted as a functional-genomics profile of pathway activity, stress response, lineage state, or mechanism-specific transcriptional change, extending earlier Connectivity Map and LINCS efforts to a higher-resolution cellular setting [5-7]. Predictive models built on these data are therefore attractive for prioritizing untested perturbation-cell combinations and for comparing representations of cellular state and chemical structure [8-22].

Benchmark validity is a central issue in this setting [23-29]. A random split can place the same drug, a near-neighbor scaffold, the same cell context, or a related mechanism of action in both the training and test partitions. Under this design, a model may score well by interpolating around familiar chemical or cellular neighborhoods rather than extrapolating to genuinely unseen perturbation settings. This problem is especially relevant for functional-genomics benchmarks because the biological interpretation of a prediction depends on whether the evaluation task resembles the intended use case.

Several model families now compete in single-cell perturbation prediction, including variational perturbation models, chemical-structure-aware models, neural optimal-transport methods, and single-cell foundation-model representations. These methods differ in input contract: some require raw single-cell profiles, some require matched controls, and some can operate on pseudobulk response matrices. A benchmark resource should therefore record not only performance but also the adaptation level and feasibility status of each model. Otherwise, hidden compatibility choices can become another source of irreproducibility.

Here we position the study as a leakage-audit protocol and resource rather than as a claim that all modern perturbation models fail under strict validation. The incremental contribution is fourfold: a joint cell-context-by-scaffold split for single-cell chemical perturbation prediction; a MoA-held-out audit showing that scaffold independence is not equivalent to mechanism independence; a rank-transfer analysis showing that random validation does not reliably predict strict-split model ordering; and a public-resource structure with split manifests, leakage audits, model-feasibility/status rows, and source data. We first audit OpenProblems NeurIPS 2023 perturbation prediction under random, cell held-out, scaffold held-out, and joint cell-plus-scaffold held-out splits. We then test whether the same audit patterns appear in Sci-Plex 3, a large single-cell chemical transcriptomics study aggregated here to drug-cell-line-dose pseudobulk responses.

## Methods

### Study design

This study audited benchmark leakage and model-ranking stability in single-cell chemical perturbation transcriptomics. The primary dataset was the OpenProblems NeurIPS 2023 perturbation response matrix. External validation used a Sci-Plex 3 24 h pseudobulk response matrix aggregated at the drug-cell-line-dose level. The analysis compared random validation with cell or cell-line held-out, scaffold held-out, and joint cell-plus-scaffold held-out validation. An OpenProblems MoA-held-out split was added as a secondary mechanism-level stress test.

### OpenProblems dataset

OpenProblems treated response records were merged across the local training and test AnnData files. Records were annotated with sample identifiers, drug names, SMILES strings, cell context, and response vectors. Drug structures were converted to Bemis-Murcko scaffolds with RDKit. MoA annotations were merged through normalized drug names where available. Response targets used the `clipped_sign_log10_pval` matrix, and DE-gene masks were taken from the corresponding `is_de` layer.

### Sci-Plex 3 external validation dataset

Sci-Plex 3 raw AnnData was processed into a 24 h pseudobulk response matrix. Records were aggregated at the drug-cell-line-dose level after filtering to treated groups with sufficient cells and to the top 2,000 genes used in the current analysis. The 2,000 genes were selected by total raw-count support in the processed Sci-Plex 3 matrix before model evaluation; this selection was split-agnostic and did not use model predictions or train-test assignments. Drug names were mapped to PubChem SMILES, and Bemis-Murcko scaffolds were computed from those structures. The resulting response matrix contained 2,200 drug-cell-line-dose records.

### Split design

For OpenProblems, four primary split types were used: random, cell held-out, scaffold held-out, and joint cell-plus-scaffold held-out. The joint split assigned records to the test set only when both the cell context and drug scaffold were held out, assigned records to the training set only when neither was held out, and excluded intermediate records. The MoA-held-out secondary split sampled held-out MoA labels and assigned all records with those MoA labels to the test set.

For Sci-Plex 3, the same structure was applied with cell line replacing cell context. Existing split manifests were treated as fixed benchmark inputs. Model adapters were not allowed to modify the assignment of training, test, or excluded records.

### Leakage audit

For each split and seed, we computed overlap between test records and training records for drug identity, scaffold identity, cell context or cell line, drug-cell pair, and MoA where available. OpenProblems MoA overlap was evaluated among annotated records; in the current MoA audit tables, all test records in the primary splits had MoA annotations. We also computed nearest training-drug Tanimoto similarity from Morgan fingerprints. In Sci-Plex 3, the audit further included same drug-dose, same drug-cell-line-dose, and same numeric dose overlap.

### Completed baseline and stress-test models

Completed predictors included global training-set mean, cell-context mean, drug mean, cell-only ridge, drug-fingerprint ridge, cell-plus-drug-fingerprint ridge, nearest-drug retrieval, and SVD-ridge stress-test variants. OpenProblems main ridge results used 100 split seeds where available; nearest-drug and SVD-ridge stress tests used existing 10-seed outputs. Sci-Plex 3 external validation used 30 seeds. MoA-held-out OpenProblems baselines used 100 seeds.

### Control-conditioned neural predictor and PRnet-interface adapter in Sci-Plex 3

The main neural arm in Sci-Plex 3 was a control-conditioned response regressor. The internal model identifier remains `transigen_adapted_sciplex3` in result files for provenance, but the benchmark run was an auditable local architecture rather than a line-by-line reproduction of the published TranSiGen workflow. The model used matched cell-line 24 h control expression, Morgan drug fingerprints, log-transformed dose, and cell-line one-hot encodings to predict treated-minus-control response vectors across the same 2,000 genes used in the Sci-Plex 3 pseudobulk benchmark (Table 7). Training used only records assigned to the training split; early stopping validation was sampled from the training records only. Test records were used only for final evaluation.

PRnet was checked against the official repository and demo files. The current Sci-Plex 3 benchmark matrix was not a drop-in replacement for PRnet's original AnnData input contract, which expects profile-level rows with SMILES, dose, paired control indices, and split annotations. We therefore created `prnet_adapted_sciplex3`, a pseudobulk adapter that imports the official PRnet PGM module and trains on basal/control expression plus dose-scaled drug fingerprints. The repeated-seed benchmark used fixed manifest seeds 1-30, with model initialization seed matched to the split seed. Validation records were sampled only from the training set; test records were used only for final evaluation.

Both neural analyses are benchmark-compatible pseudobulk adaptations. They are not full reproductions of the original model training and preprocessing workflows.

### Model feasibility and status reporting

The resource records whether each model requires raw single-cell input, supports pseudobulk response matrices, uses drug structure, uses dose, and uses cell context. Planned adapters for chemCPA/CPA and frozen foundation-model response heads are retained as reproducible extension points. They do not generate performance values unless the required external implementation and input builder are available. A skipped or failed model remains represented as a status row rather than disappearing from the benchmark.

### Metrics and statistical summaries

Completed models were summarized with all-gene row-wise Pearson, DE-gene-only row-wise Pearson where available, all-gene RMSE, and DE-gene RMSE where available. Top-k response gene recovery and direction agreement were retained from existing OpenProblems analyses. Model ranking was computed within split and seed, using DE-gene Pearson when available and all-gene Pearson otherwise. Random-to-strict contrasts were computed for random minus scaffold held-out and random minus joint held-out performance. Bootstrap confidence intervals used seed-level resampling for repeated-seed split means, paired random-to-strict contrasts, and paired model comparisons where available. Confidence intervals were treated as descriptive uncertainty summaries rather than a family-wise multiple-testing procedure; no multiplicity-adjusted inferential claims were made.

### Pathway-level response recovery analysis

Pathway-level recovery was evaluated as a gene-set-level prediction metric rather than as a pathway discovery analysis. We used the local MSigDB Hallmark GMT file provided in the workspace and did not download gene sets during this revision. Gene symbols were uppercased and matched directly to each response matrix. A gene set was retained for a dataset only when at least 10 member genes were present in the evaluated response space. This retained 47 of 50 Hallmark sets in OpenProblems and 39 of 50 Hallmark sets in the Sci-Plex 3 top-2,000-gene matrix.

For each retained Hallmark set, the true pathway response score was defined as the mean response value across overlapping member genes. The predicted pathway score was computed in the same way from the model-predicted response vector. We then computed record-level Pearson correlation and RMSE across Hallmark pathway scores, per-gene-set absolute error, seed-level summaries, paired random-to-strict contrasts, and bootstrap 95% confidence intervals using seed as the resampling unit. This analysis included only models with full prediction vectors or reconstructable train-only prediction vectors. OpenProblems global mean and ridge cell-plus-drug-fingerprint baselines were reconstructed from local matrices and fixed split manifests. OpenProblems nearest-drug and SVD-ridge outputs were not projected because the available local files contained per-row metrics but not full prediction vectors. Sci-Plex 3 global mean, ridge drug-fingerprint, ridge cell-plus-dose-plus-drug-fingerprint, control-conditioned response regressor, and PRnet-interface adapter predictors were evaluated from reconstructed baseline vectors or completed prediction-vector manifests.

### Quality control

The resource includes automated checks for metric behavior on toy data, split integrity, no-test-tuning configuration, and output completeness. Additional neural-adapter and PRnet-specific tests check prediction-target shape, record alignment, split integrity, excluded-record handling, and the no-test-leakage rule. Scaffold held-out splits must have zero same-drug and same-scaffold overlap. Joint held-out splits must have zero same-drug, same-scaffold, and same-cell overlap. MoA-held-out splits must have zero same-MoA overlap. Failed models must remain represented by status rows rather than disappearing from the output table.

## Results

### Leakage-aware split construction for single-cell perturbation transcriptomics

The OpenProblems dataset contained 545 treated records, 140 drugs, and 136 valid Bemis-Murcko scaffolds. Across 100 repeated seeds, random splits contained a mean of 444.6 training records and 108.4 test records. Cell held-out splits contained 416.7 training records and 136.3 test records. Scaffold held-out splits contained 446.3 training records and 106.7 test records. Joint cell-plus-scaffold held-out splits contained 333.4 training records, 26.8 test records, and 192.8 excluded records. The Sci-Plex 3 external validation matrix contained 2,200 drug-cell-line-dose response records (Fig. 1; Table 1; Table 2).

### Random validation overestimates transcriptomic perturbation prediction

In OpenProblems, the ridge cell-plus-drug fingerprint baseline achieved mean all-gene row-wise Pearson values of 0.362, 0.234, and 0.118 under random, scaffold held-out, and joint cell-plus-scaffold held-out validation. The corresponding DE-gene-only Pearson values were 0.568, 0.392, and 0.209. The random-to-scaffold and random-to-joint DE-gene Pearson differences were 0.176 and 0.359. Size-matched and composition-matched random controls remained near the original random-split performance, whereas actual joint held-out performance remained substantially lower (Fig. 2; Table 3). Because the OpenProblems joint split contained only 26.8 test records on average, this result is interpreted together with repeated seeds and matched random controls rather than as a single stable correlation estimate.

### Chemical-neighbor shortcuts account for random-split advantage

Random OpenProblems splits retained same-drug training neighbors for 98.9% of test observations and same-scaffold training neighbors for 97.5%. Scaffold held-out and joint held-out splits reduced same-drug and same-scaffold overlap to zero. Nearest training-drug Tanimoto similarity showed the same pattern: random validation preserved local chemical neighborhoods, whereas scaffold-strict validation forced extrapolation away from close scaffold neighbors (Fig. 3).

### Mechanism-held-out validation separates scaffold independence from MoA-level extrapolation

Scaffold independence did not imply mechanism independence. Scaffold-strict and joint-strict OpenProblems splits still retained same-MoA training neighbors for approximately 42.7% of test records. This estimate was not driven by sparse annotation coverage: the MoA audit tables reported a mean annotated fraction of 1.0 for test records in the primary splits. The MoA-held-out secondary split removed same-MoA training overlap by construction. Across 100 seeds, this split retained a mean of 434.2 training records and 110.8 test records. In this task, cell-only ridge achieved all-gene Pearson of 0.301 and DE-gene Pearson of 0.468, whereas ridge cell-plus-drug fingerprint achieved 0.212 and 0.363. Drug-fingerprint-only ridge achieved 0.141 and 0.249. Thus, simple drug fingerprints did not provide a transferable advantage when the held-out task crossed annotated mechanism boundaries.

### Sci-Plex 3 external validation reproduces random-to-strict performance decay

Sci-Plex 3 external validation showed the same qualitative pattern. The ridge cell-plus-dose-plus-drug-fingerprint baseline reached all-gene Pearson values of 0.292, 0.189, and 0.060 under random, scaffold held-out, and joint cell-line-plus-scaffold held-out validation. The random-to-scaffold and random-to-joint differences were 0.103 and 0.232. Unlike OpenProblems, cell-line held-out was harder than scaffold held-out in Sci-Plex 3, consistent with stronger baseline-expression separation among cancer cell lines than among the OpenProblems immune cell contexts.

### Control-conditioned neural prediction remains vulnerable under joint extrapolation

The Sci-Plex 3 control-conditioned response regressor and PRnet-interface adapter were run as model-integration stress tests rather than as claims about the original model workflows (Fig. 4; Table 3; Table 4; Table 7). The control-conditioned response regressor retained measurable random-split signal, with mean all-gene Pearson of 0.233 under random validation and 0.212 under scaffold held-out validation, but performance dropped to 0.017 under joint cell-line-plus-scaffold held-out validation. The paired random-minus-scaffold difference was 0.021 (bootstrap 95% CI, 0.013 to 0.029), and the paired random-minus-joint difference was 0.216 (0.208 to 0.224). In paired comparison with the ridge cell-plus-dose-plus-drug-fingerprint baseline, the neural regressor was lower under random validation by 0.059 (95% CI, 0.056 to 0.062), higher under scaffold held-out validation by 0.023 (0.019 to 0.027), and lower under joint held-out validation by 0.043 (0.032 to 0.054).

The PRnet-interface adapter reached 0.017 under random validation (bootstrap 95% CI, 0.013 to 0.022), 0.015 under scaffold held-out validation (0.012 to 0.018), and 0.001 under joint cell-line-plus-scaffold held-out validation (-0.013 to 0.014) across 30 seeds. Because the random-split score was already close to zero, this arm is treated as an executable adapter and feasibility result rather than as primary evidence for random-to-strict performance decay. The PRnet-style pseudobulk adaptation was executable but did not improve strict extrapolation under this benchmark-compatible adaptation. This result should not be interpreted as evidence that the original PRnet workflow fails.

A 10-seed input ablation of the control-conditioned response regressor suggested that fingerprint features were most useful when chemical neighborhoods were still partly available. The drug-only arm decreased from 0.232 under random validation to 0.166 under scaffold held-out validation, a paired drop of 0.067 (95% CI, 0.051 to 0.084). In contrast, the full basal-plus-drug-plus-dose-plus-cell arm had a smaller random-minus-scaffold drop of 0.019 (0.005 to 0.033), and basal-only performance was similar between random and scaffold held-out validation. The ablation was retained as a diagnostic analysis because it used 10 seeds rather than the 30-seed main neural run.

### Pathway-level recovery deteriorates under strict extrapolation

Hallmark gene-set-level response recovery showed that the random-to-strict decay was not limited to gene-wise Pearson summaries (Fig. 5). In OpenProblems, ridge cell-plus-drug-fingerprint pathway Pearson decreased from 0.552 under random validation to 0.317 under scaffold held-out validation and 0.209 under joint cell-plus-scaffold held-out validation across 10 available vector-level split seeds. The paired random-minus-joint pathway difference was 0.343, with bootstrap 95% CI from 0.293 to 0.389. The global mean baseline showed smaller pathway-level differences, consistent with its lower dependence on drug-neighbor features.

In Sci-Plex 3, ridge cell-plus-dose-plus-drug-fingerprint pathway Pearson decreased from 0.313 under random validation to 0.187 under scaffold held-out validation and 0.070 under joint cell-line-plus-scaffold held-out validation across 30 seeds. The paired random-minus-joint pathway difference was 0.243 (95% CI, 0.223 to 0.262). The control-conditioned response regressor also showed pathway-level decay, from 0.245 under random validation to 0.228 under scaffold held-out validation and 0.027 under joint held-out validation; the random-minus-joint pathway difference was 0.218 (0.202 to 0.235). The PRnet-interface adapter had low pathway-level recovery across all splits, with 0.046 under random validation, 0.032 under scaffold held-out validation, and 0.047 under joint held-out validation. Therefore, the pathway-level decay conclusion is anchored in the ridge and control-conditioned neural predictors, which had clearer random-split signal. These results indicate that leakage-aware validation affects recovery of aggregated perturbation programs as well as individual-gene response profiles.

### Random validation does not preserve model ranking

Model ranking transfer was weak across strict validation settings (Fig. 6; Table 5). In OpenProblems, aggregate random-to-joint rank Spearman correlation across completed models was 0.060, and the random-best model dropped to rank 6 under joint held-out validation. After adding the 30-seed PRnet run to a five-model Sci-Plex 3 seed-level panel, the mean random-to-joint Spearman correlation was 0.143 across seeds. This rank instability is one of the most practically important outputs of the audit because it shows that random validation can mislead model selection, not only absolute performance estimates.

## Discussion

### Principal findings

This study provides a leakage-aware benchmark resource for single-cell chemical perturbation transcriptomics. The two main findings are that random validation can overestimate response prediction and that random-split model rankings may not transfer to strict cellular-chemical extrapolation tasks. Three secondary findings support this interpretation. First, the random-split advantage aligned with same-drug and same-scaffold training overlap. Second, MoA-held-out validation showed that scaffold independence is not mechanism independence. Third, Hallmark gene-set-level recovery also deteriorated under strict extrapolation for predictors with clear random-split signal. The pseudobulk neural adapters show that modern model interfaces can be connected to the resource, but they should be interpreted as benchmark-compatible adaptations rather than full reproductions of third-party workflows.

### Random validation as local interpolation rather than extrapolation

Random validation can preserve several layers of shortcut information. Same-drug overlap gives the model direct access to perturbation identity. Same-scaffold overlap gives access to close chemical-neighbor structure. Same-cell-context overlap preserves cellular baseline state. These shortcuts are not equivalent, but each can turn a nominal test record into local interpolation. This matters for perturbation transcriptomics because a high random-split score may reflect familiarity with local chemical or cellular neighborhoods rather than recovery of a transferable transcriptional response program.

The pathway-level analysis strengthens this interpretation and connects the benchmark to established gene-set-level interpretation frameworks [30-32]. Hallmark scores aggregate genes into broad perturbation programs, but they did not remove the random-to-strict gap for models that depended on chemical or cellular-neighbor features. Thus, reporting only pathway-level recovery would not by itself protect a benchmark from local-interpolation artifacts. Instead, pathway recovery should be interpreted alongside the same split-integrity and leakage-audit fields used for gene-level metrics.

### Scaffold independence is not mechanism independence

Scaffold-held-out validation is necessary for evaluating chemical extrapolation, but it is not sufficient for mechanism-level extrapolation [26-29]. Drugs with different Bemis-Murcko scaffolds can still share annotated MoA or converge on related stress and regulatory programs. The OpenProblems MoA-held-out split therefore provides a secondary functional-genomics stress test. In this setting, simple drug fingerprints did not improve over a cell-only ridge baseline, suggesting that the tested linear fingerprint representation did not transfer across annotated mechanism boundaries in this dataset.

### Why model feasibility is part of benchmark validity

Modern perturbation models do not share a single input contract. Some require raw single-cell counts, some require matched controls, some require dose covariates, and some require pretrained embeddings. Treating these requirements as informal preprocessing choices can make a benchmark difficult to reproduce. We therefore record model adaptation level and feasibility status as first-class benchmark outputs. The control-conditioned response regressor and PRnet-interface adapter show that the resource can admit neural perturbation interfaces, but they are not sufficient to support broad claims about modern model families. Their main role here is to stress-test the audit framework and to document what was executable under the current pseudobulk benchmark contract.

### Implications for single-cell foundation model evaluation

The current draft does not claim that foundation models fail under scaffold-strict validation, because scGPT, scFoundation, Geneformer, and related embedding arms were not completed in this benchmark. The relevant implication is narrower: any foundation-model response head should be evaluated under the same leakage audit as simpler baselines. Frozen embeddings may improve cellular-state representation, but they do not by themselves remove same-drug, same-scaffold, same-MoA, or nearest-neighbor chemical shortcuts.

### Practical reporting recommendations for perturbation transcriptomics benchmarks

Future benchmark reports should include the split type, same-drug overlap, same-scaffold overlap, same-cell-context overlap, MoA-neighbor overlap, nearest training-drug Tanimoto, random-to-strict performance drop, random-to-strict model-rank transfer, model adaptation level, and failed or not-run status rows (Table 6), following broader expectations for reproducible dataset, model, and benchmark reporting [33-40]. These fields help separate biological generalization from convenient interpolation and make the benchmark more useful as a reusable functional-genomics resource.

### Limitations

This study has limitations. The control-conditioned response regressor and PRnet-interface adapter are pseudobulk adaptations in Sci-Plex 3, not full reproductions of the original TranSiGen or PRnet training workflows. OpenProblems was not used for either adaptation because the current OpenProblems response matrix lacks matched basal/control expression profiles. chemCPA/CPA and foundation-model arms remain feasibility/status rows, so the manuscript should be read as a leakage-audit resource rather than a definitive comparison of modern perturbation models. The PRnet-interface adapter was executable but had near-zero random-split performance, and it is therefore not used as primary evidence for random-to-strict decay. Sci-Plex 3 DE-gene-only metrics were not reported because the current pseudobulk artifact did not contain per-record DE masks aligned to the 2,000-gene response vectors. Training-loss and validation-loss curves were also not persisted for the completed neural runs, so convergence evidence is limited to completion status, validation-based early stopping in code, and seed-level metric stability. Pathway-level analysis was limited to local MSigDB Hallmark gene sets and direct gene-symbol overlap. It should therefore be interpreted as a pathway-level recovery metric rather than pathway discovery or mechanistic validation. OpenProblems nearest-drug and SVD-ridge predictions were not included in the pathway analysis because full prediction vectors were not available locally. The OpenProblems joint held-out split is small, although matched random controls reduce the likelihood that the main conclusion is caused by test-set size alone. The MoA-held-out analysis depends on current MoA annotation granularity. Sci-Plex 3 results use a 24 h top-2,000-gene pseudobulk matrix and a split-agnostic top-2,000-gene filter based on raw-count support, which may underrepresent low-expression but strongly responsive genes and should be stress-tested with alternative gene-selection rules in future releases.

## Conclusions

Leakage-aware validation changes the interpretation of single-cell chemical perturbation transcriptomics benchmarks. Random validation can overestimate response-prediction performance and can produce model rankings that do not transfer to scaffold-strict or joint held-out tasks. The benchmark resource provides fixed leakage-aware splits, audit tables, long-format model outputs, adapter interfaces, feasibility reporting, pathway-level source data, and QC checks. Future claims about perturbation model performance should be interpreted together with chemical-neighbor leakage, mechanism-neighbor leakage, model adaptation level, pathway-level recovery, and rank-transfer stability.

## Supplementary Information

Supplementary materials should include complete split manifests, seed-level model summaries, top-k response gene recovery, leakage overlap audits, MoA-held-out results, dose leakage audits, matched random controls, neural input-ablation diagnostics, rank-transfer tables, pathway-level recovery tables, pathway gene-set overlap checks, model feasibility reports, model registry files, benchmark report schemas, figure source data, and environment records.

## Acknowledgements

The authors have no specific acknowledgements to declare.

## Authors' contributions

Da Lin contributed to conceptualization, methodology, data curation, formal analysis, software implementation, visualization, validation, and writing of the original draft. Yu Zhang contributed to conceptualization, supervision, project administration, interpretation of results, and critical revision of the manuscript. Both authors read and approved the final manuscript.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Availability of data and materials

The local benchmark resource package contains split manifests, leakage-audit tables, model-result summaries, report schemas, model feasibility reports, data manifests, source data, and configuration files. Code and benchmark resources are available at https://github.com/seefreewind/single-cell-perturbation-leakage-benchmark. A Zenodo archive will be created before submission and the DOI will be added after deposition. MSigDB Hallmark gene sets are not redistributed in this repository and should be obtained by users from MSigDB according to its terms of use.

## Ethics approval and consent to participate

Not applicable. This study used publicly available or locally mirrored benchmark datasets and did not involve new human participant recruitment.

## Consent for publication

Not applicable.

## Competing interests

The authors declare that they have no competing interests. This statement should be confirmed before submission.

## References

1. Dixit A, Parnas O, Li B, et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. Cell. 2016;167:1853-1866.e17. doi:10.1016/j.cell.2016.11.038.
2. Adamson B, Norman TM, Jost M, et al. A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response. Cell. 2016;167:1867-1882.e21. doi:10.1016/j.cell.2016.11.048.
3. Datlinger P, Rendeiro AF, Schmidl C, et al. Pooled CRISPR screening with single-cell transcriptome readout. Nature Methods. 2017;14:297-301. doi:10.1038/nmeth.4177.
4. Srivatsan SR, McFaline-Figueroa JL, Ramani V, et al. Massively multiplex chemical transcriptomics at single-cell resolution. Science. 2020;367:45-51. doi:10.1126/science.aax6234.
5. Lamb J, Crawford ED, Peck D, et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science. 2006;313:1929-1935. doi:10.1126/science.1132939.
6. Subramanian A, Narayan R, Corsello SM, et al. A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles. Cell. 2017;171:1437-1452.e17. doi:10.1016/j.cell.2017.10.049.
7. Koleti A, Terryn R, Stathias V, et al. Data Portal for the Library of Integrated Network-based Cellular Signatures (LINCS) program: integrated access to diverse large-scale cellular perturbation response data. Nucleic Acids Research. 2018;46:D558-D566. doi:10.1093/nar/gkx1063.
8. OpenProblems Single-Cell Analysis. Perturbation Prediction task. https://openproblems.bio/
9. OpenProblems bio task repository. task_perturbation_prediction. https://github.com/openproblems-bio/task_perturbation_prediction
10. A benchmark for prediction of transcriptomic responses to chemical perturbations. NeurIPS Datasets and Benchmarks Track. 2024. https://openreview.net/forum?id=WTI4RJYSVm
11. Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. Nature Methods. 2019;16:715-721. doi:10.1038/s41592-019-0494-8.
12. Lotfollahi M, Klimovskaia Susmelj A, De Donno C, et al. Predicting cellular responses to complex perturbations in high-throughput screens. Molecular Systems Biology. 2023;19:e11517. doi:10.15252/msb.202211517.
13. Hetzel L, Boehm S, Kilbertus N, Guennemann S, Lotfollahi M, Theis FJ. Predicting cellular responses to novel drug perturbations at a single-cell resolution. Advances in Neural Information Processing Systems. 2022. https://openreview.net/forum?id=vRrFVHxFiXJ
14. Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. Nature Biotechnology. 2024;42:927-935. doi:10.1038/s41587-023-01905-6.
15. Bunne C, Stark SG, Gut G, et al. Learning single-cell perturbation responses using neural optimal transport. Nature Methods. 2023;20:1759-1768. doi:10.1038/s41592-023-01969-x.
16. Tong X, Qu N, Kong X, et al. Deep representation learning of chemical-induced transcriptional profile for phenotype-based drug discovery. Nature Communications. 2024. doi:10.1038/s41467-024-49620-3.
17. Qi X, Zhao L, Tian C, et al. Predicting transcriptional responses to novel chemical perturbations using deep generative model for drug discovery. Nature Communications. 2024. doi:10.1038/s41467-024-53457-1.
18. Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0.
19. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. Nature. 2023;618:616-624. doi:10.1038/s41586-023-06139-9.
20. Hao M, Gong J, Zeng X, et al. Large scale foundation model on single-cell transcriptomics. bioRxiv. 2023. doi:10.1101/2023.05.29.542705.
21. Lotfollahi M, Naghipourfar M, Luecken MD, et al. Mapping single-cell data to reference atlases by transfer learning. Nature Biotechnology. 2022;40:121-130. doi:10.1038/s41587-021-01001-7.
22. Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. Nature Methods. 2022;19:41-50. doi:10.1038/s41592-021-01336-8.
23. Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. Patterns. 2023;4:100804. doi:10.1016/j.patter.2023.100804.
24. Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics. 2006;7:91. doi:10.1186/1471-2105-7-91.
25. Yarkoni T, Westfall J. Choosing prediction over explanation in psychology: lessons from machine learning. Perspectives on Psychological Science. 2017;12:1100-1122. doi:10.1177/1745691617693393.
26. Wu Z, Ramsundar B, Feinberg EN, et al. MoleculeNet: a benchmark for molecular machine learning. Chemical Science. 2018;9:513-530. doi:10.1039/C7SC02664A.
27. Bemis GW, Murcko MA. The properties of known drugs. 1. Molecular frameworks. Journal of Medicinal Chemistry. 1996;39:2887-2893. doi:10.1021/jm9602928.
28. Rogers D, Hahn M. Extended-connectivity fingerprints. Journal of Chemical Information and Modeling. 2010;50:742-754. doi:10.1021/ci100050t.
29. Landrum G. RDKit: Open-source cheminformatics software. https://www.rdkit.org/
30. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Systems. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004.
31. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America. 2005;102:15545-15550. doi:10.1073/pnas.0506580102.
32. Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the unification of biology. Nature Genetics. 2000;25:25-29. doi:10.1038/75556.
33. Wilkinson MD, Dumontier M, Aalbersberg IJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data. 2016;3:160018. doi:10.1038/sdata.2016.18.
34. Gebru T, Morgenstern J, Vecchione B, et al. Datasheets for datasets. Communications of the ACM. 2021;64:86-92. doi:10.1145/3458723.
35. Mitchell M, Wu S, Zaldivar A, et al. Model cards for model reporting. Proceedings of the Conference on Fairness, Accountability, and Transparency. 2019:220-229. doi:10.1145/3287560.3287596.
36. Pineau J, Vincent-Lamarre P, Sinha K, et al. Improving reproducibility in machine learning research. Journal of Machine Learning Research. 2021;22:1-20. https://www.jmlr.org/papers/v22/20-303.html
37. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods. 2020;17:261-272. doi:10.1038/s41592-019-0686-2.
38. Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020;585:357-362. doi:10.1038/s41586-020-2649-2.
39. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. Journal of Machine Learning Research. 2011;12:2825-2830. https://jmlr.org/papers/v12/pedregosa11a.html
40. Hunter JD. Matplotlib: a 2D graphics environment. Computing in Science & Engineering. 2007;9:90-95. doi:10.1109/MCSE.2007.55.
