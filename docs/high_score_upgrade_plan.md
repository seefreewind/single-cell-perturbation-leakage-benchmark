# High-score Upgrade Plan

## Strategic Repositioning

Current manuscript status:

- One-dataset leakage audit on OpenProblems NeurIPS 2023.
- Strong baseline story: random validation is optimistic; scaffold-strict and joint held-out settings reduce performance.
- Main weakness: evidence is still single-dataset and simple-baseline heavy.

Target high-score framing:

> Chemical-neighbor leakage inflates apparent generalization in single-cell drug perturbation prediction, and leakage-aware scaffold-strict validation should become a standard audit component for single-cell drug perturbation benchmarks.

## Version 1: Statistical and Audit Hardening

Goal: make the current OpenProblems manuscript substantially more defensible.

Priority tasks:

1. Increase repeated splits from 10 to 100 seeds.
2. Add leakage overlap diagnostics:
   - same drug in train
   - same scaffold in train
   - same cell context in train
   - same drug-cell pair in train
   - nearest training-drug Tanimoto
3. Add model-ranking instability:
   - model ranks per split
   - random-vs-strict rank shifts
   - Spearman rank correlation across split regimes
4. Add hyperparameter sensitivity:
   - ridge alpha: 0.1, 1, 10, 100, 1000
   - SVD components: 10, 25, 50, 100
   - fingerprint length: 1024 vs 2048
5. Update figures and tables:
   - Move leakage audit earlier.
   - Add overlap leakage table.
   - Add ranking instability supplementary figure/table.
6. Update manuscript:
   - Sharpen OP3 distinction.
   - Clarify that foundation models are future audit targets unless real FM baselines are added.
   - Replace seed-level bootstrap language with 100-seed estimates.

## Version 2: Multi-dataset and Complex-model Expansion

Goal: move toward a stronger benchmark paper.

Priority tasks:

1. OP3 feasibility precheck:
   - availability of SMILES
   - drug-cell response matrix
   - cell context
   - scaffold counts
   - joint split size
2. scPerturb/Sci-Plex drug subset precheck:
   - at least 50 drugs
   - at least 30 valid scaffolds
   - at least 2 cell contexts
   - computable response matrix
3. Add chemCPA or CPA baseline on OpenProblems.
4. Add one practical foundation-model-adjacent baseline:
   - scGPT/Geneformer embedding + ridge, or
   - scVI/scArches latent + drug fingerprint ridge/MLP.

## Version 3: Full High-impact Benchmark

Goal: position for Genome Biology, Patterns, Cell Genomics, or similar.

Requirements:

1. At least two datasets, ideally three.
2. Multiple model families:
   - mean/linear/nearest/SVD
   - task-specific perturbation model
   - foundation model embedding baseline
3. Leakage-aware benchmark package:
   - reusable split generation
   - leakage audit checklist
   - metrics
   - source-data and figure scripts
4. Archived code and data package with DOI.

## Immediate Execution Order

1. Run 100 OpenProblems seeds using existing lightweight baselines.
2. Add leakage-overlap diagnostics.
3. Add ranking instability.
4. Add hyperparameter sensitivity.
5. In parallel, inspect OP3 data availability.

