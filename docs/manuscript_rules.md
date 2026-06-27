# Manuscript Rules

This file preserves project-level writing and audit rules for the manuscript.

## Article framing

The manuscript should be framed as a leakage-resistant benchmark and failure analysis of single-cell foundation models for drug perturbation prediction. It should not be framed primarily as a new complex model paper unless a genuinely new method is later developed and validated.

## Default structure

Use a research-article structure:

1. Title
2. Abstract
3. Keywords
4. Background or Introduction
5. Methods
6. Results
7. Discussion
8. Conclusions
9. Supplementary Information
10. Acknowledgements
11. Authors' contributions
12. Funding
13. Availability of data and materials
14. Ethics approval and consent to participate
15. Consent for publication
16. Competing interests
17. References

## Claim strength

- Use restrained language: "was associated with", "showed evidence of", "may contribute to", "suggests", "requires further validation".
- Avoid unsupported claims such as "proved", "confirmed the mechanism", "clinically validated", "highly accurate biomarker", or "therapeutic target" unless direct validation supports them.
- Discovery-cohort, public-data, machine-learning, and benchmark findings must not be described as clinically validated.

## Citation placement

- Use in-text references for background, rationale, prior-work comparison, interpretation, and limitations.
- Avoid literature citations in Methods and Results unless a method, database, or software citation is required.

## Methods and Results

- Methods must report data source, accession or URL, sample size, inclusion criteria, preprocessing, software versions, parameters, statistical thresholds, and multiple-testing correction.
- Results should follow the analysis workflow: dataset characteristics, primary benchmark, leakage audit, ablations, external validation.
- Results should be descriptive and statistical. Avoid speculative phrases such as "this may be due to" in Results.

## Figures and tables

- Main figures must be self-contained.
- Figure legends should define axes, colors, thresholds, abbreviations, sample sizes, and statistical tests.
- Tables should include effect estimates, standard errors or confidence intervals where applicable, P values, adjusted P values or FDR, and threshold notes.

