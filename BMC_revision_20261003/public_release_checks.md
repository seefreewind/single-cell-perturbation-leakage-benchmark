# Public-release checks

The audit CLI example test passed (1 test, 3.51 seconds) on 2026-10-03 using the bundled Python 3.12 interpreter and the revision-local dependency directory. An initial attempt selected system Python 3.9 in a subprocess and failed on an incompatible RDKit extension. The test now uses `sys.executable`, keeping the subprocess in the same interpreter. This test exercises the audit CLI, not the full model-training workflow.

No model fitting or scientific summary regeneration was performed during this publication pass. Frozen CSVs were copied byte-for-byte. The existing nonlinear adapter change records training/validation loss and returns history; it is the previously used local implementation, not a new model experiment.

GitHub's selected connector reported read-only access. The existing authenticated `gh` account had repository write permission and was used for publication. No token was exported into the release.

Zenodo browser navigation and HTTPS API access failed certificate validation. No certificate warning was bypassed and no Zenodo publication is claimed. The earlier DOI remains the citation for the previous version until a new version is published.
