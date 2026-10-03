#!/usr/bin/env python3
"""Generate BMC manuscript tables from the BMC master result source."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT / "results/bmc_master"
OUT = ROOT / "BMC_Bioinformatics_submission/assets/tables"


def write(name: str, lines: list[str]) -> None:
    (OUT / name).write_text("\n".join(lines) + "\n", encoding="utf-8")


def dataset_table() -> None:
    d = pd.read_csv(MASTER / "master_dataset_characteristics.csv")
    rows = []
    for _, x in d.iterrows():
        scaff = "--" if pd.isna(x.unique_scaffolds_nonmissing) else str(int(x.unique_scaffolds_nonmissing))
        rows.append(f"{x.dataset} & {x.stage} & {int(x.records)} & {int(x.unique_drugs)} & {int(x.unique_cell_contexts)} & {scaff} & {int(x.unique_doses)} \\\\")
    write("table_dataset_characteristics.tex", rows)


def split_table() -> None:
    s = pd.read_csv(MASTER / "master_split_summary.csv")
    s = s[s.count_stage.eq("canonical_manuscript_reporting")].copy()
    rows = []
    for (dataset, split), g in s.groupby(["dataset", "split"]):
        values = g.set_index("quantity")
        def cell(key: str) -> str:
            x = values.loc[key]
            return f"{x['mean']:.1f} $\\pm$ {x['sd']:.1f}"
        rows.append(f"{dataset} & {values.iloc[0]['split_label']} & {int(values.loc['manifest_test', 'n_seeds'])} & {cell('manifest_train')} & {cell('manifest_test')} & {cell('manifest_excluded')} & {cell('evaluable_test')} \\\\")
    write("table_split_accounting.tex", rows)


def probe_table() -> None:
    p = pd.read_csv(MASTER / "master_performance_summary.csv")
    configurations = [
        ("OpenProblems", ["global_train_mean", "cell_context_mean", "drug_mean", "ridge_cell", "ridge_drug_fp", "ridge_cell_drug_fp"]),
        ("Sci-Plex 3", ["global_train_mean", "cell_context_mean", "drug_mean", "ridge_drug_fp", "ridge_cell_dose_drug_fp", "control_conditioned_nonlinear_benchmark"]),
    ]
    labels = {
        "global_train_mean": "Global mean", "cell_context_mean": "Cell context/line mean", "drug_mean": "Drug mean",
        "ridge_cell": "Cell-context ridge", "ridge_drug_fp": "Drug-fingerprint ridge",
        "ridge_cell_drug_fp": "Cell-context + fingerprint ridge", "ridge_cell_dose_drug_fp": "Cell-line + dose + fingerprint ridge",
        "control_conditioned_nonlinear_benchmark": "Control-conditioned nonlinear benchmark",
    }
    rows = []
    for dataset, models in configurations:
        d = p[(p.dataset.eq(dataset)) & p.metric.eq("rowwise_pearson")]
        for model in models:
            m = d[d.model.eq(model)].set_index("split_label")
            if not {"random", "scaffold held-out", "joint held-out"}.issubset(m.index):
                continue
            value = lambda split: f"{m.loc[split, 'mean']:.3f} $\\pm$ {m.loc[split, 'sd']:.3f}"
            rows.append(f"{dataset} & {labels[model]} & {value('random')} & {value('scaffold held-out')} & {value('joint held-out')} \\\\")
    write("table_multiprobe_performance.tex", rows)


def control_table() -> None:
    c = pd.read_csv(MASTER / "master_contrast_summary.csv")
    c = c[c.comparison.eq("train_test_size_matched_random minus joint held-out")].copy()
    labels = {"ridge_cell_drug_fp": "Cell-context + fingerprint ridge", "control_conditioned_nonlinear_benchmark": "Control-conditioned nonlinear benchmark"}
    rows = []
    for _, x in c.iterrows():
        rows.append(f"{x.dataset} & {labels[x.model]} & {int(x.n_paired_seeds)} & {x.mean_difference:.3f} & {x.bootstrap_95ci_low:.3f} to {x.bootstrap_95ci_high:.3f} \\\\")
    write("table_train_size_matched_contrasts.tex", rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dataset_table(); split_table(); probe_table(); control_table()


if __name__ == "__main__":
    main()
