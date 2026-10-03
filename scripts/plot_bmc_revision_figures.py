#!/usr/bin/env python3
"""Create BMC revision figures only from results/bmc_master tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT / "results/bmc_master"
OUT = ROOT / "BMC_Bioinformatics_submission/assets"
FIG = OUT / "figures"
SRC = OUT / "source_data"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "legend.fontsize": 6,
    "figure.dpi": 150,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
})
COLORS = {"random": "#277DA1", "scaffold held-out": "#E8A13A", "joint held-out": "#C84B57"}
NEUTRAL = "#4D4D4D"
GRID = "#E8E8E8"


def panel_label(ax: plt.Axes, letter: str, x: float = -0.08, y: float = 1.06) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, fontweight="bold")


def label(model: str, dataset: str) -> str:
    names = {
        "global_train_mean": "Global mean",
        "cell_context_mean": "Cell-line mean" if dataset == "Sci-Plex 3" else "Cell-context mean",
        "drug_mean": "Drug mean",
        "ridge_drug_fp": "Drug-fingerprint ridge",
        "ridge_cell": "Cell-context ridge",
        "ridge_cell_drug_fp": "Cell-context + fingerprint ridge",
        "ridge_cell_dose_drug_fp": "Cell-line + dose + fingerprint ridge",
        "control_conditioned_nonlinear_benchmark": "Control-conditioned nonlinear benchmark",
    }
    return names.get(model, model.replace("_", " "))


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIG / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def fig1() -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.15))
    ax.axis("off")
    ax.text(.5, .98, "Leakage-aware benchmark evaluation defines what each score can claim", ha="center", va="top", weight="bold", fontsize=10)
    ax.text(.5, .91, "Generalization is interpreted as a model--data--split--claim system", ha="center", va="top", fontsize=7.2, color=NEUTRAL)

    datasets = [
        ("OpenProblems", "545 evaluable treated records\n138 stored drugs, 4 cell contexts\n100 fixed seeds"),
        ("Sci-Plex 3", "2,200 pseudobulk records\n188 drugs, 3 cell lines, 4 doses\n30 fixed seeds"),
    ]
    for i, (name, text) in enumerate(datasets):
        x0 = .07 + i * .45
        ax.add_patch(plt.Rectangle((x0, .745), .39, .12, facecolor="#F2F5F7", edgecolor="#8A99A8", linewidth=.8))
        ax.text(x0 + .02, .825, name, weight="bold", fontsize=7.6, va="center")
        ax.text(x0 + .02, .775, text, fontsize=6.0, va="center")

    layers = [
        ("A. Fixed tasks", "random local interpolation\ncell/context held-out\nscaffold held-out\njoint cell + scaffold held-out", "#E8F1FA"),
        ("B. Neighborhood audit", "same drug / scaffold / cell / dose / MoA\nnearest training-drug Tanimoto\nmissingness + integrity checks", "#EEF6E9"),
        ("C. Transparent probes", "mean and ridge baselines\ncontrol-conditioned nonlinear benchmark\ntrain-only validation and status rows", "#FDF2DF"),
        ("D. Portable evidence", "seed-level metrics\npaired matched-random contrasts\nsource data + reporting checklist", "#F7EAF4"),
    ]
    y0 = .405
    xs = [.04, .275, .51, .745]
    for (title, text, color), x0 in zip(layers, xs):
        ax.add_patch(plt.Rectangle((x0, y0), .205, .275, facecolor=color, edgecolor="#4A5568", linewidth=.8))
        ax.text(x0 + .012, y0 + .238, title, weight="bold", fontsize=6.8, va="center")
        ax.text(x0 + .012, y0 + .126, text, fontsize=5.15, va="center")
        if x0 < xs[-1]:
            ax.annotate("", xy=(x0 + .228, y0 + .137), xytext=(x0 + .205, y0 + .137), arrowprops={"arrowstyle": "->", "lw": .8, "color": "#4A5568"})

    tasks = [
        ("Random", "familiar drug/cell neighborhoods", COLORS["random"]),
        ("Scaffold held-out", "new chemical framework", COLORS["scaffold held-out"]),
        ("Joint held-out", "new cell + new scaffold intersection", COLORS["joint held-out"]),
    ]
    for i, (name, text, color) in enumerate(tasks):
        x0 = .12 + i * .28
        ax.add_patch(plt.Rectangle((x0, .14), .22, .12, facecolor="white", edgecolor=color, linewidth=1.1))
        ax.text(x0 + .11, .22, name, ha="center", va="center", weight="bold", color=color, fontsize=6.8)
        ax.text(x0 + .11, .175, text, ha="center", va="center", fontsize=5.5)
    ax.text(.5, .055, "Output: interpretable scores, not a single undifferentiated claim of model generalization", ha="center", fontsize=6.6, color=NEUTRAL)
    pd.concat([
        pd.DataFrame({"component_type": "dataset", "component": [x[0] for x in datasets], "description": [x[1].replace("\n", "; ") for x in datasets]}),
        pd.DataFrame({"component_type": "workflow_layer", "component": [x[0] for x in layers], "description": [x[1].replace("\n", "; ") for x in layers]}),
        pd.DataFrame({"component_type": "evaluation_task", "component": [x[0] for x in tasks], "description": [x[1] for x in tasks]}),
    ], ignore_index=True).to_csv(SRC / "Figure_1_source_data.csv", index=False)
    save(fig, "Figure_1_framework")


def fig2() -> None:
    leak = pd.read_csv(MASTER / "master_leakage_summary.csv")
    order = ["random", "scaffold held-out", "joint held-out"]
    d = leak[(leak.dataset.eq("OpenProblems")) & (leak.split_label.isin(order))].copy()
    d["split_label"] = pd.Categorical(d.split_label, order, ordered=True)
    d = d.sort_values("split_label")
    fields = [("same_drug_overlap", "same drug"), ("same_scaffold_overlap", "same scaffold"), ("same_cell_overlap", "same cell"), ("same_moa_overlap_among_annotated", "same MoA\n(annotated)")]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.1, 3.05), gridspec_kw={"width_ratios": [1.9, 1]})
    x = np.arange(len(fields)); width = .23
    for j, split in enumerate(order):
        values = d.loc[d.split_label.eq(split), [name for name, _ in fields]].iloc[0].to_numpy(float) * 100
        bars = ax.bar(x + (j - 1) * width, values, width, label=split, color=COLORS[split], edgecolor="white", linewidth=.4)
        for bar, val, field in zip(bars, values, [name for name, _ in fields]):
            if field == "same_cell_overlap":
                continue
            if val >= 94 or (split != "random" and val == 0):
                ax.text(bar.get_x() + bar.get_width() / 2, val + 2.1, f"{val:.0f}", ha="center", va="bottom", fontsize=5.3, color=NEUTRAL)
    ax.set_ylabel("Test records with a training neighbor (%)")
    ax.set_xticks(x, [x[1] for x in fields]); ax.set_ylim(0, 108)
    ax.yaxis.grid(True, color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.52, -0.20), ncol=3)
    ax.set_title("Neighborhood overlap", loc="left", weight="bold")
    panel_label(ax, "a", -0.07, 1.04)
    tan = d.set_index("split_label").loc[order, "mean_nearest_train_tanimoto"]
    bars = ax2.bar(range(3), tan, color=[COLORS[x] for x in order], edgecolor="white", linewidth=.4)
    for bar, val in zip(bars, tan):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + .035, f"{val:.3f}", ha="center", va="bottom", fontsize=5.8, color=NEUTRAL)
    ax2.set_ylim(0, 1.05); ax2.set_ylabel("Mean nearest training-drug\nTanimoto")
    ax2.set_xticks(range(3), ["random", "scaffold", "joint"], rotation=22, ha="right")
    ax2.yaxis.grid(True, color=GRID, linewidth=.6)
    ax2.set_axisbelow(True)
    ax2.set_title("Chemical proximity", loc="left", weight="bold")
    panel_label(ax2, "b", -0.16, 1.04)
    ax2.text(.02, -.32, "MoA overlap uses only annotated records; missingness is stored separately.", transform=ax2.transAxes, fontsize=6.2)
    d.to_csv(SRC / "Figure_2_source_data.csv", index=False)
    save(fig, "Figure_2_neighborhood_audit")


def fig3() -> None:
    perf = pd.read_csv(MASTER / "master_performance_seed_level.csv")
    contrast = pd.read_csv(MASTER / "master_contrast_summary.csv")
    settings = [
        ("OpenProblems", "ridge_cell_drug_fp", ["split_cell_scaffold_heldout", "test_size_matched_random", "train_test_size_matched_random", "composition_train_test_size_matched_random"], ["joint", "test-size\nmatched", "train + test\nmatched", "composition +\ntrain/test matched"]),
        ("Sci-Plex 3", "control_conditioned_nonlinear_benchmark", ["split_cell_scaffold_heldout", "test_size_matched_random", "train_test_size_matched_random"], ["joint", "test-size\nmatched", "train + test\nmatched"]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.25), sharey=False)
    out = []
    for ax, (dataset, model, order, xticks) in zip(axes, settings):
        d = perf[(perf.dataset.eq(dataset)) & (perf.model.eq(model)) & (perf.metric.eq("rowwise_pearson")) & perf.split.isin(order)]
        pivot = d.pivot(index="seed", columns="split", values="value").reindex(columns=order).dropna()
        sample = pivot.sample(n=min(35, len(pivot)), random_state=17)
        xx = np.arange(len(order))
        for _, row in sample.iterrows():
            ax.plot(xx, row.to_numpy(), color="#A9B7C5", lw=.45, alpha=.45)
        means, sd = pivot.mean(), pivot.std(ddof=1)
        ax.errorbar(xx, means, yerr=sd, fmt="o", color="#17202A", capsize=2.5, zorder=3, lw=1.0)
        ax.set_xticks(xx, xticks, fontsize=6.7); ax.set_ylabel("Mean row-wise Pearson")
        ax.set_title(dataset, loc="left", weight="bold")
        ax.set_ylim(-.03, max(.42, float(pivot.to_numpy().max()) + .06))
        ax.yaxis.grid(True, color=GRID, linewidth=.6)
        ax.set_axisbelow(True)
        panel_label(ax, "a" if dataset == "OpenProblems" else "b", -0.08, 1.04)
        out.append(d.assign(panel=dataset))
    contrast_keep = contrast[contrast.comparison.str.contains("matched", na=False) & contrast.model.isin(["ridge_cell_drug_fp", "control_conditioned_nonlinear_benchmark"])].copy()
    pd.concat(out, ignore_index=True).to_csv(SRC / "Figure_3_seed_level_source_data.csv", index=False)
    contrast_keep.to_csv(SRC / "Figure_3_paired_contrast_source_data.csv", index=False)
    save(fig, "Figure_3_train_size_matched_controls")


def fig4() -> None:
    perf = pd.read_csv(MASTER / "master_performance_summary.csv")
    configs = [
        ("OpenProblems", ["global_train_mean", "cell_context_mean", "drug_mean", "ridge_cell", "ridge_drug_fp", "ridge_cell_drug_fp"], .42),
        ("Sci-Plex 3", ["global_train_mean", "cell_context_mean", "drug_mean", "ridge_drug_fp", "ridge_cell_dose_drug_fp", "control_conditioned_nonlinear_benchmark"], .35),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.25))
    out = []
    for ax, (dataset, models, ymax) in zip(axes, configs):
        d = perf[(perf.dataset.eq(dataset)) & (perf.metric.eq("rowwise_pearson")) & perf.model.isin(models) & perf.split_label.isin(["random", "scaffold held-out", "joint held-out"])].copy()
        for model in models:
            m = d[d.model.eq(model)].set_index("split_label").reindex(["random", "scaffold held-out", "joint held-out"])
            if m.empty: continue
            ax.plot(["random", "scaffold", "joint"], m["mean"], marker="o", lw=1.15, markersize=3.2, label=label(model, dataset))
        ax.set_ylim(-.03, ymax); ax.set_ylabel("Mean row-wise Pearson")
        ax.set_title(dataset, loc="left", weight="bold")
        ax.yaxis.grid(True, color=GRID, linewidth=.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, ncol=1, loc="upper right", fontsize=5.5)
        panel_label(ax, "a" if dataset == "OpenProblems" else "b", -0.08, 1.04)
        out.append(d)
    pd.concat(out, ignore_index=True).to_csv(SRC / "Figure_4_multiprobe_source_data.csv", index=False)
    save(fig, "Figure_4_multiprobe_consistency")


def supp_loss() -> None:
    path = ROOT / "results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30/training_loss_by_epoch.csv"
    d = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for split, color in [("test_size_matched_random", "#277DA1"), ("train_test_size_matched_random", "#D1495B")]:
        x = d[d.split_type.eq(split)].groupby("epoch")[["train_loss", "val_loss"]].agg(["mean", "std"])
        epochs = x.index.to_numpy()
        for metric, style in [("train_loss", "-"), ("val_loss", "--")]:
            y = x[(metric, "mean")].to_numpy(); sd = x[(metric, "std")].to_numpy()
            ax.plot(epochs, y, color=color, linestyle=style, label=f"{split.replace('_', ' ')}: {metric.replace('_', ' ')}")
            ax.fill_between(epochs, y - sd, y + sd, color=color, alpha=.12)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Standardized MSE loss")
    ax.set_title("Control-conditioned nonlinear benchmark: matched-control training traces")
    ax.legend(frameon=False, fontsize=5.8, ncol=2)
    d.to_csv(SRC / "Supplementary_Figure_S1_loss_source_data.csv", index=False)
    save(fig, "Supplementary_Figure_S1_nonlinear_training_loss")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)
    fig1(); fig2(); fig3(); fig4(); supp_loss()


if __name__ == "__main__":
    main()
