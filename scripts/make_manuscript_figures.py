#!/usr/bin/env python3
"""Create manuscript-ready main figures for the leakage-resistant benchmark."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTDIR = Path("figures/manuscript_main")
SOURCE_DIR = OUTDIR / "source_data"

SPLIT_ORDER = [
    "split_random",
    "split_cell_heldout",
    "split_scaffold_heldout",
    "split_cell_scaffold_heldout",
]
SPLIT_LABELS = {
    "split_random": "Random",
    "split_cell_heldout": "Cell\nheld-out",
    "split_scaffold_heldout": "Scaffold\nheld-out",
    "split_cell_scaffold_heldout": "Cell + scaffold\nheld-out",
}
BASELINE_LABELS = {
    "global_train_mean": "Global mean",
    "ridge_cell_drug_fp": "Ridge: cell + drug FP",
    "nearest_drug_any_cell": "Nearest drug",
    "svd50_ridge_cell_drug_fp": "SVD50 + ridge",
}
COLORS = {
    "global_train_mean": "#7a7f87",
    "ridge_cell_drug_fp": "#2f6f73",
    "nearest_drug_any_cell": "#b56576",
    "svd50_ridge_cell_drug_fp": "#5b6ea6",
    "split_random": "#6f8fb4",
    "split_cell_heldout": "#72a88d",
    "split_scaffold_heldout": "#d39b5f",
    "split_cell_scaffold_heldout": "#9b6a8f",
}


def configure() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "legend.frameon": False,
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def add_box(ax: plt.Axes, xy: tuple[float, float], w: float, h: float, text: str, fc: str) -> None:
    patch = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=0.8,
        edgecolor="#2f3437",
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=7, color="#1f2427")


def add_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.8,
            color="#4c5358",
            shrinkA=2,
            shrinkB=2,
        )
    )


def make_figure1() -> None:
    fig = plt.figure(figsize=(7.2, 4.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.72, 1.0], width_ratios=[1.1, 0.9], hspace=0.30, wspace=0.34)
    ax_flow = fig.add_subplot(gs[0, :])
    ax_cov = fig.add_subplot(gs[1, 0])
    ax_split = fig.add_subplot(gs[1, 1])

    ax_flow.set_axis_off()
    panel_label(ax_flow, "a")
    boxes = [
        ((0.03, 0.42), 0.15, 0.28, "OpenProblems\nDGE matrix", "#e9eef3"),
        ((0.24, 0.42), 0.15, 0.28, "SMILES to\nscaffold", "#edf4ed"),
        ((0.45, 0.42), 0.15, 0.28, "Leakage-aware\nsplits", "#f4eee5"),
        ((0.66, 0.42), 0.15, 0.28, "Baselines\nand metrics", "#eeeaf3"),
        ((0.84, 0.42), 0.13, 0.28, "Similarity\naudit", "#f3ecec"),
    ]
    for xy, w, h, text, fc in boxes:
        add_box(ax_flow, xy, w, h, text, fc)
    for x0, x1 in [(0.18, 0.24), (0.39, 0.45), (0.60, 0.66), (0.81, 0.84)]:
        add_arrow(ax_flow, (x0, 0.56), (x1, 0.56))
    ax_flow.text(
        0.50,
        0.18,
        "Core contrast: interpolation-like random validation versus unseen scaffold and unseen cell-context validation",
        ha="center",
        va="center",
        fontsize=7,
        color="#4c5358",
    )
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)

    panel_label(ax_cov, "b")
    cov = pd.read_csv("results/feasibility_openproblems_full/feasibility_summary.csv").iloc[0]
    labels = ["Records", "Treated", "Drugs", "Valid\nscaffolds", "Cell\ncontexts"]
    values = [
        cov["n_records"],
        545,
        cov["n_drugs_in_perturbation_metadata"],
        cov["n_valid_scaffolds"],
        cov["n_cell_contexts"],
    ]
    bars = ax_cov.bar(labels, values, color=["#8aa4bd", "#8aa4bd", "#75a38b", "#d39b5f", "#9b6a8f"], width=0.68)
    ax_cov.set_ylabel("Count")
    ax_cov.set_title("Benchmark coverage", loc="left", fontsize=8, pad=3)
    for b, v in zip(bars, values):
        ax_cov.text(b.get_x() + b.get_width() / 2, b.get_height() + max(values) * 0.02, f"{int(v)}", ha="center", va="bottom", fontsize=6.5)
    ax_cov.set_ylim(0, max(values) * 1.22)

    panel_label(ax_split, "c")
    split = pd.read_csv("results/openproblems_multiseed_10_de/all_split_summary.csv")
    summary = split.groupby("split")[["train_records", "test_records", "excluded_records"]].mean().loc[SPLIT_ORDER]
    x = np.arange(len(summary))
    bottom = np.zeros(len(summary))
    stack_cols = [("train_records", "Train", "#6f8fb4"), ("test_records", "Test", "#d39b5f"), ("excluded_records", "Excluded", "#c8c8c8")]
    for col, lab, color in stack_cols:
        ax_split.bar(x, summary[col], bottom=bottom, color=color, width=0.72, label=lab)
        bottom += summary[col].to_numpy()
    ax_split.set_xticks(x, [SPLIT_LABELS[s] for s in SPLIT_ORDER], rotation=0)
    ax_split.set_ylabel("Mean records")
    ax_split.set_title("Repeated split sizes", loc="left", fontsize=8, pad=3)
    ax_split.legend(loc="upper left", bbox_to_anchor=(1.01, 1.02), ncols=1, handlelength=1.0)

    fig.suptitle("Leakage-resistant benchmark design", x=0.02, y=0.995, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    save_figure(fig, "figure1_benchmark_design")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"label": labels, "count": values}).to_csv(SOURCE_DIR / "figure1b_coverage.csv", index=False)
    summary.reset_index().to_csv(SOURCE_DIR / "figure1c_split_sizes.csv", index=False)


def plot_metric_panel(ax: plt.Axes, df: pd.DataFrame, metric: str, ylabel: str, baselines: list[str], title: str) -> None:
    x = np.arange(len(SPLIT_ORDER))
    offsets = np.linspace(-0.27, 0.27, len(baselines))
    for off, baseline in zip(offsets, baselines):
        sub = (
            df[df["baseline"].eq(baseline)]
            .groupby("split")[metric]
            .agg(mean="mean", sd="std")
            .reindex(SPLIT_ORDER)
        )
        y = sub["mean"].to_numpy()
        yerr = sub["sd"].to_numpy()
        ax.errorbar(
            x + off,
            y,
            yerr=yerr,
            marker="o",
            markersize=3.2,
            linewidth=1.2,
            capsize=2,
            color=COLORS.get(baseline, "#555555"),
            label=BASELINE_LABELS.get(baseline, baseline),
        )
    ax.set_xticks(x, [SPLIT_LABELS[s] for s in SPLIT_ORDER])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=8, pad=3)
    ax.grid(axis="y", color="#e8e8e8", linewidth=0.6)


def make_figure2() -> None:
    df = pd.read_csv("results/openproblems_multiseed_10_de/all_baseline_with_nearest_summary.csv")
    keep = ["global_train_mean", "ridge_cell_drug_fp", "nearest_drug_any_cell"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    panel_label(axes[0], "a")
    plot_metric_panel(axes[0], df, "mean_rowwise_pearson", "Mean row-wise Pearson", keep, "All genes")
    panel_label(axes[1], "b")
    plot_metric_panel(axes[1], df, "mean_de_rowwise_pearson", "DE-gene row-wise Pearson", keep, "DE genes")
    axes[1].legend(loc="upper right", bbox_to_anchor=(1.02, 1.02))
    fig.suptitle("Baseline performance declines under stricter validation", x=0.01, y=1.03, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "figure2_baseline_decay")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    df[df["baseline"].isin(keep)].to_csv(SOURCE_DIR / "figure2_baseline_decay.csv", index=False)


def make_figure3() -> None:
    df = pd.read_csv("results/topk_de_overlap_multiseed_10/all_topk_de_overlap_summary.csv")
    keep = ["global_train_mean", "ridge_cell_drug_fp", "nearest_drug_any_cell"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True)
    panel_specs = [
        (axes[0, 0], 100, "mean_topk_overlap", "Top-100 overlap", "a"),
        (axes[0, 1], 100, "mean_direction_agreement", "Top-100 direction", "b"),
        (axes[1, 0], 200, "mean_topk_overlap", "Top-200 overlap", "c"),
        (axes[1, 1], 200, "mean_direction_agreement", "Top-200 direction", "d"),
    ]
    for ax, k, metric, title, lab in panel_specs:
        panel_label(ax, lab)
        subdf = df[df["top_k"].eq(k)].copy()
        plot_metric_panel(ax, subdf, metric, title, keep, title)
        if "direction" in metric:
            ax.set_ylim(0.55, 1.02)
        else:
            ax.set_ylim(0, 0.20)
    axes[0, 1].legend(loc="lower left", bbox_to_anchor=(0.02, 0.02))
    fig.suptitle("Response-gene recovery deteriorates in strict held-out tasks", x=0.01, y=1.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "figure3_topk_response_gene_recovery")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    df[df["baseline"].isin(keep)].to_csv(SOURCE_DIR / "figure3_topk_response_gene_recovery.csv", index=False)


def make_figure4() -> None:
    ci = pd.read_csv("results/drug_similarity_bootstrap_ci/similarity_metric_seed_bootstrap_ci.csv")
    ci = ci[ci["metric"].eq("mean_max_train_tanimoto")].set_index("split").loc[SPLIT_ORDER].reset_index()
    contrast = pd.read_csv("results/drug_similarity_bootstrap_ci/similarity_metric_seed_bootstrap_contrasts.csv")
    contrast = contrast[contrast["metric"].eq("mean_max_train_tanimoto")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.25, 0.9]})
    panel_label(axes[0], "a")
    x = np.arange(len(ci))
    axes[0].bar(x, ci["mean"], color=[COLORS[s] for s in ci["split"]], width=0.65)
    yerr = np.vstack([ci["mean"] - ci["ci_low"], ci["ci_high"] - ci["mean"]])
    axes[0].errorbar(x, ci["mean"], yerr=yerr, fmt="none", ecolor="#2f3437", capsize=2, linewidth=0.8)
    axes[0].set_xticks(x, [SPLIT_LABELS[s] for s in ci["split"]])
    axes[0].set_ylabel("Nearest train-drug Tanimoto")
    axes[0].set_ylim(0, 1.08)
    axes[0].set_title("Chemical-neighbor availability", loc="left", fontsize=8, pad=3)
    axes[0].grid(axis="y", color="#e8e8e8", linewidth=0.6)

    panel_label(axes[1], "b")
    labels = ["Random -\nscaffold", "Random -\njoint"]
    y = contrast["mean_difference"].to_numpy()
    yerr = np.vstack([y - contrast["ci_low"].to_numpy(), contrast["ci_high"].to_numpy() - y])
    axes[1].bar(np.arange(len(y)), y, color=["#d39b5f", "#9b6a8f"], width=0.62)
    axes[1].errorbar(np.arange(len(y)), y, yerr=yerr, fmt="none", ecolor="#2f3437", capsize=2, linewidth=0.8)
    axes[1].axhline(0, color="#333333", linewidth=0.7)
    axes[1].set_xticks(np.arange(len(y)), labels)
    axes[1].set_ylabel("Paired difference")
    axes[1].set_ylim(0, 0.82)
    axes[1].set_title("Seed-level bootstrap contrasts", loc="left", fontsize=8, pad=3)
    axes[1].grid(axis="y", color="#e8e8e8", linewidth=0.6)
    fig.suptitle("Scaffold-strict splits remove chemical-neighbor shortcuts", x=0.01, y=1.03, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "figure4_chemical_similarity_audit")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    ci.to_csv(SOURCE_DIR / "figure4a_similarity_ci.csv", index=False)
    contrast.to_csv(SOURCE_DIR / "figure4b_similarity_contrasts.csv", index=False)


def make_figure5() -> None:
    baseline = pd.read_csv("results/openproblems_multiseed_10_de/all_baseline_with_nearest_summary.csv")
    svd = pd.read_csv("results/svd_ridge_multiseed_10/all_svd_ridge_baseline_summary.csv")
    df = pd.concat([baseline, svd], ignore_index=True, sort=False)
    keep = ["global_train_mean", "ridge_cell_drug_fp", "nearest_drug_any_cell", "svd50_ridge_cell_drug_fp"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharex=True)
    panel_label(axes[0], "a")
    plot_metric_panel(axes[0], df, "mean_rowwise_pearson", "Mean row-wise Pearson", keep, "All genes")
    panel_label(axes[1], "b")
    plot_metric_panel(axes[1], df, "mean_de_rowwise_pearson", "DE-gene row-wise Pearson", keep, "DE genes")
    axes[1].legend(loc="upper right", bbox_to_anchor=(1.02, 1.02))
    fig.suptitle("Low-rank SVD-ridge improves easier splits but not joint generalization", x=0.01, y=1.03, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "figure5_svd_ridge_stress_test")
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    df[df["baseline"].isin(keep)].to_csv(SOURCE_DIR / "figure5_svd_ridge_stress_test.csv", index=False)


def main() -> None:
    configure()
    make_figure1()
    make_figure2()
    make_figure3()
    make_figure4()
    make_figure5()
    print(f"Wrote figures to {OUTDIR}")


if __name__ == "__main__":
    main()
