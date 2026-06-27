#!/usr/bin/env python3
"""Pathway-level response recovery audit for the JGG revision.

This script uses a local MSigDB Hallmark GMT file only. It reconstructs
available benchmark-compatible prediction vectors, projects true and predicted
response vectors onto Hallmark gene sets, and writes source data, summary
tables, figures, and a short feasibility/check report.
"""

from __future__ import annotations

import csv
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
GMT = ROOT / "resources/gene_sets/hallmark_symbols.gmt"
OUTDIR = ROOT / "results/pathway_level"
FIGDIR = ROOT / "figures"
SRCDIR = ROOT / "figures/source_data"
DOCSDIR = ROOT / "docs"
MIN_OVERLAP = 10
ALPHA = 10.0
N_BITS = 1024
RNG = np.random.default_rng(20260625)

OPENPROBLEMS_SEEDS = list(range(1, 11))
SCIPLEX_SEEDS = list(range(1, 31))
SPLITS = {
    "split_random": "random",
    "split_scaffold_heldout": "scaffold_heldout",
    "split_cell_scaffold_heldout": "cell_scaffold_heldout",
}
SCIPLEX_SPLITS = {
    "random": "random",
    "scaffold_heldout": "scaffold_heldout",
    "cell_scaffold_heldout": "cell_scaffold_heldout",
}
MODEL_LABELS = {
    "global_train_mean": "Global mean",
    "ridge_cell_drug_fp": "Ridge cell+drug FP",
    "ridge_drug_fp": "Ridge drug FP",
    "ridge_cell_dose_drug_fp": "Ridge cell+dose+drug FP",
    "transigen_adapted_sciplex3": "TranSiGen-style",
    "prnet_adapted_sciplex3": "PRnet-style",
}
REPRESENTATIVE_SETS = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
]


@dataclass(frozen=True)
class GeneSet:
    name: str
    description: str
    genes: tuple[str, ...]


def ensure_dirs() -> None:
    for p in [OUTDIR, FIGDIR, SRCDIR, DOCSDIR]:
        p.mkdir(parents=True, exist_ok=True)


def parse_gmt(path: Path) -> list[GeneSet]:
    if not path.exists():
        raise FileNotFoundError(path)
    gene_sets: list[GeneSet] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise ValueError(f"Invalid GMT row {line_no}: expected at least 3 tab-delimited fields.")
            name, desc, *genes = parts
            genes = [g.strip().upper() for g in genes if g.strip()]
            if not name or not genes:
                raise ValueError(f"Invalid GMT row {line_no}: missing name or genes.")
            gene_sets.append(GeneSet(name=name.strip(), description=desc.strip(), genes=tuple(dict.fromkeys(genes))))
    if not gene_sets:
        raise ValueError("No gene sets parsed from GMT.")
    return gene_sets


def gene_set_indices(gene_sets: list[GeneSet], gene_names: list[str], dataset: str) -> tuple[list[dict], list[tuple[GeneSet, np.ndarray]]]:
    index = {g.upper(): i for i, g in enumerate(gene_names)}
    rows: list[dict] = []
    usable: list[tuple[GeneSet, np.ndarray]] = []
    for gs in gene_sets:
        overlap = [g for g in gs.genes if g in index]
        status = "usable" if len(overlap) >= MIN_OVERLAP else "low_overlap"
        rows.append(
            {
                "dataset": dataset,
                "gene_set_name": gs.name,
                "n_genes_in_set": len(gs.genes),
                "n_genes_overlap": len(overlap),
                "min_overlap": MIN_OVERLAP,
                "status": status,
            }
        )
        if status == "usable":
            usable.append((gs, np.array([index[g] for g in overlap], dtype=np.int64)))
    return rows, usable


def pathway_scores(y: np.ndarray, usable_sets: list[tuple[GeneSet, np.ndarray]]) -> np.ndarray:
    scores = np.empty((y.shape[0], len(usable_sets)), dtype=np.float64)
    for j, (_, idx) in enumerate(usable_sets):
        scores[:, j] = y[:, idx].mean(axis=1)
    return scores


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ac = a - a.mean(axis=1, keepdims=True)
    bc = b - b.mean(axis=1, keepdims=True)
    num = (ac * bc).sum(axis=1)
    den = np.sqrt((ac**2).sum(axis=1) * (bc**2).sum(axis=1))
    return np.divide(num, den, out=np.full(a.shape[0], np.nan), where=den != 0)


def rmse(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sqrt(((a - b) ** 2).mean(axis=1))


def ridge_predict(train_x: np.ndarray, test_x: np.ndarray, train_y: np.ndarray, alpha: float = ALPHA) -> np.ndarray:
    train_design = np.column_stack([np.ones(train_x.shape[0], dtype=np.float64), train_x.astype(np.float64)])
    test_design = np.column_stack([np.ones(test_x.shape[0], dtype=np.float64), test_x.astype(np.float64)])
    penalty = np.eye(train_design.shape[1], dtype=np.float64) * alpha
    penalty[0, 0] = 0.0
    coef = np.linalg.solve(train_design.T @ train_design + penalty, train_design.T @ train_y)
    return test_design @ coef


def one_hot_train_test(train_values: pd.Series, test_values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    levels = sorted(map(str, pd.unique(train_values)))
    pos = {v: i for i, v in enumerate(levels)}
    train = np.zeros((len(train_values), len(levels)), dtype=np.float64)
    test = np.zeros((len(test_values), len(levels)), dtype=np.float64)
    for i, v in enumerate(map(str, train_values)):
        train[i, pos[v]] = 1.0
    for i, v in enumerate(map(str, test_values)):
        if v in pos:
            test[i, pos[v]] = 1.0
    return train, test


def morgan_fingerprint(smiles: str, n_bits: int = N_BITS) -> np.ndarray:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    arr = np.zeros((n_bits,), dtype=np.float64)
    if not isinstance(smiles, str) or not smiles:
        return arr
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return arr
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    arr[list(fp.GetOnBits())] = 1.0
    return arr


def append_run(
    metrics_path: Path,
    seed_summary_rows: list[dict],
    gene_error_rows: list[dict],
    dataset: str,
    model: str,
    split: str,
    seed: int,
    record_ids: np.ndarray | list[str],
    true_y: np.ndarray,
    pred_y: np.ndarray,
    usable_sets: list[tuple[GeneSet, np.ndarray]],
    notes: str = "",
) -> None:
    if true_y.shape != pred_y.shape:
        raise ValueError(f"{dataset}/{model}/{split}/seed{seed}: true/pred shape mismatch {true_y.shape} vs {pred_y.shape}")
    true_scores = pathway_scores(true_y, usable_sets)
    pred_scores = pathway_scores(pred_y, usable_sets)
    pearson = rowwise_pearson(true_scores, pred_scores)
    record_rmse = rmse(true_scores, pred_scores)
    abs_error = np.abs(true_scores - pred_scores)

    seed_summary_rows.append(
        {
            "dataset": dataset,
            "model": model,
            "model_label": MODEL_LABELS.get(model, model),
            "split": split,
            "seed": seed,
            "n_records": int(true_y.shape[0]),
            "n_gene_sets": int(len(usable_sets)),
            "mean_pathway_pearson": float(np.nanmean(pearson)),
            "sd_pathway_pearson": float(np.nanstd(pearson, ddof=1)) if np.isfinite(pearson).sum() > 1 else np.nan,
            "median_pathway_pearson": float(np.nanmedian(pearson)),
            "iqr_pathway_pearson": float(np.nanpercentile(pearson, 75) - np.nanpercentile(pearson, 25)),
            "mean_pathway_rmse": float(np.nanmean(record_rmse)),
            "median_pathway_rmse": float(np.nanmedian(record_rmse)),
            "status": "completed",
            "notes": notes,
        }
    )

    long_rows = []
    record_ids = list(map(str, record_ids))
    for j, (gs, idx) in enumerate(usable_sets):
        ae = abs_error[:, j]
        gene_error_rows.append(
            {
                "dataset": dataset,
                "model": model,
                "model_label": MODEL_LABELS.get(model, model),
                "split": split,
                "seed": seed,
                "gene_set_name": gs.name,
                "n_genes_in_set": len(gs.genes),
                "n_genes_overlap": int(len(idx)),
                "mean_abs_error": float(np.nanmean(ae)),
                "median_abs_error": float(np.nanmedian(ae)),
                "mean_true_score": float(np.nanmean(true_scores[:, j])),
                "mean_pred_score": float(np.nanmean(pred_scores[:, j])),
                "status": "completed",
            }
        )
        for i, rid in enumerate(record_ids):
            long_rows.append(
                {
                    "dataset": dataset,
                    "model": model,
                    "model_label": MODEL_LABELS.get(model, model),
                    "split": split,
                    "seed": seed,
                    "record_id": rid,
                    "gene_set_name": gs.name,
                    "n_genes_in_set": len(gs.genes),
                    "n_genes_overlap": int(len(idx)),
                    "pathway_true_score": true_scores[i, j],
                    "pathway_pred_score": pred_scores[i, j],
                    "pathway_abs_error": ae[i],
                    "pathway_pearson": pearson[i],
                    "pathway_rmse": record_rmse[i],
                    "status": "completed",
                    "notes": notes,
                }
            )
    chunk = pd.DataFrame(long_rows)
    chunk.to_csv(metrics_path, mode="a", header=not metrics_path.exists(), index=False)


def load_openproblems_matrix() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    import anndata as ad

    metas = []
    ys = []
    genes: list[str] | None = None
    for path, source_split in [
        (ROOT / "data/raw/openproblems_neurips2023/de_train.h5ad", "de_train"),
        (ROOT / "data/raw/openproblems_neurips2023/de_test.h5ad", "de_test"),
    ]:
        adata = ad.read_h5ad(path)
        if genes is None:
            genes = list(map(str, adata.var_names))
        meta = adata.obs.copy().reset_index(names="sample_id")
        meta["source_split"] = source_split
        metas.append(meta)
        ys.append(np.asarray(adata.layers["clipped_sign_log10_pval"], dtype=np.float64))
    return pd.concat(metas, ignore_index=True), np.vstack(ys), genes or []


def run_openproblems(gene_sets: list[GeneSet], metrics_path: Path) -> tuple[list[dict], list[dict], list[dict], list[str]]:
    meta0, y0, genes = load_openproblems_matrix()
    overlap_rows, usable_sets = gene_set_indices(gene_sets, genes, "openproblems_neurips2023")
    if len(usable_sets) < 3:
        return overlap_rows, [], [], [f"OpenProblems usable Hallmark sets below threshold: {len(usable_sets)}"]

    seed_summary_rows: list[dict] = []
    gene_error_rows: list[dict] = []
    notes: list[str] = []
    for seed in OPENPROBLEMS_SEEDS:
        split_path = ROOT / f"results/openproblems_multiseed_10_de/seed_{seed:03d}/splits/candidate_splits.csv"
        if not split_path.exists():
            notes.append(f"OpenProblems seed {seed}: split manifest missing; skipped.")
            continue
        splits = pd.read_csv(split_path)
        meta = meta0.merge(
            splits[["sample_id", "drug_name", "cell_context", "condition", "smiles", *SPLITS.keys()]],
            on="sample_id",
            how="left",
            validate="one_to_one",
        )
        treated = meta["condition"].eq("treated").to_numpy()
        meta = meta.loc[treated].reset_index(drop=True)
        y = y0[treated]
        for split_col, split_name in SPLITS.items():
            train_mask = meta[split_col].eq("train").to_numpy()
            test_mask = meta[split_col].eq("test").to_numpy()
            if train_mask.sum() == 0 or test_mask.sum() == 0:
                notes.append(f"OpenProblems seed {seed} {split_name}: empty train/test; skipped.")
                continue
            train_meta = meta.loc[train_mask].reset_index(drop=True)
            test_meta = meta.loc[test_mask].reset_index(drop=True)
            train_y = y[train_mask]
            test_y = y[test_mask]
            record_ids = test_meta["sample_id"].astype(str).to_numpy()

            mean_pred = np.repeat(train_y.mean(axis=0, keepdims=True), len(test_meta), axis=0)
            append_run(
                metrics_path,
                seed_summary_rows,
                gene_error_rows,
                "openproblems_neurips2023",
                "global_train_mean",
                split_name,
                seed,
                record_ids,
                test_y,
                mean_pred,
                usable_sets,
                "Reconstructed from train-only mean response vectors.",
            )

            cell_train, cell_test = one_hot_train_test(train_meta["cell_context"], test_meta["cell_context"])
            train_fp = np.vstack([morgan_fingerprint(s) for s in train_meta["smiles"]])
            test_fp = np.vstack([morgan_fingerprint(s) for s in test_meta["smiles"]])
            pred = ridge_predict(np.hstack([cell_train, train_fp]), np.hstack([cell_test, test_fp]), train_y)
            append_run(
                metrics_path,
                seed_summary_rows,
                gene_error_rows,
                "openproblems_neurips2023",
                "ridge_cell_drug_fp",
                split_name,
                seed,
                record_ids,
                test_y,
                pred,
                usable_sets,
                "Reconstructed with alpha=10 using train-only cell-context one-hot and Morgan fingerprints.",
            )
    notes.append("OpenProblems nearest-drug and SVD-ridge were not projected because only per-row metrics, not full prediction vectors, were available locally.")
    return overlap_rows, seed_summary_rows, gene_error_rows, notes


def load_sciplex() -> tuple[pd.DataFrame, np.ndarray, list[str], dict[str, np.ndarray]]:
    matrices = np.load(ROOT / "data/processed/transigen/sciplex3/matrices.npz", allow_pickle=True)
    response = np.asarray(matrices["response"], dtype=np.float64)
    genes = list(map(str, matrices["gene"]))
    record_ids = list(map(str, matrices["record_id"]))
    meta = pd.read_csv(ROOT / "data/processed/transigen/sciplex3/record_metadata.csv")
    meta = meta.set_index("record_id").loc[record_ids].reset_index()
    fp = np.load(ROOT / "data/processed/transigen/sciplex3/compound_fingerprints.npz", allow_pickle=True)
    fps = {str(d): np.asarray(v, dtype=np.float64) for d, v in zip(fp["drug_name"], fp["morgan1024"])}
    return meta, response, genes, fps


def sciplex_feature_blocks(
    train_meta: pd.DataFrame,
    test_meta: pd.DataFrame,
    fps: dict[str, np.ndarray],
    feature_set: str,
) -> tuple[np.ndarray, np.ndarray]:
    train_blocks = []
    test_blocks = []
    if "cell" in feature_set:
        tr, te = one_hot_train_test(train_meta["cell_line"], test_meta["cell_line"])
        train_blocks.append(tr)
        test_blocks.append(te)
    if "dose" in feature_set:
        train_dose = np.log10(train_meta[["dose_value"]].astype(float).to_numpy())
        test_dose = np.log10(test_meta[["dose_value"]].astype(float).to_numpy())
        mu = train_dose.mean(axis=0)
        sd = train_dose.std(axis=0)
        sd[sd == 0] = 1.0
        train_blocks.append((train_dose - mu) / sd)
        test_blocks.append((test_dose - mu) / sd)
    if "drug_fp" in feature_set:
        train_blocks.append(np.vstack([fps.get(str(d), np.zeros(N_BITS)) for d in train_meta["drug_name"]]))
        test_blocks.append(np.vstack([fps.get(str(d), np.zeros(N_BITS)) for d in test_meta["drug_name"]]))
    return np.hstack(train_blocks), np.hstack(test_blocks)


def npz_prediction_runs(manifest_path: Path, model_name: str) -> dict[tuple[int, str], Path]:
    if not manifest_path.exists():
        return {}
    df = pd.read_csv(manifest_path)
    df = df[df["status"].eq("completed")]
    runs = {}
    for (seed, split), g in df.groupby(["seed", "split_type"], observed=True):
        path = Path(str(g["pred_response_vector_path"].iloc[0]))
        if path.exists():
            runs[(int(seed), str(split))] = path
    return runs


def run_sciplex(gene_sets: list[GeneSet], metrics_path: Path) -> tuple[list[dict], list[dict], list[dict], list[str]]:
    meta, y, genes, fps = load_sciplex()
    overlap_rows, usable_sets = gene_set_indices(gene_sets, genes, "sciplex3_24h_top2000")
    if len(usable_sets) < 3:
        return overlap_rows, [], [], [f"Sci-Plex 3 usable Hallmark sets below threshold: {len(usable_sets)}"]

    split_assign = pd.read_csv(ROOT / "data/processed/transigen/sciplex3/split_assignments_long.csv")
    id_to_pos = {rid: i for i, rid in enumerate(meta["record_id"].astype(str))}
    seed_summary_rows: list[dict] = []
    gene_error_rows: list[dict] = []
    notes: list[str] = []
    deep_runs = {
        "transigen_adapted_sciplex3": npz_prediction_runs(
            ROOT / "results/deep_model_panel/sciplex3/transigen_main30_predictions_manifest.csv",
            "transigen_adapted_sciplex3",
        ),
        "prnet_adapted_sciplex3": npz_prediction_runs(
            ROOT / "results/deep_model_panel/sciplex3/prnet_main30_predictions_manifest.csv",
            "prnet_adapted_sciplex3",
        ),
    }

    for seed in SCIPLEX_SEEDS:
        seed_split = split_assign[split_assign["seed"].eq(seed)]
        for split_source, split_name in SCIPLEX_SPLITS.items():
            g = seed_split[seed_split["split_name"].eq(split_source)]
            train_ids = g.loc[g["assignment"].eq("train"), "record_id"].astype(str).tolist()
            test_ids = g.loc[g["assignment"].eq("test"), "record_id"].astype(str).tolist()
            train_idx = np.array([id_to_pos[r] for r in train_ids if r in id_to_pos], dtype=np.int64)
            test_idx = np.array([id_to_pos[r] for r in test_ids if r in id_to_pos], dtype=np.int64)
            if len(train_idx) == 0 or len(test_idx) == 0:
                notes.append(f"Sci-Plex 3 seed {seed} {split_name}: empty train/test; skipped.")
                continue
            train_meta = meta.iloc[train_idx].reset_index(drop=True)
            test_meta = meta.iloc[test_idx].reset_index(drop=True)
            train_y = y[train_idx]
            test_y = y[test_idx]
            record_ids = test_meta["record_id"].astype(str).to_numpy()

            mean_pred = np.repeat(train_y.mean(axis=0, keepdims=True), len(test_meta), axis=0)
            append_run(
                metrics_path,
                seed_summary_rows,
                gene_error_rows,
                "sciplex3_24h_top2000",
                "global_train_mean",
                split_name,
                seed,
                record_ids,
                test_y,
                mean_pred,
                usable_sets,
                "Reconstructed from train-only mean response vectors.",
            )

            for model, feature_set in [
                ("ridge_drug_fp", "drug_fp"),
                ("ridge_cell_dose_drug_fp", "cell+dose+drug_fp"),
            ]:
                train_x, test_x = sciplex_feature_blocks(train_meta, test_meta, fps, feature_set)
                pred = ridge_predict(train_x, test_x, train_y)
                append_run(
                    metrics_path,
                    seed_summary_rows,
                    gene_error_rows,
                    "sciplex3_24h_top2000",
                    model,
                    split_name,
                    seed,
                    record_ids,
                    test_y,
                    pred,
                    usable_sets,
                    "Reconstructed with alpha=10 and training-only feature preprocessing.",
                )

            for model, runs in deep_runs.items():
                path = runs.get((seed, split_source))
                if path is None:
                    notes.append(f"{model} seed {seed} {split_name}: completed NPZ missing; skipped.")
                    continue
                data = np.load(path, allow_pickle=True)
                npz_ids = data["record_id"].astype(str)
                pred = np.asarray(data["prediction"], dtype=np.float64)
                target = np.asarray(data["target"], dtype=np.float64)
                if pred.shape[1] != 2000 or target.shape[1] != 2000:
                    raise ValueError(f"{path}: expected 2000 genes, got pred {pred.shape}, target {target.shape}.")
                if list(npz_ids) != list(record_ids):
                    raise ValueError(f"{path}: NPZ record_id order does not match split test record_id order.")
                append_run(
                    metrics_path,
                    seed_summary_rows,
                    gene_error_rows,
                    "sciplex3_24h_top2000",
                    model,
                    split_name,
                    seed,
                    record_ids,
                    target,
                    pred,
                    usable_sets,
                    f"Loaded completed prediction vectors from {path.name}.",
                )
    return overlap_rows, seed_summary_rows, gene_error_rows, notes


def bootstrap_ci(values: np.ndarray, n_boot: int = 10000) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    if len(values) == 1:
        return float(values[0]), float(values[0]), float(values[0])
    idx = RNG.integers(0, len(values), size=(n_boot, len(values)))
    means = values[idx].mean(axis=1)
    return float(values.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def summarize_outputs(seed_summary: pd.DataFrame, per_gene: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    contrast_rows = []
    ci_rows = []
    for (dataset, model), g in seed_summary.groupby(["dataset", "model"], observed=True):
        for split, sg in g.groupby("split", observed=True):
            vals = sg["mean_pathway_pearson"].to_numpy()
            mean, lo, hi = bootstrap_ci(vals)
            ci_rows.append(
                {
                    "dataset": dataset,
                    "model": model,
                    "model_label": MODEL_LABELS.get(model, model),
                    "quantity": "mean_pathway_pearson",
                    "split": split,
                    "contrast": "",
                    "n_seed_units": int(len(vals)),
                    "mean": mean,
                    "ci95_low": lo,
                    "ci95_high": hi,
                }
            )
        pivot = g.pivot_table(index="seed", columns="split", values="mean_pathway_pearson", aggfunc="mean")
        for strict in ["scaffold_heldout", "cell_scaffold_heldout"]:
            if {"random", strict}.issubset(pivot.columns):
                diff = (pivot["random"] - pivot[strict]).dropna().to_numpy()
                mean, lo, hi = bootstrap_ci(diff)
                contrast_rows.append(
                    {
                        "dataset": dataset,
                        "model": model,
                        "model_label": MODEL_LABELS.get(model, model),
                        "contrast": f"random_minus_{strict}",
                        "metric": "mean_pathway_pearson",
                        "n_seed_units": int(len(diff)),
                        "mean_difference": mean,
                        "ci95_low": lo,
                        "ci95_high": hi,
                    }
                )
                ci_rows.append(
                    {
                        "dataset": dataset,
                        "model": model,
                        "model_label": MODEL_LABELS.get(model, model),
                        "quantity": "paired_difference",
                        "split": "",
                        "contrast": f"random_minus_{strict}",
                        "n_seed_units": int(len(diff)),
                        "mean": mean,
                        "ci95_low": lo,
                        "ci95_high": hi,
                    }
                )
    return pd.DataFrame(contrast_rows), pd.DataFrame(ci_rows)


def add_summary_stats(seed_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, model, split), g in seed_summary.groupby(["dataset", "model", "split"], observed=True):
        vals = g["mean_pathway_pearson"].to_numpy(dtype=np.float64)
        rows.append(
            {
                "dataset": dataset,
                "model": model,
                "model_label": MODEL_LABELS.get(model, model),
                "split": split,
                "n_seed_units": int(len(vals)),
                "mean": float(np.nanmean(vals)),
                "sd": float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else np.nan,
                "median": float(np.nanmedian(vals)),
                "iqr": float(np.nanpercentile(vals, 75) - np.nanpercentile(vals, 25)),
                "metric": "mean_pathway_pearson",
            }
        )
    return pd.DataFrame(rows)


def write_figures(seed_summary: pd.DataFrame, contrasts: pd.DataFrame, per_gene: pd.DataFrame, overlap: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
        }
    )
    summary = add_summary_stats(seed_summary)
    source_rows = []
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.0))
    datasets = [("openproblems_neurips2023", "OpenProblems"), ("sciplex3_24h_top2000", "Sci-Plex 3")]
    colors = {"random": "#3C78D8", "scaffold_heldout": "#6AA84F", "cell_scaffold_heldout": "#CC4125"}
    for ax, (dataset, title) in zip(axes[0], datasets):
        sub = summary[summary["dataset"].eq(dataset)].copy()
        if sub.empty:
            ax.text(0.5, 0.5, "No vector-level runs", ha="center", va="center")
            ax.set_axis_off()
            continue
        models = list(dict.fromkeys(sub["model_label"]))
        width = 0.23
        x = np.arange(len(models))
        for k, split in enumerate(["random", "scaffold_heldout", "cell_scaffold_heldout"]):
            vals = []
            errs = []
            for model in models:
                row = sub[sub["model_label"].eq(model) & sub["split"].eq(split)]
                vals.append(float(row["mean"].iloc[0]) if not row.empty else np.nan)
                errs.append(float(row["sd"].iloc[0]) if not row.empty else 0.0)
                if not row.empty:
                    d = row.iloc[0].to_dict()
                    d["panel"] = "A" if dataset.startswith("open") else "B"
                    source_rows.append(d)
            ax.bar(x + (k - 1) * width, vals, width, yerr=errs, color=colors[split], label=split.replace("_", " "), capsize=2)
        ax.set_title(title)
        ax.set_ylabel("Pathway-level Pearson")
        ax.set_xticks(x)
        ax.set_xticklabels(models, rotation=35, ha="right")
        ax.axhline(0, color="#555555", linewidth=0.6)
        ax.legend(frameon=False, fontsize=7)

    ax = axes[1, 0]
    csub = contrasts[contrasts["contrast"].eq("random_minus_cell_scaffold_heldout")].copy()
    if not csub.empty:
        labels = [f"{r.dataset.replace('_24h_top2000','').replace('_neurips2023','')} | {MODEL_LABELS.get(r.model, r.model)}" for r in csub.itertuples()]
        y = np.arange(len(labels))
        ax.barh(y, csub["mean_difference"], color="#674EA7")
        ax.errorbar(
            csub["mean_difference"],
            y,
            xerr=[csub["mean_difference"] - csub["ci95_low"], csub["ci95_high"] - csub["mean_difference"]],
            fmt="none",
            ecolor="#222222",
            capsize=2,
            linewidth=0.8,
        )
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Random minus joint pathway Pearson")
        ax.set_title("Random-to-joint deterioration")
        for r in csub.to_dict("records"):
            r["panel"] = "C"
            source_rows.append(r)
    else:
        ax.text(0.5, 0.5, "No paired contrasts", ha="center", va="center")
        ax.set_axis_off()

    ax = axes[1, 1]
    heat = (
        per_gene[
            per_gene["dataset"].eq("sciplex3_24h_top2000")
            & per_gene["split"].eq("cell_scaffold_heldout")
            & per_gene["gene_set_name"].isin(REPRESENTATIVE_SETS)
            & per_gene["model"].isin(["ridge_cell_dose_drug_fp", "transigen_adapted_sciplex3", "prnet_adapted_sciplex3"])
        ]
        .groupby(["gene_set_name", "model_label"], observed=True)["mean_abs_error"]
        .mean()
        .reset_index()
    )
    if not heat.empty:
        matrix = heat.pivot(index="gene_set_name", columns="model_label", values="mean_abs_error")
        matrix = matrix.reindex([g for g in REPRESENTATIVE_SETS if g in matrix.index])
        im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="viridis")
        ax.set_xticks(np.arange(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=35, ha="right")
        ax.set_yticks(np.arange(matrix.shape[0]))
        ax.set_yticklabels([s.replace("HALLMARK_", "").replace("_", " ") for s in matrix.index], fontsize=7)
        ax.set_title("Joint-split Hallmark absolute error")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("Mean absolute error")
        tmp = heat.copy()
        tmp["panel"] = "D"
        source_rows.extend(tmp.to_dict("records"))
    else:
        ax.text(0.5, 0.5, "Representative sets unavailable", ha="center", va="center")
        ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(FIGDIR / "figure5_pathway_level_recovery.png", bbox_inches="tight")
    fig.savefig(FIGDIR / "figure5_pathway_level_recovery.pdf", bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(source_rows).to_csv(SRCDIR / "figure5_pathway_level_recovery.csv", index=False)

    # Supplementary overlap figure.
    ov = overlap.copy()
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    data = [ov.loc[ov["dataset"].eq(d), "n_genes_overlap"].to_numpy() for d in ov["dataset"].drop_duplicates()]
    ax.boxplot(data, labels=[d.replace("_24h_top2000", "").replace("_neurips2023", "") for d in ov["dataset"].drop_duplicates()], showfliers=False)
    ax.set_ylabel("Genes overlapping each Hallmark set")
    ax.set_title("Hallmark gene-set overlap")
    fig.tight_layout()
    fig.savefig(FIGDIR / "supp_pathway_gene_set_overlap.png", bbox_inches="tight")
    plt.close(fig)
    ov.to_csv(SRCDIR / "supp_pathway_gene_set_overlap.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    for split, g in seed_summary[seed_summary["dataset"].eq("sciplex3_24h_top2000")].groupby("split", observed=True):
        ax.scatter(
            np.repeat(split, len(g)),
            g["mean_pathway_pearson"],
            s=12,
            alpha=0.45,
            label=split.replace("_", " "),
        )
    ax.set_ylabel("Seed-level pathway Pearson")
    ax.set_title("Sci-Plex 3 pathway recovery by seed")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(FIGDIR / "supp_pathway_per_seed_distribution.png", bbox_inches="tight")
    plt.close(fig)
    seed_summary.to_csv(SRCDIR / "supp_pathway_per_seed_distribution.csv", index=False)

    if not heat.empty:
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="magma")
        ax.set_xticks(np.arange(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=35, ha="right")
        ax.set_yticks(np.arange(matrix.shape[0]))
        ax.set_yticklabels([s.replace("HALLMARK_", "").replace("_", " ") for s in matrix.index], fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Mean absolute error")
        ax.set_title("Representative Hallmark error heatmap")
        fig.tight_layout()
        fig.savefig(FIGDIR / "supp_pathway_per_gene_set_error_heatmap.png", bbox_inches="tight")
        plt.close(fig)
        heat.to_csv(SRCDIR / "supp_pathway_per_gene_set_error_heatmap.csv", index=False)


def write_check_report(gene_sets: list[GeneSet], overlap: pd.DataFrame, notes: list[str]) -> None:
    lines = [
        "# Pathway Gene Set Check",
        "",
        f"Local GMT file: `{GMT}`",
        f"Gene sets parsed: {len(gene_sets)}",
        f"Minimum overlap threshold: {MIN_OVERLAP} genes",
        "",
        "## Dataset Overlap",
        "",
    ]
    for dataset, g in overlap.groupby("dataset", observed=True):
        usable = int(g["status"].eq("usable").sum())
        lines.extend(
            [
                f"### {dataset}",
                f"- Gene sets with >= {MIN_OVERLAP} overlapping genes: {usable}/{len(g)}",
                f"- Median overlap: {float(g['n_genes_overlap'].median()):.1f}",
                f"- Minimum overlap: {int(g['n_genes_overlap'].min())}",
                f"- Maximum overlap: {int(g['n_genes_overlap'].max())}",
                "",
            ]
        )
    lines.extend(
        [
            "## Model Inclusion Notes",
            "",
            "- Included models are limited to runs with full predicted response vectors or reconstructable train-only baseline predictions.",
            "- Nearest-drug and SVD-ridge OpenProblems outputs were not projected because the available local files contain per-row metrics but not full prediction vectors.",
            "",
            "## Run Notes",
            "",
        ]
    )
    lines.extend(f"- {n}" for n in notes)
    (DOCSDIR / "pathway_gene_set_check.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ensure_dirs()
    if not GMT.exists():
        msg = f"Local GMT file not found: {GMT}"
        (DOCSDIR / "pathway_analysis_feasibility_report.md").write_text(msg + "\n", encoding="utf-8")
        (ROOT / "summary_pathway_failed.md").write_text(msg + "\n", encoding="utf-8")
        return 2

    gene_sets = parse_gmt(GMT)
    metrics_path = OUTDIR / "pathway_metrics_long.csv"
    if metrics_path.exists():
        metrics_path.unlink()

    all_overlap: list[dict] = []
    all_seed_summary: list[dict] = []
    all_gene_errors: list[dict] = []
    all_notes: list[str] = []

    overlap, seed_rows, gene_rows, notes = run_openproblems(gene_sets, metrics_path)
    all_overlap.extend(overlap)
    all_seed_summary.extend(seed_rows)
    all_gene_errors.extend(gene_rows)
    all_notes.extend(notes)

    overlap, seed_rows, gene_rows, notes = run_sciplex(gene_sets, metrics_path)
    all_overlap.extend(overlap)
    all_seed_summary.extend(seed_rows)
    all_gene_errors.extend(gene_rows)
    all_notes.extend(notes)

    overlap_df = pd.DataFrame(all_overlap)
    seed_summary = pd.DataFrame(all_seed_summary)
    per_gene = pd.DataFrame(all_gene_errors)
    if seed_summary.empty:
        msg = "No pathway-level runs completed; check gene-set overlap and prediction-vector availability."
        (DOCSDIR / "pathway_analysis_feasibility_report.md").write_text(msg + "\n", encoding="utf-8")
        (ROOT / "summary_pathway_failed.md").write_text(msg + "\n", encoding="utf-8")
        return 3

    contrasts, bootstrap = summarize_outputs(seed_summary, per_gene)
    summary_stats = add_summary_stats(seed_summary)

    overlap_df.to_csv(OUTDIR / "pathway_gene_set_overlap.csv", index=False)
    seed_summary.to_csv(OUTDIR / "pathway_seed_summary.csv", index=False)
    per_gene.to_csv(OUTDIR / "pathway_per_gene_set_errors.csv", index=False)
    contrasts.to_csv(OUTDIR / "pathway_random_to_strict_contrasts.csv", index=False)
    bootstrap.to_csv(OUTDIR / "pathway_bootstrap_ci.csv", index=False)
    summary_stats.to_csv(OUTDIR / "pathway_summary_stats.csv", index=False)

    write_figures(seed_summary, contrasts, per_gene, overlap_df)
    write_check_report(gene_sets, overlap_df, all_notes)

    print("Completed pathway-level recovery analysis.")
    print(summary_stats.to_string(index=False))
    print(contrasts.to_string(index=False))
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    raise SystemExit(main())
