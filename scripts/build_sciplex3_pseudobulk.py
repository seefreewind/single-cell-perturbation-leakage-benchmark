#!/usr/bin/env python3
"""Build a lightweight Sci-Plex 3 pseudobulk response matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def compute_gene_sums(adata: ad.AnnData, row_idx: np.ndarray, chunk_size: int) -> np.ndarray:
    sums = np.zeros(adata.n_vars, dtype=np.float64)
    for start in range(0, len(row_idx), chunk_size):
        rows = row_idx[start : start + chunk_size]
        block = adata.X[rows, :]
        sums += np.asarray(block.sum(axis=0)).ravel()
        print(f"gene sums: {min(start + chunk_size, len(row_idx))}/{len(row_idx)} cells")
    return sums


def mean_log_norm(
    adata: ad.AnnData,
    row_idx: np.ndarray,
    gene_idx: np.ndarray,
    ncounts: np.ndarray,
) -> np.ndarray:
    # CSRDataset direct row+column fancy indexing materializes the full matrix.
    # Row slicing first keeps the read bounded to the cells in this group.
    x = adata.X[row_idx, :][:, gene_idx].toarray().astype(np.float32)
    denom = ncounts[row_idx].astype(np.float32)
    denom[denom <= 0] = 1.0
    x = np.log1p(x / denom[:, None] * 10000.0)
    return x.mean(axis=0, dtype=np.float64).astype(np.float32)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=Path("data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad"))
    parser.add_argument("--smiles", type=Path, default=Path("metadata/sciplex3_pubchem_smiles_with_scaffolds.csv"))
    parser.add_argument("--out", type=Path, default=Path("data/processed/sciplex3/sciplex3_24h_top5000_response.h5ad"))
    parser.add_argument("--metadata-out", type=Path, default=Path("metadata/sciplex3_24h_top5000_response_metadata.csv"))
    parser.add_argument("--time", type=float, default=24.0)
    parser.add_argument("--min-cells", type=int, default=100)
    parser.add_argument("--top-genes", type=int, default=5000)
    parser.add_argument("--chunk-size", type=int, default=5000)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy().reset_index(names="cell_id")
    obs["_row"] = np.arange(len(obs))
    smiles = pd.read_csv(args.smiles)

    keep_base = (
        obs["time"].eq(args.time)
        & obs["cell_line"].isin(["A549", "K562", "MCF7"])
        & obs["ncounts"].notna()
    )
    base_rows = obs.loc[keep_base, "_row"].to_numpy(dtype=int)
    gene_sums = compute_gene_sums(adata, base_rows, args.chunk_size)
    selected = np.argsort(gene_sums)[::-1][: args.top_genes]
    selected = np.sort(selected)
    var = adata.var.iloc[selected].copy()
    var["raw_count_sum_for_selection"] = gene_sums[selected]

    ncounts = obs["ncounts"].to_numpy(dtype=np.float64)
    control = obs.loc[
        keep_base & obs["perturbation"].astype(str).str.lower().eq("control")
    ].copy()
    control_means: dict[str, np.ndarray] = {}
    for cell_line, frame in control.groupby("cell_line", observed=True):
        rows = frame["_row"].to_numpy(dtype=int)
        if len(rows) < args.min_cells:
            continue
        control_means[str(cell_line)] = mean_log_norm(adata, rows, selected, ncounts)
        print(f"control {cell_line}: {len(rows)} cells")

    treated = obs.loc[
        keep_base
        & obs["perturbation"].notna()
        & ~obs["perturbation"].astype(str).str.lower().eq("control")
    ].copy()
    group_cols = ["perturbation", "cell_line", "dose_value", "time"]
    rows = []
    response_blocks = []
    treated_blocks = []
    control_blocks = []
    for key, frame in treated.groupby(group_cols, observed=True):
        perturbation, cell_line, dose_value, time_value = key
        if str(cell_line) not in control_means:
            continue
        if len(frame) < args.min_cells:
            continue
        idx = frame["_row"].to_numpy(dtype=int)
        treated_mean = mean_log_norm(adata, idx, selected, ncounts)
        control_mean = control_means[str(cell_line)]
        response = treated_mean - control_mean
        rows.append(
            {
                "perturbation": perturbation,
                "cell_line": cell_line,
                "dose_value": dose_value,
                "time": time_value,
                "n_cells": len(frame),
            }
        )
        response_blocks.append(response)
        treated_blocks.append(treated_mean)
        control_blocks.append(control_mean)
        if len(rows) % 100 == 0:
            print(f"aggregated {len(rows)} treated groups")

    out_obs = pd.DataFrame(rows)
    out_obs = out_obs.merge(
        smiles[
            [
                "query_name",
                "canonical_smiles",
                "isomeric_smiles",
                "pubchem_cid",
                "scaffold",
            ]
        ],
        left_on="perturbation",
        right_on="query_name",
        how="left",
    ).drop(columns=["query_name"])
    out_obs.index = [
        f"{r.perturbation}__{r.cell_line}__{r.dose_value:g}nM__{r.time:g}h"
        for r in out_obs.itertuples()
    ]
    x = np.vstack(response_blocks).astype(np.float32)
    result = ad.AnnData(X=x, obs=out_obs, var=var)
    result.layers["treated_mean_log1p_cp10k"] = np.vstack(treated_blocks).astype(np.float32)
    result.layers["control_mean_log1p_cp10k"] = np.vstack(control_blocks).astype(np.float32)
    result.uns["response_definition"] = "treated mean log1p(CP10K) minus matched cell-line control mean log1p(CP10K)"
    result.uns["source_h5ad"] = str(args.h5ad)
    result.uns["time_filter"] = args.time
    result.uns["min_cells"] = args.min_cells
    result.uns["top_genes"] = args.top_genes
    result.write_h5ad(args.out)
    out_obs.to_csv(args.metadata_out, index=True, index_label="sample_id")
    adata.file.close()
    print(f"wrote {args.out} with shape {result.shape}")


if __name__ == "__main__":
    main()
