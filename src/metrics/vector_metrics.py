"""Shared vector metrics for perturbation-response predictions."""

from __future__ import annotations

import numpy as np


def rowwise_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: {y_true.shape} != {y_pred.shape}")
    true_centered = y_true - np.nanmean(y_true, axis=1, keepdims=True)
    pred_centered = y_pred - np.nanmean(y_pred, axis=1, keepdims=True)
    numerator = np.nansum(true_centered * pred_centered, axis=1)
    denominator = np.sqrt(np.nansum(true_centered**2, axis=1) * np.nansum(pred_centered**2, axis=1))
    out = numerator / denominator
    out[denominator == 0] = np.nan
    return out


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: {y_true.shape} != {y_pred.shape}")
    return np.sqrt(np.nanmean((y_true - y_pred) ** 2, axis=1))


def topk_overlap(y_true: np.ndarray, y_pred: np.ndarray, k: int) -> np.ndarray:
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: {y_true.shape} != {y_pred.shape}")
    k = min(k, y_true.shape[1])
    true_idx = np.argpartition(np.abs(y_true), -k, axis=1)[:, -k:]
    pred_idx = np.argpartition(np.abs(y_pred), -k, axis=1)[:, -k:]
    return np.array([len(set(t).intersection(set(p))) / k for t, p in zip(true_idx, pred_idx)])


def topk_direction_agreement(y_true: np.ndarray, y_pred: np.ndarray, k: int) -> np.ndarray:
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: {y_true.shape} != {y_pred.shape}")
    k = min(k, y_true.shape[1])
    pred_idx = np.argpartition(np.abs(y_pred), -k, axis=1)[:, -k:]
    rows = []
    for i, idx in enumerate(pred_idx):
        true_sign = np.sign(y_true[i, idx])
        pred_sign = np.sign(y_pred[i, idx])
        nonzero = (true_sign != 0) & (pred_sign != 0)
        rows.append(float(np.mean(true_sign[nonzero] == pred_sign[nonzero])) if np.any(nonzero) else np.nan)
    return np.asarray(rows)
