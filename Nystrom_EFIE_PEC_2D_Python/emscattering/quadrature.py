"""Colton-Kress trigonometric quadrature weights."""

from __future__ import annotations

import numpy as np


def colton_kress_weight(n_intervals: int, t: float, tj: float) -> float:
    """Scalar Colton-Kress weight for the periodic logarithmic kernel."""
    if n_intervals <= 0:
        raise ValueError("n_intervals must be positive")
    out = 0.0
    for m in range(1, n_intervals):
        out += np.cos(m * (t - tj)) / m
    return -(2.0 * np.pi / n_intervals) * out - (np.pi / n_intervals**2) * np.cos(
        n_intervals * (t - tj)
    )


def colton_kress_weights(n_intervals: int, t: np.ndarray) -> np.ndarray:
    """Dense matrix of Colton-Kress weights for all pairs of nodes in t."""
    if n_intervals <= 0:
        raise ValueError("n_intervals must be positive")
    t = np.asarray(t, dtype=float)
    diff = t[:, None] - t[None, :]
    weights = np.zeros_like(diff)
    for m in range(1, n_intervals):
        weights += np.cos(m * diff) / m
    weights = -(2.0 * np.pi / n_intervals) * weights
    weights -= (np.pi / n_intervals**2) * np.cos(n_intervals * diff)
    return weights
