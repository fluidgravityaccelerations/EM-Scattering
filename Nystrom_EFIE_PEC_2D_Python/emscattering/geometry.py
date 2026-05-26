"""Boundary parametrizations for 2-D scatterers."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BoundaryGeometry:
    """Samples of a smooth 2*pi-periodic boundary."""

    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    dxdt: np.ndarray
    dydt: np.ndarray
    dsdt: np.ndarray


def midpoint_nodes(n_intervals: int) -> np.ndarray:
    """Return the 2N midpoint nodes used by the Colton-Kress rule."""
    if n_intervals <= 0:
        raise ValueError("n_intervals must be positive")
    m = 2 * n_intervals
    return np.pi * (2 * np.arange(m) + 1) / (2 * n_intervals)


def ellipse_geometry(a: float, b: float, n_intervals: int) -> BoundaryGeometry:
    """Sample the ellipse x=a*cos(t), y=b*sin(t) at 2N midpoint nodes."""
    if a <= 0 or b <= 0:
        raise ValueError("ellipse semi-axes a and b must be positive")

    t = midpoint_nodes(n_intervals)
    x = a * np.cos(t)
    y = b * np.sin(t)
    dxdt = -a * np.sin(t)
    dydt = b * np.cos(t)
    dsdt = np.sqrt(dxdt**2 + dydt**2)
    return BoundaryGeometry(t=t, x=x, y=y, dxdt=dxdt, dydt=dydt, dsdt=dsdt)
