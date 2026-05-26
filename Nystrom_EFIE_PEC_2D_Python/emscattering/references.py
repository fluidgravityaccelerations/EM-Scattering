"""Analytical circular-cylinder reference solutions."""

from __future__ import annotations

import numpy as np
from scipy import special
from .constants import ETA0


def current_circle(phi: np.ndarray, radius: float, k: float, eta0: float = ETA0, e0: complex = 1.0) -> np.ndarray:
    """Exact TMz surface current on a circular PEC cylinder.

    Plane wave: E_inc = E0 * exp(-1j*k*x).
    Time convention follows the MATLAB example and outgoing waves use H_n^(2).
    """
    phi = np.asarray(phi, dtype=float)
    if radius <= 0 or k <= 0:
        raise ValueError("radius and k must be positive")

    nmax = int(np.ceil(k * radius + 10))
    out = np.zeros_like(phi, dtype=np.complex128)
    for n in range(-nmax, nmax + 1):
        out += (1j ** (-n)) * np.exp(1j * n * phi) / special.hankel2(n, k * radius)
    return (2.0 * e0 / (np.pi * radius * k * eta0)) * out


def farfield_circle(phi: np.ndarray, radius: float, k: float, e0: complex = 1.0) -> np.ndarray:
    """Reference asymptotic far-field pattern for a circular PEC cylinder."""
    phi = np.asarray(phi, dtype=float)
    if radius <= 0 or k <= 0:
        raise ValueError("radius and k must be positive")

    nmax = int(np.ceil(k * radius + 20))
    out = np.zeros_like(phi, dtype=np.complex128)
    for n in range(-nmax, nmax + 1):
        out += special.jv(n, k * radius) / special.hankel2(n, k * radius) * np.exp(1j * n * phi)
    return -e0 * np.sqrt(2j / (np.pi * k)) * out
