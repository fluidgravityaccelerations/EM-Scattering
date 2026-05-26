"""Nyström EFIE solver for 2-D TM scattering by a smooth PEC cylinder."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import special

from .constants import EPS0, MU0, C0, ETA0, EULER_GAMMA
from .geometry import BoundaryGeometry, ellipse_geometry
from .quadrature import colton_kress_weights


@dataclass(frozen=True)
class NystromResult:
    """Output of the EFIE Nyström solve."""

    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    dsdt: np.ndarray
    current: np.ndarray
    matrix: np.ndarray
    rhs: np.ndarray
    k: float
    omega: float
    mu: float
    epsilon: float
    eta: float


def assemble_single_layer_operator(
    geometry: BoundaryGeometry,
    k: float,
    n_intervals: int,
    *,
    matlab_convention: bool = True,
) -> np.ndarray:
    """Assemble the scaled EFIE matrix used in the MATLAB example.

    The MATLAB code writes M2 with H^(1), then takes a complex conjugate of the
    full quadrature matrix. For real k and real boundary points this is equivalent
    to an H^(2) outgoing-wave convention. The default keeps the MATLAB convention
    exactly, which makes the numerical output directly comparable.
    """
    if k <= 0:
        raise ValueError("k must be positive")

    t = geometry.t
    x = geometry.x
    y = geometry.y
    ds = geometry.dsdt
    m = t.size

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    distance = np.sqrt(dx**2 + dy**2)
    diff = t[:, None] - t[None, :]

    r_weights = colton_kress_weights(n_intervals, t)
    source_ds = ds[None, :]

    m1 = -(1.0 / (2.0 * np.pi)) * special.jv(0, k * distance) * source_ds
    m2 = np.empty((m, m), dtype=np.complex128)
    offdiag = ~np.eye(m, dtype=bool)

    logterm = np.empty((m, m), dtype=float)
    logterm[offdiag] = np.log(4.0 * np.sin(diff[offdiag] / 2.0) ** 2)

    if matlab_convention:
        m2[offdiag] = (
            0.5j * special.hankel1(0, k * distance[offdiag]) * np.broadcast_to(source_ds, (m, m))[offdiag]
            - m1[offdiag] * logterm[offdiag]
        )
        m2[np.eye(m, dtype=bool)] = (
            0.5j - EULER_GAMMA / np.pi - (1.0 / np.pi) * np.log((k / 2.0) * ds)
        ) * ds
        q = r_weights * m1 + (np.pi / n_intervals) * m2
        return np.conj(q)

    # Direct H^(2) form, algebraically equivalent to the MATLAB convention.
    m2[offdiag] = (
        -0.5j * special.hankel2(0, k * distance[offdiag]) * np.broadcast_to(source_ds, (m, m))[offdiag]
        - m1[offdiag] * logterm[offdiag]
    )
    m2[np.eye(m, dtype=bool)] = (
        -0.5j - EULER_GAMMA / np.pi - (1.0 / np.pi) * np.log((k / 2.0) * ds)
    ) * ds
    return r_weights * m1 + (np.pi / n_intervals) * m2


def solve_efie_nystrom(
    *,
    wavelength: float = 1.0,
    a: float = 2.0,
    b: float = 2.0,
    n_intervals: int = 150,
    e0: complex = 1.0,
    epsilon_r: float = 1.0,
    mu_r: float = 1.0,
    propagation_angle: float = 0.0,
    matlab_convention: bool = True,
) -> NystromResult:
    """Solve the TM EFIE for a PEC elliptical cylinder.

    Parameters
    ----------
    wavelength:
        Host-medium wavelength.
    a, b:
        Ellipse semi-axes. Use a=b for a circular cylinder.
    n_intervals:
        N in the 2N-node midpoint discretization.
    propagation_angle:
        Incident plane-wave direction angle. 0 means propagation along +x.
    matlab_convention:
        Keep the exact MATLAB assembly convention when True.
    """
    if wavelength <= 0:
        raise ValueError("wavelength must be positive")
    if epsilon_r <= 0 or mu_r <= 0:
        raise ValueError("relative material parameters must be positive")

    epsilon = EPS0 * epsilon_r
    mu = MU0 * mu_r
    c = 1.0 / np.sqrt(epsilon * mu)
    k = 2.0 * np.pi / wavelength
    omega = k * c
    eta = np.sqrt(mu / epsilon)

    geom = ellipse_geometry(a, b, n_intervals)
    q = assemble_single_layer_operator(geom, k, n_intervals, matlab_convention=matlab_convention)
    matrix = 0.5j * omega * mu * q

    khat_x = np.cos(propagation_angle)
    khat_y = np.sin(propagation_angle)
    rhs = e0 * np.exp(-1j * k * (geom.x * khat_x + geom.y * khat_y))

    current = np.linalg.solve(matrix, rhs)
    return NystromResult(
        t=geom.t,
        x=geom.x,
        y=geom.y,
        dsdt=geom.dsdt,
        current=current,
        matrix=matrix,
        rhs=rhs,
        k=k,
        omega=omega,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
    )


def farfield_pattern(result: NystromResult, phi: np.ndarray) -> np.ndarray:
    """Compute the asymptotic angular far-field pattern from a Nyström current."""
    phi = np.asarray(phi, dtype=float)
    phase = result.k * (
        np.outer(np.cos(phi), result.x) + np.outer(np.sin(phi), result.y)
    )
    return 0.25j * result.omega * result.mu * np.sum(
        np.exp(1j * phase) * result.current[None, :] * result.dsdt[None, :], axis=1
    )


def rms_relative_error(numerical: np.ndarray, reference: np.ndarray) -> float:
    """Return 100*||numerical-reference||_2/||reference||_2."""
    numerical = np.asarray(numerical)
    reference = np.asarray(reference)
    denom = np.sum(np.abs(reference) ** 2)
    if denom == 0:
        raise ValueError("reference norm is zero")
    return float(100.0 * np.sqrt(np.sum(np.abs(numerical - reference) ** 2) / denom))
