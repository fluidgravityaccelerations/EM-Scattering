"""Small Python port of the 2-D PEC EFIE Nyström MATLAB example."""

from .constants import EPS0, MU0, C0, ETA0
from .geometry import ellipse_geometry
from .quadrature import colton_kress_weights
from .references import current_circle, farfield_circle
from .nystrom import solve_efie_nystrom, farfield_pattern, rms_relative_error

__all__ = [
    "EPS0", "MU0", "C0", "ETA0",
    "ellipse_geometry", "colton_kress_weights",
    "current_circle", "farfield_circle",
    "solve_efie_nystrom", "farfield_pattern", "rms_relative_error",
]
