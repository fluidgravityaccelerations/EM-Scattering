import numpy as np
from emscattering import current_circle, farfield_circle, farfield_pattern, rms_relative_error, solve_efie_nystrom


def test_current_circle_accuracy():
    result = solve_efie_nystrom(wavelength=1.0, a=2.0, b=2.0, n_intervals=60)
    ref = current_circle(result.t, 2.0, result.k, result.eta)
    assert rms_relative_error(result.current, ref) < 0.05


def test_farfield_circle_accuracy():
    result = solve_efie_nystrom(wavelength=1.0, a=2.0, b=2.0, n_intervals=60)
    phi = np.linspace(0.0, 2.0*np.pi, 180, endpoint=False)
    f = farfield_pattern(result, phi)
    ref = farfield_circle(phi, 2.0, result.k)
    f = f / np.max(np.abs(f))
    ref = ref / np.max(np.abs(ref))
    assert rms_relative_error(np.abs(f), np.abs(ref)) < 0.2
