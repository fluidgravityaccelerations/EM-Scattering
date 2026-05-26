from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

# Allow running the example directly from a source checkout without installation.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from emscattering import (  # noqa: E402
    current_circle,
    farfield_circle,
    farfield_pattern,
    rms_relative_error,
    solve_efie_nystrom,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Nyström EFIE solver for a circular or elliptical PEC cylinder")
    parser.add_argument("--wavelength", type=float, default=1.0)
    parser.add_argument("--a", type=float, default=2.0, help="semi-axis along x")
    parser.add_argument("--b", type=float, default=2.0, help="semi-axis along y")
    parser.add_argument("--n", type=int, default=150, help="N; the solver uses 2N nodes")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    result = solve_efie_nystrom(
        wavelength=args.wavelength,
        a=args.a,
        b=args.b,
        n_intervals=args.n,
    )

    print(f"Solved with 2N = {2*args.n} boundary nodes")

    is_circle = np.isclose(args.a, args.b)
    if is_circle:
        j_ref = current_circle(result.t, args.a, result.k, result.eta)
        err_j = rms_relative_error(result.current, j_ref)
        print(f"RMS current error against circular-cylinder series: {err_j:.6g} %")
    else:
        j_ref = None
        print("Current reference skipped because a != b")

    nphi = 360
    phi = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    f = farfield_pattern(result, phi)

    if is_circle:
        f_ref = farfield_circle(phi, args.a, result.k)
        f_norm = f / np.max(np.abs(f))
        f_ref_norm = f_ref / np.max(np.abs(f_ref))
        err_f = rms_relative_error(np.abs(f_norm), np.abs(f_ref_norm))
        print(f"RMS normalized far-field magnitude error against series: {err_f:.6g} %")
    else:
        f_ref = None

    if args.no_plots:
        return

    plt.figure()
    plt.plot(result.t, np.abs(result.current), "k", linewidth=1.2, label="Nyström")
    if j_ref is not None:
        plt.plot(result.t, np.abs(j_ref), "r--", linewidth=1.2, label="Reference series")
    plt.xlabel("t (rad)")
    plt.ylabel("|J(t)|")
    plt.title("Induced surface current density - magnitude")
    plt.grid(True)
    plt.legend()

    plt.figure()
    plt.plot(result.t, np.angle(result.current), "k", linewidth=1.2, label="Nyström")
    if j_ref is not None:
        plt.plot(result.t, np.angle(j_ref), "r--", linewidth=1.2, label="Reference series")
    plt.xlabel("t (rad)")
    plt.ylabel("Phase[J(t)] (rad)")
    plt.title("Induced surface current density - phase")
    plt.grid(True)
    plt.legend()

    plt.figure()
    if f_ref is not None:
        plt.plot(phi, np.abs(f) / np.max(np.abs(f)), "k", linewidth=1.2, label="Nyström")
        plt.plot(phi, np.abs(f_ref) / np.max(np.abs(f_ref)), "r--", linewidth=1.2, label="Reference series")
        plt.ylabel("Normalized |F(phi)|")
    else:
        plt.plot(phi, np.abs(f), "k", linewidth=1.2, label="Nyström")
        plt.ylabel("|F(phi)|")
    plt.xlabel("phi (rad)")
    plt.title("Far-field pattern")
    plt.grid(True)
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
