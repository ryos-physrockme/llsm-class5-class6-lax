"""Reconstruct the documented Class-6 nilpotent-polynomial coefficient search.

The optimization uses the finite matrix library {I, N, N^2}, off-shell local
jets, and the F=A E certificate. It is a clean reimplementation from the
research note rather than a byte-for-byte historical training script.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ml"
OUT.mkdir(parents=True, exist_ok=True)

IDENTITY = np.eye(3, dtype=complex)
N = np.array([[0, 0, -1], [0, 0, -1j], [-1, -1j, 0]], dtype=complex)
N2 = N @ N

RUN_PARAMETERS = [
    (0.70, 0.15),
    (0.90, 0.25),
    (1.10, 0.35),
    (1.30, 0.45),
    (1.50, 0.55),
    (1.80, 0.70),
]


def local_jets(seed: int, size: int) -> list[tuple[np.ndarray, ...]]:
    rng = np.random.default_rng(seed)
    jets = []
    for _ in range(size):
        spin = rng.normal(size=3)
        spin /= np.linalg.norm(spin)
        spin_x = rng.normal(size=3)
        spin_x -= spin * np.dot(spin, spin_x)
        spin_xx = rng.normal(size=3)
        spin_xx += spin * (-np.dot(spin_x, spin_x) - np.dot(spin, spin_xx))
        spin_t = rng.normal(size=3)
        spin_t -= spin * np.dot(spin, spin_t)
        jets.append((spin, spin_x, spin_xx, spin_t))
    return jets


def matrices(parameters: np.ndarray, a0: float) -> tuple[np.ndarray, ...]:
    a1, a2, b0, b1, b2, c0, c1, c2 = parameters
    A = a0 * IDENTITY + a1 * N + a2 * N2
    B = b0 * IDENTITY + b1 * N + b2 * N2
    C = c0 * IDENTITY + c1 * N + c2 * N2
    return A, B, C


def residual_vector(
    parameters: np.ndarray,
    a0: float,
    alpha: float,
    jets: list[tuple[np.ndarray, ...]],
) -> np.ndarray:
    A, B, C = matrices(parameters, a0)
    residuals: list[float] = []
    for spin, spin_x, spin_xx, spin_t in jets:
        spatial = A @ spin
        spatial_t = A @ spin_t
        temporal = B @ spin + C @ np.cross(spin, spin_x)
        temporal_x = B @ spin_x + C @ np.cross(spin, spin_xx)
        curvature = spatial_t - temporal_x + np.cross(spatial, temporal)
        eom = spin_t - np.cross(spin, spin_xx - alpha * N @ spin)
        delta = curvature - A @ eom
        residuals.extend(delta.real)
        residuals.extend(delta.imag)
    return np.asarray(residuals, dtype=float)


def exact_coefficients(a0: float, alpha: float) -> np.ndarray:
    return np.array(
        [
            alpha / (2 * a0),
            -alpha**2 / (8 * a0**3),
            -a0**2,
            alpha / 2,
            -3 * alpha**2 / (8 * a0**2),
            a0,
            alpha / (2 * a0),
            -alpha**2 / (8 * a0**3),
        ],
        dtype=float,
    )


def run() -> dict:
    records = []
    for run_index, (a0, alpha) in enumerate(RUN_PARAMETERS):
        train = local_jets(40_000 + run_index, 64)
        validation = local_jets(50_000 + run_index, 64)
        rng = np.random.default_rng(60_000 + run_index)
        initial = rng.normal(scale=0.2, size=8)
        solution = least_squares(
            residual_vector,
            initial,
            args=(a0, alpha, train),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=3000,
        )
        exact = exact_coefficients(a0, alpha)
        validation_residual = residual_vector(solution.x, a0, alpha, validation)
        A, _, _ = matrices(solution.x, a0)
        square_target = a0**2 * IDENTITY + alpha * N
        records.append(
            {
                "run": run_index,
                "a0": a0,
                "alpha": alpha,
                "coefficients": solution.x.tolist(),
                "exact_coefficients": exact.tolist(),
                "max_coefficient_error": float(np.max(np.abs(solution.x - exact))),
                "validation_rms": float(np.sqrt(np.mean(validation_residual**2))),
                "square_root_residual": float(np.linalg.norm(A @ A - square_target)),
                "nfev": int(solution.nfev),
            }
        )

    summary = {
        "implementation": "clean reconstruction from documented equations",
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "number_of_runs": len(records),
        "max_coefficient_error": max(r["max_coefficient_error"] for r in records),
        "max_validation_rms": max(r["validation_rms"] for r in records),
        "max_square_root_residual": max(r["square_root_residual"] for r in records),
        "runs": records,
    }
    path = OUT / "class6_reproduced_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "runs"}, indent=2))
    if summary["max_coefficient_error"] > 1e-8:
        raise SystemExit("Class-6 coefficient recovery did not reach the expected branch.")
    return summary


if __name__ == "__main__":
    run()
