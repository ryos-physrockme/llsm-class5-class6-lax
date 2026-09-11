"""Reconstruct the documented Class-5 symmetry-informed coefficient search.

This is a clean reimplementation from the equations preserved in the full
research note, not a byte-for-byte copy of the historical training script.
It samples off-shell local jets, optimizes the finite feature-library
coefficients, evaluates held-out validation residuals, and writes a JSON
summary to results/ml/class5_reproduced_summary.json.
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

M5 = np.array([-1j, -1.0, 0.0], dtype=complex)
K5 = np.outer(M5, M5)

RUN_PARAMETERS = [
    (0.30, 0.70),
    (0.50, 0.80),
    (0.70, 1.00),
    (0.40, 1.20),
    (0.80, 1.30),
    (0.60, 1.50),
    (0.90, 1.70),
    (0.55, 2.00),
    (0.75, 2.30),
]


def local_jets(seed: int, size: int) -> list[tuple[np.ndarray, ...]]:
    """Sample unit-spin local jets without imposing the equation of motion."""
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


def residual_vector(
    parameters: np.ndarray,
    c: float,
    spectral_parameter: float,
    jets: list[tuple[np.ndarray, ...]],
) -> np.ndarray:
    """Off-shell F-QE residual for the over-complete local feature library."""
    beta = parameters[0]
    coefficients = parameters[1:]
    residuals: list[float] = []

    for spin, spin_x, spin_xx, spin_t in jets:
        chi = M5 @ spin
        xi = M5 @ spin_x
        varrho = M5 @ np.cross(spin, spin_x)
        varrho_x = M5 @ np.cross(spin, spin_xx)
        xi_x = M5 @ spin_xx

        features = [
            np.cross(spin, spin_x),
            varrho * M5,
            spin,
            chi * M5,
            spin_x,
            chi * spin_x,
            xi * spin,
            np.cross(M5, spin_x),
        ]
        feature_x = [
            np.cross(spin, spin_xx),
            varrho_x * M5,
            spin_x,
            xi * M5,
            spin_xx,
            xi * spin_x + chi * spin_xx,
            xi_x * spin + xi * spin_x,
            np.cross(M5, spin_xx),
        ]

        spatial = spectral_parameter * spin + beta * chi * M5
        spatial_t = spectral_parameter * spin_t + beta * (M5 @ spin_t) * M5
        temporal = sum(ca * fa for ca, fa in zip(coefficients, features))
        temporal_x = sum(ca * fa for ca, fa in zip(coefficients, feature_x))
        curvature = spatial_t - temporal_x + np.cross(spatial, temporal)

        eom = spin_t - np.cross(spin, spin_xx - c * c * M5 * chi)
        predicted = spectral_parameter * eom + beta * M5 * (M5 @ eom)
        delta = curvature - predicted
        residuals.extend(delta.real)
        residuals.extend(delta.imag)

    return np.asarray(residuals, dtype=float)


def exact_coefficients(c: float, spectral_parameter: float) -> np.ndarray:
    return np.array(
        [
            c * c / (2 * spectral_parameter),
            spectral_parameter,
            c * c / (2 * spectral_parameter),
            -spectral_parameter**2,
            c * c / 2,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )


def run() -> dict:
    records = []
    for run_index, (c, spectral_parameter) in enumerate(RUN_PARAMETERS):
        train = local_jets(10_000 + run_index, 64)
        validation = local_jets(20_000 + run_index, 64)
        rng = np.random.default_rng(30_000 + run_index)
        initial = rng.normal(scale=0.2, size=9)
        solution = least_squares(
            residual_vector,
            initial,
            args=(c, spectral_parameter, train),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=3000,
        )
        exact = exact_coefficients(c, spectral_parameter)
        validation_residual = residual_vector(
            solution.x, c, spectral_parameter, validation
        )
        validation_mse = float(np.mean(validation_residual**2))
        records.append(
            {
                "run": run_index,
                "c": c,
                "spectral_parameter": spectral_parameter,
                "coefficients": solution.x.tolist(),
                "exact_coefficients": exact.tolist(),
                "max_coefficient_error": float(np.max(np.abs(solution.x - exact))),
                "validation_mse": validation_mse,
                "distractor_max_abs": float(np.max(np.abs(solution.x[5:]))),
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
        "max_validation_mse": max(r["validation_mse"] for r in records),
        "median_validation_mse": float(np.median([r["validation_mse"] for r in records])),
        "max_distractor_abs": max(r["distractor_max_abs"] for r in records),
        "runs": records,
    }
    path = OUT / "class5_reproduced_summary.json"
    path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "runs"}, indent=2))
    if summary["max_coefficient_error"] > 1e-8:
        raise SystemExit("Class-5 coefficient recovery did not reach the expected branch.")
    return summary


if __name__ == "__main__":
    run()
