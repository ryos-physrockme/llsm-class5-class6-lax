#!/usr/bin/env python3
"""Gradient-based discovery of the Class-5 Lax coefficients.

The script follows the coefficient-learning philosophy of Integrability Ex Machina:
for each spectral/deformation point, simple ansatz coefficients are optimized from
zero-curvature samples. The analytic dependence is not supplied to the optimizer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch

from class5_lax import (
    RDTYPE,
    exact_coefficients,
    rms,
    sample_complex_jets,
    transformed_curvature,
    transformed_eom_rhs,
)


def train_one(
    a: float,
    k: float,
    seed: int,
    n_train: int,
    n_valid: int,
    adam_steps: int,
) -> dict[str, float | int]:
    torch.manual_seed(seed)
    s, sx, sxx = sample_complex_jets(n_train, seed=100_000 + seed)
    st = transformed_eom_rhs(s, sx, sxx, k)
    sv, sxv, sxxv = sample_complex_jets(n_valid, seed=200_000 + seed)
    stv = transformed_eom_rhs(sv, sxv, sxxv, k)

    # Deliberately broad initialization; no exact coefficient is supplied.
    params = torch.nn.Parameter(0.6 * torch.randn(5, dtype=RDTYPE))
    optimizer = torch.optim.Adam([params], lr=3.0e-2)
    best_loss = float("inf")
    best = None
    t0 = time.perf_counter()
    for step in range(adam_steps):
        optimizer.zero_grad(set_to_none=True)
        f = transformed_curvature(s, sx, sxx, st, a, params)
        loss = torch.mean(torch.abs(f) ** 2)
        loss.backward()
        optimizer.step()
        value = float(loss.detach())
        if value < best_loss:
            best_loss = value
            best = params.detach().clone()
        if step > 400 and value < 1.0e-22:
            break
    if best is not None:
        params.data.copy_(best)

    # Deterministic quasi-Newton refinement on the same training set.
    lbfgs = torch.optim.LBFGS(
        [params],
        lr=0.8,
        max_iter=120,
        tolerance_grad=1.0e-13,
        tolerance_change=1.0e-15,
        line_search_fn="strong_wolfe",
    )

    def closure() -> torch.Tensor:
        lbfgs.zero_grad(set_to_none=True)
        f = transformed_curvature(s, sx, sxx, st, a, params)
        loss = torch.mean(torch.abs(f) ** 2)
        loss.backward()
        return loss

    lbfgs.step(closure)
    elapsed = time.perf_counter() - t0

    with torch.no_grad():
        f_train = transformed_curvature(s, sx, sxx, st, a, params)
        f_valid = transformed_curvature(sv, sxv, sxxv, stv, a, params)
        learned = params.detach().cpu().numpy()
    exact = exact_coefficients(a, k)
    row: dict[str, float | int] = {
        "a": a,
        "k": k,
        "seed": seed,
        "steps_requested": adam_steps,
        "elapsed_s": elapsed,
        "train_rms": rms(f_train),
        "valid_rms": rms(f_valid),
        "max_coefficient_error": float(np.max(np.abs(learned - exact))),
    }
    names = ["b", "d", "e", "g", "h"]
    for name, value, truth in zip(names, learned, exact, strict=True):
        row[f"{name}_learned"] = float(value)
        row[f"{name}_exact"] = float(truth)
        row[f"{name}_error"] = float(value - truth)
    return row


def fit_symbolic_library(df: pd.DataFrame) -> pd.DataFrame:
    """Fit each learned coefficient to a small interpretable feature library."""
    a = df["a"].to_numpy(float)
    k = df["k"].to_numpy(float)
    features = {
        "1": np.ones_like(a),
        "a": a,
        "a^2": a**2,
        "1/a": 1.0 / a,
        "k": k,
        "k/a": k / a,
        "k*a": k * a,
        "k/a^2": k / a**2,
        "k^2": k**2,
    }
    x = np.column_stack(list(features.values()))
    names = list(features)
    rows: list[dict[str, float | str]] = []
    # A simple exhaustive one-feature and two-feature search is more transparent
    # here than an opaque symbolic-regression package.
    for target in ["b_learned", "d_learned", "e_learned", "g_learned", "h_learned"]:
        y = df[target].to_numpy(float)
        candidates: list[tuple[float, str, np.ndarray]] = []
        for j, name in enumerate(names):
            c, *_ = np.linalg.lstsq(x[:, [j]], y, rcond=None)
            pred = x[:, [j]] @ c
            candidates.append((float(np.sqrt(np.mean((pred - y) ** 2))), f"{c[0]:+.12g}*{name}", pred))
        for j in range(len(names)):
            for l in range(j + 1, len(names)):
                c, *_ = np.linalg.lstsq(x[:, [j, l]], y, rcond=None)
                pred = x[:, [j, l]] @ c
                formula = f"{c[0]:+.12g}*{names[j]} {c[1]:+.12g}*{names[l]}"
                candidates.append((float(np.sqrt(np.mean((pred - y) ** 2))), formula, pred))
        candidates.sort(key=lambda t: (t[0], len(t[1])))
        err, formula, _ = candidates[0]
        rows.append({"target": target, "best_formula": formula, "fit_rms": err})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n-train", type=int, default=1024)
    parser.add_argument("--n-valid", type=int, default=4096)
    parser.add_argument("--adam-steps", type=int, default=1800)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    spectral_values = [0.55, 0.8, 1.15, 1.6]
    alpha_values = [0.45, 0.8, 1.25]
    seeds = [0, 1]
    rows = []
    for alpha in alpha_values:
        k = alpha * alpha / 4.0
        for a in spectral_values:
            for seed in seeds:
                row = train_one(a, k, seed, args.n_train, args.n_valid, args.adam_steps)
                row["alpha"] = alpha
                rows.append(row)
                print(
                    f"alpha={alpha:.3f} a={a:.3f} seed={seed} "
                    f"valid={row['valid_rms']:.3e} maxerr={row['max_coefficient_error']:.3e}",
                    flush=True,
                )
    df = pd.DataFrame(rows)
    df.to_csv(args.out / "ml_coefficients.csv", index=False)
    # Average seeds before symbolic distillation.
    mean_df = df.groupby(["alpha", "a", "k"], as_index=False).mean(numeric_only=True)
    symbolic = fit_symbolic_library(mean_df)
    symbolic.to_csv(args.out / "symbolic_distillation.csv", index=False)
    summary = {
        "runs": len(df),
        "max_validation_rms": float(df["valid_rms"].max()),
        "median_validation_rms": float(df["valid_rms"].median()),
        "max_coefficient_error": float(df["max_coefficient_error"].max()),
        "median_coefficient_error": float(df["max_coefficient_error"].median()),
        "n_train": args.n_train,
        "n_valid": args.n_valid,
        "adam_steps": args.adam_steps,
    }
    (args.out / "ml_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(symbolic.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
