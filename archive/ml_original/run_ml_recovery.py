#!/usr/bin/env python3
"""ML recovery of a Class-5 companion matrix from an overcomplete local ansatz.

The exact coefficients are used only after training for evaluation.  The
training objective is the off-shell certificate ||F-QE||^2 plus a uniform, very light
sparsity penalty on all coefficients in the overcomplete library.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, NamedTuple, Tuple

import numpy as np
import pandas as pd
import torch


torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)
M = torch.tensor([-1j, -1.0, 0.0], dtype=torch.complex128)
FEATURE_NAMES = [
    "S_cross_Sx", "r_m", "S", "q_m",
    "Sx", "q_Sx", "p_S", "m_cross_Sx",
    "q_m_cross_S", "r_S", "m", "q2_m",
]


class Prepared(NamedTuple):
    S: torch.Tensor
    St: torch.Tensor
    q: torch.Tensor
    qt: torch.Tensor
    features: torch.Tensor
    features_x: torch.Tensor
    E: torch.Tensor
    m_dot_E: torch.Tensor


def dot_m(x: torch.Tensor) -> torch.Tensor:
    return torch.einsum("i,bi->b", M, x)


def sample_jets(n: int, seed: int) -> Tuple[torch.Tensor, ...]:
    g = torch.Generator().manual_seed(seed)
    S = torch.randn((n, 3), generator=g)
    S = S / torch.linalg.vector_norm(S, dim=1, keepdim=True)
    Sx = torch.randn((n, 3), generator=g)
    Sx -= torch.sum(S * Sx, dim=1, keepdim=True) * S
    St = torch.randn((n, 3), generator=g)
    St -= torch.sum(S * St, dim=1, keepdim=True) * S
    Sxx_tan = torch.randn((n, 3), generator=g)
    Sxx_tan -= torch.sum(S * Sxx_tan, dim=1, keepdim=True) * S
    Sxx = Sxx_tan - torch.sum(Sx * Sx, dim=1, keepdim=True) * S
    return tuple(x.to(torch.complex128) for x in (S, Sx, Sxx, St))


def prepare(n: int, seed: int, c: float) -> Prepared:
    S, Sx, Sxx, St = sample_jets(n, seed)
    q, p, px, qt = dot_m(S), dot_m(Sx), dot_m(Sxx), dot_m(St)
    Sx_cross = torch.linalg.cross(S, Sx)
    Sxx_cross = torch.linalg.cross(S, Sxx)
    r, rx = dot_m(Sx_cross), dot_m(Sxx_cross)
    mxS = torch.linalg.cross(M.expand_as(S), S)
    mxSx = torch.linalg.cross(M.expand_as(S), Sx)
    mxSxx = torch.linalg.cross(M.expand_as(S), Sxx)

    f = torch.stack([
        Sx_cross,
        r[:, None] * M,
        S,
        q[:, None] * M,
        Sx,
        q[:, None] * Sx,
        p[:, None] * S,
        mxSx,
        q[:, None] * mxS,
        r[:, None] * S,
        M.expand_as(S),
        q[:, None] ** 2 * M,
    ], dim=1)
    fx = torch.stack([
        Sxx_cross,
        rx[:, None] * M,
        Sx,
        p[:, None] * M,
        Sxx,
        p[:, None] * Sx + q[:, None] * Sxx,
        px[:, None] * S + p[:, None] * Sx,
        mxSxx,
        p[:, None] * mxS + q[:, None] * mxSx,
        rx[:, None] * S + r[:, None] * Sx,
        torch.zeros_like(S),
        2.0 * q[:, None] * p[:, None] * M,
    ], dim=1)
    E = St - torch.linalg.cross(S, Sxx - c * c * q[:, None] * M)
    return Prepared(S, St, q, qt, f, fx, E, dot_m(E))


def residual(data: Prepared, z: float, beta: torch.Tensor, coeff: torch.Tensor) -> torch.Tensor:
    cc = coeff.to(torch.complex128)
    bc = beta.to(torch.complex128)
    V = torch.einsum("k,bki->bi", cc, data.features)
    Vx = torch.einsum("k,bki->bi", cc, data.features_x)
    U = z * data.S + bc * data.q[:, None] * M
    Ut = z * data.St + bc * data.qt[:, None] * M
    F = Ut - Vx + torch.linalg.cross(U, V)
    QE = z * data.E + bc * data.m_dot_E[:, None] * M
    return F - QE


def objective(data: Prepared, z: float, beta: torch.Tensor, coeff: torch.Tensor, l1: float) -> torch.Tensor:
    R = residual(data, z, beta, coeff)
    mse = torch.mean(torch.abs(R) ** 2)
    # Uniform sparsity prior: no feature is labelled as physical during training.
    penalty = torch.sum(torch.sqrt(coeff ** 2 + 1e-20)) + torch.sqrt(beta ** 2 + 1e-20)
    return mse + l1 * penalty


def train_one(c: float, z: float, seed: int, train_n: int, val_n: int,
              adam_steps: int, lbfgs_steps: int, l1: float) -> Tuple[Dict[str, float], pd.DataFrame]:
    train = prepare(train_n, 1000 + seed, c)
    val = prepare(val_n, 100000 + seed, c)
    gen = torch.Generator().manual_seed(seed)
    beta = torch.nn.Parameter(0.25 * torch.randn((), generator=gen))
    coeff = torch.nn.Parameter(0.25 * torch.randn((len(FEATURE_NAMES),), generator=gen))
    opt = torch.optim.Adam([beta, coeff], lr=3e-2)
    history: List[Dict[str, float]] = []
    t0 = time.time()

    for step in range(adam_steps):
        opt.zero_grad(set_to_none=True)
        value = objective(train, z, beta, coeff, l1)
        value.backward()
        torch.nn.utils.clip_grad_norm_([beta, coeff], 100.0)
        opt.step()
        if step % 10 == 0 or step == adam_steps - 1:
            with torch.no_grad():
                vmse = torch.mean(torch.abs(residual(val, z, beta, coeff)) ** 2)
            history.append({"stage": "adam", "step": step,
                            "train_objective": float(value.detach()),
                            "validation_mse": float(vmse), "beta": float(beta)})

    solver = torch.optim.LBFGS([beta, coeff], lr=0.8, max_iter=lbfgs_steps,
                               tolerance_grad=1e-12, tolerance_change=1e-14,
                               line_search_fn="strong_wolfe")
    calls = 0
    def closure() -> torch.Tensor:
        nonlocal calls
        solver.zero_grad(set_to_none=True)
        value = objective(train, z, beta, coeff, l1 * 0.05)
        value.backward()
        calls += 1
        return value
    solver.step(closure)

    with torch.no_grad():
        tr = residual(train, z, beta, coeff)
        vr = residual(val, z, beta, coeff)
        train_mse = float(torch.mean(torch.abs(tr) ** 2))
        val_mse = float(torch.mean(torch.abs(vr) ** 2))
        point_norm = torch.linalg.vector_norm(vr, dim=1)
        learned = coeff.detach().numpy()
        beta_val = float(beta)

    beta_exact = c * c / (2.0 * z)
    exact = np.zeros(len(FEATURE_NAMES))
    exact[:4] = [z, beta_exact, -z * z, c * c / 2.0]
    row: Dict[str, float] = {
        "c": c, "alpha": 2*c, "z": z, "seed": seed,
        "beta_learned": beta_val, "beta_exact": beta_exact,
        "beta_abs_error": abs(beta_val-beta_exact),
        "train_mse": train_mse, "validation_mse": val_mse,
        "validation_rms_vector": float(torch.sqrt(torch.mean(point_norm**2))),
        "validation_max_vector": float(torch.max(point_norm)),
        "max_physical_coeff_error": float(np.max(np.abs(learned[:4]-exact[:4]))),
        "max_distractor_abs": float(np.max(np.abs(learned[4:]))),
        "elapsed_seconds": time.time()-t0, "lbfgs_closure_calls": calls,
    }
    for name, x, y in zip(FEATURE_NAMES, learned, exact):
        row[f"coef_{name}"] = float(x)
        row[f"exact_{name}"] = float(y)
    history.append({"stage": "lbfgs_final", "step": adam_steps,
                    "train_objective": train_mse,
                    "validation_mse": val_mse, "beta": beta_val})
    h = pd.DataFrame(history)
    h["c"], h["z"], h["seed"] = c, z, seed
    return row, h


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--train-n", type=int, default=768)
    p.add_argument("--val-n", type=int, default=4096)
    p.add_argument("--adam-steps", type=int, default=600)
    p.add_argument("--lbfgs-steps", type=int, default=50)
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--l1", type=float, default=1e-9)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    cases = [(0.4, 0.7), (0.4, 1.2), (0.7, 1.2)]
    rows, histories = [], []
    for c, z in cases:
        for seed in range(args.seeds):
            row, hist = train_one(c, z, seed, args.train_n, args.val_n,
                                  args.adam_steps, args.lbfgs_steps, args.l1)
            rows.append(row); histories.append(hist)
            print(f"c={c:.2f} z={z:.2f} seed={seed} "
                  f"beta_err={row['beta_abs_error']:.2e} "
                  f"coef_err={row['max_physical_coeff_error']:.2e} "
                  f"val={row['validation_mse']:.2e} "
                  f"distr={row['max_distractor_abs']:.2e}", flush=True)

    df = pd.DataFrame(rows)
    hist = pd.concat(histories, ignore_index=True)
    df.to_csv(args.out / "ml_recovery_summary.csv", index=False)
    hist.to_csv(args.out / "ml_training_history.csv", index=False)
    aggregate = {
        "runs": len(df),
        "max_beta_abs_error": float(df.beta_abs_error.max()),
        "median_beta_abs_error": float(df.beta_abs_error.median()),
        "max_physical_coeff_error": float(df.max_physical_coeff_error.max()),
        "median_physical_coeff_error": float(df.max_physical_coeff_error.median()),
        "max_validation_mse": float(df.validation_mse.max()),
        "median_validation_mse": float(df.validation_mse.median()),
        "max_distractor_abs": float(df.max_distractor_abs.max()),
        "median_distractor_abs": float(df.max_distractor_abs.median()),
    }
    (args.out / "ml_aggregate.json").write_text(json.dumps(aggregate, indent=2), encoding="utf-8")
    print(json.dumps(aggregate, indent=2), flush=True)


if __name__ == "__main__":
    main()
