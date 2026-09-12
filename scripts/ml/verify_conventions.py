"""Check the restored searches against original code and manuscript matrices."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp
import torch

import class5_ml_recovery as current5
import paper_conventions as pc

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/ml/convention_verification.json"
checks = {}


def original(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "archive/ml_original" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def close(name, actual, expected, tolerance=2e-11):
    error = float(torch.max(torch.abs(actual - expected)).detach())
    assert error < tolerance, (name, error)
    checks[name] = error


def check_original_losses():
    old5 = original("original_class5_recovery", "run_ml_recovery.py")
    old6 = original("original_class56_lax", "class5_lax.py")
    assert len(old5.FEATURE_NAMES) == len(current5.FEATURE_NAMES) == 12
    c, z = 0.4, 0.7
    alpha5, lambda_ = 2 * c, 1j / (2 * z)
    old_data = old5.prepare(32, 1701, c)
    data = current5.prepare(32, 1701, alpha5)
    close("class5_sampler_spin", data.S, pc.apply_matrix(pc.ROTATION5, old_data.S))
    close("class5_sampler_time", data.St, -2 * pc.apply_matrix(pc.ROTATION5, old_data.St))
    generator = torch.Generator().manual_seed(984)
    beta = torch.randn((), generator=generator, requires_grad=True)
    coeff = torch.randn(12, generator=generator, requires_grad=True)
    old_r = old5.residual(old_data, z, beta, coeff)
    new_r = current5.residual(data, lambda_, beta, coeff)
    close("class5_residual_for_arbitrary_coefficients", new_r,
          pc.apply_matrix(pc.ROTATION5, old_r))
    loss_old = old5.objective(old_data, z, beta, coeff, 1e-9)
    loss_new = current5.objective(data, lambda_, beta, coeff, 1e-9)
    close("class5_loss", loss_new, loss_old)
    gradients_old = torch.autograd.grad(loss_old, (beta, coeff))
    gradients_new = torch.autograd.grad(loss_new, (beta, coeff))
    for i, (left, right) in enumerate(zip(gradients_new, gradients_old)):
        close(f"class5_loss_gradient_{i}", left, right)

    spin, sx, sxx = old6.sample_complex_jets(32, 234)
    spin_new, sx_new, sxx_new = pc.sample_complex_jets(32, 234)
    for name, left, right in zip(("spin", "spin_x", "spin_xx"),
                                 (spin_new, sx_new, sxx_new), (spin, sx, sxx)):
        close("class6_sampler_" + name, left, pc.apply_matrix(pc.ROTATION6, right))
    alpha6, leading = 0.35, 0.65
    lambda_ = 1j / (2 * leading)
    st = old6.class6_eom_rhs(spin, sx, sxx, alpha6)
    st_new = pc.class6_eom_rhs(spin_new, sx_new, sxx_new, alpha6)
    close("class6_paper_eom", st_new, -2 * pc.apply_matrix(pc.ROTATION6, st))
    coeff = torch.randn(8, generator=generator, requires_grad=True)
    old_f = old6.class6_curvature(spin, sx, sxx, st, leading, coeff)
    new_f = pc.class6_training_curvature(spin_new, sx_new, sxx_new, st_new, lambda_, coeff)
    close("class6_training_curvature_for_arbitrary_coefficients", new_f,
          pc.apply_matrix(pc.ROTATION6, old_f), 2e-10)
    old_loss, new_loss = torch.mean(abs(old_f)**2), torch.mean(abs(new_f)**2)
    close("class6_loss", new_loss, old_loss, 2e-10)
    close("class6_loss_gradient", torch.autograd.grad(new_loss, coeff)[0],
          torch.autograd.grad(old_loss, coeff)[0], 2e-10)

    # Independent time data exercise the off-shell extension against the
    # archived on-shell implementation for arbitrary, untrained coefficients.
    st_off = pc.sample_class6_time_derivative(spin_new, 951)
    close("class6_independent_time_is_tangent", pc.bilinear_dot(spin_new, st_off),
          torch.zeros(spin_new.shape[0], dtype=pc.CDTYPE))
    assert torch.max(abs(st_off-st_new)) > 1
    off_r = pc.class6_training_residual(spin_new, sx_new, sxx_new, st_off,
                                       lambda_, coeff, alpha6)
    old_f = old6.class6_curvature(spin, sx, sxx, st, leading, coeff)
    close("class6_offshell_matches_original_onshell", off_r,
          pc.apply_matrix(pc.ROTATION6, old_f), 2e-10)
    off_r2 = pc.class6_training_residual(spin_new, sx_new, sxx_new, 3*st_off,
                                        lambda_, coeff, alpha6)
    close("class6_offshell_time_independence", off_r2, off_r, 2e-10)
    off_loss, old_loss = torch.mean(abs(off_r)**2), torch.mean(abs(old_f)**2)
    close("class6_offshell_loss", off_loss, old_loss, 2e-10)
    close("class6_offshell_loss_gradient", torch.autograd.grad(off_loss, coeff)[0],
          torch.autograd.grad(old_loss, coeff)[0], 2e-10)


def check_exact_matrices():
    ii = sp.I
    alpha, lam = sp.symbols("alpha lambda", nonzero=True)
    p, q, z, px, qx, zx = sp.symbols("p q z px qx zx")
    spin = sp.Matrix([(p + q)/2, (p - q)/(2*ii), z])
    sx = sp.Matrix([(px + qx)/2, (px - qx)/(2*ii), zx])
    sigma = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -ii], [ii, 0]]), sp.diag(1, -1)]
    iota = lambda v: -ii/2 * sum((v[i] * sigma[i] for i in range(3)), sp.zeros(2))
    constraints = sp.groebner([p*q + z*z - 1, q*px + p*qx + 2*z*zx],
                              px, qx, zx, p, q, z, domain=sp.QQ_I.frac_field(alpha, lam))

    def zero(name, matrix):
        remainders = [sp.factor(constraints.reduce(sp.expand(x))[1]) for x in matrix]
        assert all(x == 0 for x in remainders), (name, remainders)
        checks[name] = "exact_zero"

    m = sp.Matrix([ii, -1, 0])
    matrix5 = sp.Matrix([[0, 0, -1], [0, 0, -ii], [1, ii, 0]])
    rotation = sp.eye(3) + alpha*lam*matrix5 + (alpha*lam)**2*matrix5**2/2
    leading = ii/(2*lam)
    beta = alpha**2/(8*leading)
    spatial5 = leading*sp.eye(3) + beta*(m*m.T)
    ds = sx + alpha/2*matrix5*spin
    u5 = rotation*spatial5*spin - alpha/2*m
    v5 = -2*rotation*(spatial5*spin.cross(ds) - leading**2*spin
                      + alpha**2/8*m*(m.dot(spin)))
    mu5 = 1/(4*lam) + alpha*p/4
    nu5 = ii*alpha/2*(p*zx-z*px) + ii*alpha**2*p*p/4
    U5 = iota(u5) + mu5*sp.eye(2)
    V5 = iota(v5) + nu5*sp.eye(2)
    paper_u5 = sp.Matrix([[1+z+2*alpha*lam*p, q-2*alpha*lam*(1+z)], [p, 1-z]])/(4*lam)
    paper_v5 = ii/(4*lam**2)*sp.Matrix([
        [z+alpha*lam*p*(1+z)+2*alpha**2*lam**2*p*p+lam*(q*px-p*qx)+4*alpha*lam**2*(p*zx-z*px),
         q-alpha*lam*(1+z)**2-2*alpha**2*lam**2*p*(1+z)+2*lam*(z*qx-q*zx)+2*alpha*lam**2*(p*qx-q*px)],
        [p+alpha*lam*p*p+2*lam*(p*zx-z*px), -z-alpha*lam*p*(1+z)+lam*(p*qx-q*px)]])
    zero("class5_exact_spatial_matrix_matches_paper", U5-paper_u5)
    zero("class5_exact_time_matrix_matches_paper", V5-paper_v5)
    zero("class5_rotation_orthogonality", rotation.T*rotation-sp.eye(3))
    assert sp.factor(rotation.det()) == 1
    checks["class5_rotation_determinant"] = 1

    n6 = sp.Matrix([[0, 0, -1], [0, 0, -ii], [-1, -ii, 0]])
    spatial6 = leading*sp.eye(3) + alpha/(2*leading)*n6 - alpha**2/(8*leading**3)*n6**2
    temporal6 = -leading**2*sp.eye(3) + alpha/2*n6 - 3*alpha**2/(8*leading**2)*n6**2
    U6 = iota(spatial6*spin) + sp.eye(2)/(4*lam)
    V6 = iota(-2*(temporal6*spin + spatial6*spin.cross(sx)))
    paper_u6 = sp.Matrix([[1+z+2*alpha*lam**2*p, q+4*alpha*lam**2*z-4*alpha**2*lam**4*p],
                         [p, 1-z-2*alpha*lam**2*p]])/(4*lam)
    v11 = z-2*alpha*lam**2*p+lam*(q*px-p*qx)+4*alpha*lam**3*(p*zx-z*px)
    v12 = (q-4*alpha*lam**2*z+12*alpha**2*lam**4*p+2*lam*(z*qx-q*zx)
           +4*alpha*lam**3*(q*px-p*qx)+8*alpha**2*lam**5*(z*px-p*zx))
    paper_v6 = ii/(4*lam**2)*sp.Matrix([[v11, v12], [p+2*lam*(p*zx-z*px), -v11]])
    zero("class6_exact_spatial_matrix_matches_paper", U6-paper_u6)
    zero("class6_exact_time_matrix_matches_paper", V6-paper_v6)
    zero("class6_square_root", spatial6**2-leading**2*sp.eye(3)-alpha*n6)


if __name__ == "__main__":
    check_original_losses()
    check_exact_matrices()
    OUT.write_text(json.dumps({"checks": checks, "all_passed": True}, indent=2) + "\n")
    print(json.dumps({"checks": len(checks), "all_passed": True}, indent=2))
