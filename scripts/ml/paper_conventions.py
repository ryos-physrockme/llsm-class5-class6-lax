"""Original finite-coefficient Lax searches in the manuscript conventions.

The spin is S, time is t, continuum couplings are alpha5/alpha6, and
lambda_ is the spectral parameter of the lattice-derived manuscript matrices.
The historical real coefficient domain is retained by sampling lambda_ on
the positive imaginary axis, so i/(2 lambda_) is positive and real.
"""
from __future__ import annotations

import numpy as np
import torch

CDTYPE = torch.complex128
RDTYPE = torch.float64
M5 = torch.tensor([[0, 0, -1], [0, 0, -1j], [1, 1j, 0]], dtype=CDTYPE)
m5 = torch.tensor([1j, -1, 0], dtype=CDTYPE)
N6 = torch.tensor([[0, 0, -1], [0, 0, -1j], [-1, -1j, 0]], dtype=CDTYPE)
ROTATION5 = torch.diag(torch.tensor([-1, 1, -1], dtype=CDTYPE))
ROTATION6 = torch.diag(torch.tensor([1, -1, -1], dtype=CDTYPE))
IDENTITY3 = torch.eye(3, dtype=CDTYPE)
PAULI = torch.tensor([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                      [[1, 0], [0, -1]]], dtype=CDTYPE)


def apply_matrix(matrix, vector):
    return torch.einsum("ij,...j->...i", matrix, vector)


def spinor(vector):
    """iota(v)=-i v^i sigma^i/2; commutators represent cross products."""
    return -0.5j * torch.einsum("...i,ijk->...jk", vector, PAULI)


def spectral_coefficient(lambda_):
    value = 1j / (2 * complex(lambda_))
    if abs(value.imag) > 1e-13 or value.real <= 0:
        raise ValueError("The original real-parameter runs use positive imaginary lambda.")
    return value.real


def class5_rotation(alpha5, lambda_, inverse=False):
    coefficient = alpha5 * complex(lambda_) * (-1 if inverse else 1)
    return IDENTITY3 + coefficient * M5 + coefficient**2 * (M5 @ M5) / 2


def class5_scalar_parts(spin, spin_x, alpha5, lambda_):
    plus = spin[..., 0] + 1j * spin[..., 1]
    plus_x = spin_x[..., 0] + 1j * spin_x[..., 1]
    spatial = 1 / (4 * lambda_) + alpha5 * plus / 4
    temporal = (0.5j * alpha5 * (plus * spin_x[..., 2] - spin[..., 2] * plus_x)
                + 0.25j * alpha5**2 * plus**2)
    return spatial, temporal


def class6_matrix_torch(device="cpu"):
    return N6.to(device)


def bilinear_dot(left, right):
    return torch.sum(left * right, dim=-1)


def sample_complex_jets(n_samples, seed, device="cpu", scale_x=0.8, scale_xx=0.8):
    """Original Class 6 sampler, followed by its constant spin-basis rotation."""
    gen = torch.Generator(device="cpu").manual_seed(seed)
    out_s, out_x, out_xx = [], [], []
    remaining = n_samples
    while remaining > 0:
        batch = max(2 * remaining, 256)
        z = torch.randn(batch, 3, generator=gen, dtype=RDTYPE)
        z = z + 1j * torch.randn(batch, 3, generator=gen, dtype=RDTYPE)
        q = bilinear_dot(z, z)
        mask = torch.abs(q) > 0.25
        z, q = z[mask], q[mask]
        if z.shape[0] == 0:
            continue
        spin = z / torch.sqrt(q)[:, None]
        r1 = torch.randn(spin.shape[0], 3, generator=gen, dtype=RDTYPE)
        r1 = r1 + 1j * torch.randn(spin.shape[0], 3, generator=gen, dtype=RDTYPE)
        spin_x = scale_x * (r1 - spin * bilinear_dot(spin, r1)[:, None])
        r2 = torch.randn(spin.shape[0], 3, generator=gen, dtype=RDTYPE)
        r2 = r2 + 1j * torch.randn(spin.shape[0], 3, generator=gen, dtype=RDTYPE)
        tangent = r2 - spin * bilinear_dot(spin, r2)[:, None]
        spin_xx = scale_xx * tangent - spin * bilinear_dot(spin_x, spin_x)[:, None]
        take = min(remaining, spin.shape[0])
        for output, value in zip((out_s, out_x, out_xx), (spin, spin_x, spin_xx)):
            output.append(value[:take])
        remaining -= take
    return tuple(apply_matrix(ROTATION6, torch.cat(values)).to(device=device)
                 for values in (out_s, out_x, out_xx))


def class6_eom_rhs(spin, spin_x, spin_xx, alpha6):
    return -2 * torch.linalg.cross(spin, spin_xx - alpha6 * apply_matrix(N6, spin))


def class6_coefficient_matrices(lambda_, coefficients):
    a1, a2, b0, b1, b2, c0, c1, c2 = coefficients.to(CDTYPE)
    n2 = N6 @ N6
    spatial = spectral_coefficient(lambda_) * IDENTITY3 + a1 * N6 + a2 * n2
    temporal = b0 * IDENTITY3 + b1 * N6 + b2 * n2
    derivative = c0 * IDENTITY3 + c1 * N6 + c2 * n2
    return spatial, temporal, derivative


def class6_vectors(spin, spin_x, lambda_, coefficients):
    spatial, temporal, derivative = class6_coefficient_matrices(lambda_, coefficients)
    return (apply_matrix(spatial, spin),
            -2 * (apply_matrix(temporal, spin)
                  + apply_matrix(derivative, torch.linalg.cross(spin, spin_x))))


def class6_curvature(spin, spin_x, spin_xx, spin_t, lambda_, coefficients):
    spatial, temporal, derivative = class6_coefficient_matrices(lambda_, coefficients)
    u, v = class6_vectors(spin, spin_x, lambda_, coefficients)
    v_x = -2 * (apply_matrix(temporal, spin_x)
                + apply_matrix(derivative, torch.linalg.cross(spin, spin_xx)))
    return apply_matrix(spatial, spin_t) - v_x + torch.linalg.cross(u, v)


def class6_training_curvature(spin, spin_x, spin_xx, spin_t, lambda_, coefficients):
    # Restores the original loss normalization after t_old=-2t.
    return -0.5 * class6_curvature(spin, spin_x, spin_xx, spin_t, lambda_, coefficients)


def sample_class6_time_derivative(spin, seed):
    """Independent tangent time derivative; spatial sampling is unchanged."""
    gen = torch.Generator(device="cpu").manual_seed(seed)
    raw = torch.randn(spin.shape, generator=gen, dtype=RDTYPE)
    raw = raw + 1j * torch.randn(spin.shape, generator=gen, dtype=RDTYPE)
    return -2 * (raw - spin * bilinear_dot(spin, raw)[:, None])


def class6_training_residual(spin, spin_x, spin_xx, spin_t, lambda_,
                             coefficients, alpha6, mode="off-shell"):
    """Normalized F-A E. For this ansatz it equals the on-shell curvature.

    A has no field dependence, so the A S_t terms cancel for arbitrary
    coefficients. Off-shell mode supplies S_t independently of the EOM.
    """
    curvature = class6_curvature(spin, spin_x, spin_xx, spin_t, lambda_, coefficients)
    if mode == "off-shell":
        spatial, _, _ = class6_coefficient_matrices(lambda_, coefficients)
        eom = spin_t - class6_eom_rhs(spin, spin_x, spin_xx, alpha6)
        return -0.5 * (curvature - apply_matrix(spatial, eom))
    if mode == "on-shell":
        return -0.5 * curvature
    raise ValueError(f"Unknown training mode: {mode}")


def class6_exact_coefficients(lambda_, alpha6):
    leading = spectral_coefficient(lambda_)
    return np.array([alpha6/(2*leading), -alpha6**2/(8*leading**3),
                     -leading**2, alpha6/2, -3*alpha6**2/(8*leading**2),
                     leading, alpha6/(2*leading), -alpha6**2/(8*leading**3)])


def rms(value):
    return float(torch.sqrt(torch.mean(torch.abs(value)**2)).detach().cpu())
