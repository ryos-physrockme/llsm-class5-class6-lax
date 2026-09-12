"""Utilities for the Class-5 deformed Landau--Lifshitz Lax search.

All vector products and dot products are complex bilinear (no conjugation),
matching the complexified so(3) algebra and the constraint S^T S = 1.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import torch

CDTYPE = torch.complex128
RDTYPE = torch.float64


def class5_matrix_torch(device: str | torch.device = "cpu") -> torch.Tensor:
    return torch.tensor(
        [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0j], [1.0, -1.0j, 0.0]],
        dtype=CDTYPE,
        device=device,
    )


def class5_matrix_numpy() -> np.ndarray:
    return np.array(
        [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0j], [1.0, -1.0j, 0.0]],
        dtype=np.complex128,
    )


def bilinear_dot(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    return torch.sum(a * b, dim=-1)


def apply_matrix(mat: torch.Tensor, vec: torch.Tensor) -> torch.Tensor:
    return torch.einsum("ij,...j->...i", mat, vec)


def sample_complex_jets(
    n_samples: int,
    seed: int,
    device: str | torch.device = "cpu",
    scale_x: float = 0.8,
    scale_xx: float = 0.8,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Sample (S, S_x, S_xx) satisfying S^T S=1 and differentiated constraints.

    The complex quadric is sampled by normalizing generic complex Gaussian vectors
    with the bilinear norm. Samples too close to the null cone are rejected.
    """
    gen = torch.Generator(device="cpu")
    gen.manual_seed(seed)
    out_s: list[torch.Tensor] = []
    out_x: list[torch.Tensor] = []
    out_xx: list[torch.Tensor] = []
    remaining = n_samples
    while remaining > 0:
        batch = max(2 * remaining, 256)
        z = torch.randn(batch, 3, generator=gen, dtype=RDTYPE)
        z = z + 1.0j * torch.randn(batch, 3, generator=gen, dtype=RDTYPE)
        q = bilinear_dot(z, z)
        mask = torch.abs(q) > 0.25
        z = z[mask]
        q = q[mask]
        if z.shape[0] == 0:
            continue
        s = z / torch.sqrt(q)[:, None]
        r1 = torch.randn(s.shape[0], 3, generator=gen, dtype=RDTYPE)
        r1 = r1 + 1.0j * torch.randn(s.shape[0], 3, generator=gen, dtype=RDTYPE)
        sx = r1 - s * bilinear_dot(s, r1)[:, None]
        sx = scale_x * sx
        r2 = torch.randn(s.shape[0], 3, generator=gen, dtype=RDTYPE)
        r2 = r2 + 1.0j * torch.randn(s.shape[0], 3, generator=gen, dtype=RDTYPE)
        tangent = r2 - s * bilinear_dot(s, r2)[:, None]
        sxx = scale_xx * tangent - s * bilinear_dot(sx, sx)[:, None]
        take = min(remaining, s.shape[0])
        out_s.append(s[:take])
        out_x.append(sx[:take])
        out_xx.append(sxx[:take])
        remaining -= take
    return (
        torch.cat(out_s).to(device=device, dtype=CDTYPE),
        torch.cat(out_x).to(device=device, dtype=CDTYPE),
        torch.cat(out_xx).to(device=device, dtype=CDTYPE),
    )


def transformed_eom_rhs(
    s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor, k: float
) -> torch.Tensor:
    m = class5_matrix_torch(s.device)
    p = m @ m
    return torch.linalg.cross(s, sxx - k * apply_matrix(p, s), dim=-1)


def transformed_curvature(
    s: torch.Tensor,
    sx: torch.Tensor,
    sxx: torch.Tensor,
    st: torch.Tensor,
    a: float,
    coeffs: torch.Tensor,
) -> torch.Tensor:
    """Curvature for a five-coefficient ansatz.

    coeffs=(b,d,e,g,h),
      u=(a I+b P)S,
      v=(d I+e P)S+(g I+h P)(S x S_x).
    """
    b, d, e, g, h = coeffs
    eye = torch.eye(3, dtype=CDTYPE, device=s.device)
    m = class5_matrix_torch(s.device)
    p = m @ m
    amat = a * eye + b.to(CDTYPE) * p
    bmat = d.to(CDTYPE) * eye + e.to(CDTYPE) * p
    cmat = g.to(CDTYPE) * eye + h.to(CDTYPE) * p
    u = apply_matrix(amat, s)
    ut = apply_matrix(amat, st)
    cross_s_sx = torch.linalg.cross(s, sx, dim=-1)
    v = apply_matrix(bmat, s) + apply_matrix(cmat, cross_s_sx)
    vx = apply_matrix(bmat, sx) + apply_matrix(
        cmat, torch.linalg.cross(s, sxx, dim=-1)
    )
    return ut - vx + torch.linalg.cross(u, v, dim=-1)


def exact_coefficients(a: float, k: float) -> np.ndarray:
    return np.array([k / (2.0 * a), -a * a, k / 2.0, a, k / (2.0 * a)])


def transformed_exact_certificate(
    s: torch.Tensor,
    sx: torch.Tensor,
    sxx: torch.Tensor,
    st: torch.Tensor,
    a: float,
    k: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    coeffs = torch.tensor(exact_coefficients(a, k), dtype=RDTYPE, device=s.device)
    curvature = transformed_curvature(s, sx, sxx, st, a, coeffs)
    eye = torch.eye(3, dtype=CDTYPE, device=s.device)
    m = class5_matrix_torch(s.device)
    p = m @ m
    amat = a * eye + (k / (2.0 * a)) * p
    eom = st - transformed_eom_rhs(s, sx, sxx, k)
    certificate = apply_matrix(amat, eom)
    return curvature, certificate


def original_exact_curvature_and_certificate(
    s: torch.Tensor,
    sx: torch.Tensor,
    sxx: torch.Tensor,
    st: torch.Tensor,
    a: float,
    alpha: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Exact Class-5 curvature and A times the original EOM residual.

    Spatial and temporal Lax vectors are
      u=A S + c m,
      v=B S + A[S x (S_x-c M S)],
    where c=alpha/2, P=M^2=m m^T,
      A=a I+c^2/(2a)P, B=-a^2 I+c^2/2 P.
    """
    c = alpha / 2.0
    k = c * c
    eye = torch.eye(3, dtype=CDTYPE, device=s.device)
    mmat = class5_matrix_torch(s.device)
    p = mmat @ mmat
    mvec = torch.tensor([-1.0j, -1.0, 0.0], dtype=CDTYPE, device=s.device)
    amat = a * eye + (k / (2.0 * a)) * p
    bmat = -a * a * eye + (k / 2.0) * p
    ds = sx - c * apply_matrix(mmat, s)
    dsx = sxx - c * apply_matrix(mmat, sx)
    u = apply_matrix(amat, s) + c * mvec
    ut = apply_matrix(amat, st)
    v = apply_matrix(bmat, s) + apply_matrix(
        amat, torch.linalg.cross(s, ds, dim=-1)
    )
    vx = apply_matrix(bmat, sx) + apply_matrix(
        amat,
        torch.linalg.cross(sx, ds, dim=-1)
        + torch.linalg.cross(s, dsx, dim=-1),
    )
    curvature = ut - vx + torch.linalg.cross(u, v, dim=-1)
    rhs = torch.linalg.cross(
        s, sxx - alpha * apply_matrix(mmat, sx), dim=-1
    )
    certificate = apply_matrix(amat, st - rhs)
    return curvature, certificate


def rms(z: torch.Tensor) -> float:
    return float(torch.sqrt(torch.mean(torch.abs(z) ** 2)).detach().cpu())


def max_abs(z: torch.Tensor) -> float:
    return float(torch.max(torch.abs(z)).detach().cpu())


@dataclass(frozen=True)
class ConstraintErrors:
    s2: float
    s_sx: float
    s_sxx_plus_sx2: float


def constraint_errors(s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor) -> ConstraintErrors:
    return ConstraintErrors(
        s2=max_abs(bilinear_dot(s, s) - 1.0),
        s_sx=max_abs(bilinear_dot(s, sx)),
        s_sxx_plus_sx2=max_abs(bilinear_dot(s, sxx) + bilinear_dot(sx, sx)),
    )


def original_trainable_curvature(
    s: torch.Tensor,
    sx: torch.Tensor,
    sxx: torch.Tensor,
    st: torch.Tensor,
    a: float,
    coeffs: torch.Tensor,
) -> torch.Tensor:
    """Curvature for a direct Class-5 ansatz in the original variables.

    coeffs=(b,d,e,g,h,r,j),
      u=(a I+bP)S+r m,
      v=(d I+eP)S+(g I+hP)[S x (S_x-j M S)].
    The optimizer is not told that r=j=alpha/2 or that the two maps coincide.
    """
    b, d, e, g, h, r, j = coeffs
    eye = torch.eye(3, dtype=CDTYPE, device=s.device)
    mmat = class5_matrix_torch(s.device)
    p = mmat @ mmat
    mvec = torch.tensor([-1.0j, -1.0, 0.0], dtype=CDTYPE, device=s.device)
    amat = a * eye + b.to(CDTYPE) * p
    bmat = d.to(CDTYPE) * eye + e.to(CDTYPE) * p
    cmat = g.to(CDTYPE) * eye + h.to(CDTYPE) * p
    ds = sx - j.to(CDTYPE) * apply_matrix(mmat, s)
    dsx = sxx - j.to(CDTYPE) * apply_matrix(mmat, sx)
    u = apply_matrix(amat, s) + r.to(CDTYPE) * mvec
    ut = apply_matrix(amat, st)
    v = apply_matrix(bmat, s) + apply_matrix(
        cmat, torch.linalg.cross(s, ds, dim=-1)
    )
    vx = apply_matrix(bmat, sx) + apply_matrix(
        cmat,
        torch.linalg.cross(sx, ds, dim=-1)
        + torch.linalg.cross(s, dsx, dim=-1),
    )
    return ut - vx + torch.linalg.cross(u, v, dim=-1)


def original_eom_rhs(
    s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor, alpha: float
) -> torch.Tensor:
    mmat = class5_matrix_torch(s.device)
    return torch.linalg.cross(s, sxx - alpha * apply_matrix(mmat, sx), dim=-1)


def original_exact_coefficients(a: float, alpha: float) -> np.ndarray:
    c = alpha / 2.0
    k = c * c
    return np.array([k / (2.0 * a), -a * a, k / 2.0, a, k / (2.0 * a), c, c])


def class6_matrix_torch(device: str | torch.device = "cpu") -> torch.Tensor:
    return torch.tensor(
        [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0j], [1.0, -1.0j, 0.0]],
        dtype=CDTYPE,
        device=device,
    )


def class6_eom_rhs(s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor, alpha: float) -> torch.Tensor:
    nmat = class6_matrix_torch(s.device)
    return torch.linalg.cross(s, sxx - alpha * apply_matrix(nmat, s), dim=-1)


def class6_curvature(
    s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor, st: torch.Tensor,
    a: float, coeffs: torch.Tensor,
) -> torch.Tensor:
    """Class-6 curvature for polynomial-in-N maps.

    coeffs=(b,c,d,e,f,g,h,j):
      A=aI+bN+cN^2, B=dI+eN+fN^2, C=gI+hN+jN^2,
      u=A S, v=B S+C(S x S_x).
    """
    b,c,d,e,f,g,h,j=coeffs
    eye=torch.eye(3,dtype=CDTYPE,device=s.device)
    nmat=class6_matrix_torch(s.device); n2=nmat@nmat
    amat=a*eye+b.to(CDTYPE)*nmat+c.to(CDTYPE)*n2
    bmat=d.to(CDTYPE)*eye+e.to(CDTYPE)*nmat+f.to(CDTYPE)*n2
    cmat=g.to(CDTYPE)*eye+h.to(CDTYPE)*nmat+j.to(CDTYPE)*n2
    u=apply_matrix(amat,s); ut=apply_matrix(amat,st)
    v=apply_matrix(bmat,s)+apply_matrix(cmat,torch.linalg.cross(s,sx,dim=-1))
    vx=apply_matrix(bmat,sx)+apply_matrix(cmat,torch.linalg.cross(s,sxx,dim=-1))
    return ut-vx+torch.linalg.cross(u,v,dim=-1)


def class6_exact_coefficients(a: float, alpha: float) -> np.ndarray:
    return np.array([
        alpha/(2*a), -alpha**2/(8*a**3),
        -a*a, alpha/2, -3*alpha**2/(8*a*a),
        a, alpha/(2*a), -alpha**2/(8*a**3),
    ])


def class6_exact_certificate(
    s: torch.Tensor, sx: torch.Tensor, sxx: torch.Tensor, st: torch.Tensor,
    a: float, alpha: float,
) -> tuple[torch.Tensor,torch.Tensor]:
    coeffs=torch.tensor(class6_exact_coefficients(a,alpha),dtype=RDTYPE,device=s.device)
    f=class6_curvature(s,sx,sxx,st,a,coeffs)
    eye=torch.eye(3,dtype=CDTYPE,device=s.device)
    nmat=class6_matrix_torch(s.device); n2=nmat@nmat
    amat=a*eye+(alpha/(2*a))*nmat-(alpha**2/(8*a**3))*n2
    eom=st-class6_eom_rhs(s,sx,sxx,alpha)
    return f,apply_matrix(amat,eom)
