#!/usr/bin/env python3
"""Exact comparison of the manuscript Class 6 Lax pair with the known
rational eleven-vertex Landau--Lifshitz pair.

Literature conventions used here are the explicit rational LL formulas
quoted in Atalikov--Zotov, arXiv:2010.14297, Sec. 4, which refer back to
the eleven-vertex construction of Levin--Olshanetsky--Zotov,
arXiv:1406.2995.

Define beta by beta**2 = -2*alpha6, z = beta*lambda and
Sigma = [[S3,Splus],[Sminus,-S3]].  The comparison checks

  U6_traceless = (beta/4) U_11v(z,Sigma)^T,
  V6           = (i*beta**2/2) V_11v(z,Sigma)^T,

with the LL spatial coefficient k=-4/beta.  The time variables then obey
 t_11v = (i*beta**2/2) t = -i*alpha6*t.

The scalar I/(4*lambda) in the manuscript U6 has zero curvature and is
therefore kept separate from the comparison.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import sympy as sp

I = sp.I
ID2 = sp.eye(2)
BASE_COMMIT = "ac97ec22c21c5f0c8feb5ca3aaf73e1a72ea2325"


def comm(a: sp.MatrixBase, b: sp.MatrixBase) -> sp.Matrix:
    return a*b-b*a


def zero(name: str, matrix: sp.MatrixBase, checks: dict[str, str]) -> None:
    result = matrix.applyfunc(lambda entry: sp.simplify(sp.expand(entry)))
    if any(entry != 0 for entry in result):
        raise AssertionError(f"{name}: {result}")
    checks[name] = "exact_zero"
    print(f"{name}: exact_zero")


def u_eleven_vertex(z: sp.Expr, s: sp.MatrixBase) -> sp.Matrix:
    """Atalikov--Zotov Eq. (4.3), equivalent to tr_2(r_12(z) S_2)."""
    s11, s12, s21 = s[0, 0], s[0, 1], s[1, 0]
    return sp.Matrix([
        [s11/z-z*s12, s12/z],
        [s21/z-2*z*s11-z**3*s12, -s11/z+z*s12],
    ])


def m_top(z: sp.Expr, s: sp.MatrixBase) -> sp.Matrix:
    """Atalikov--Zotov Eq. (4.5)."""
    s11, s12 = s[0, 0], s[0, 1]
    return sp.Matrix([[s12, 0], [2*s11+2*z**2*s12, -s12]])


def main() -> None:
    beta, lam = sp.symbols("beta lambda", nonzero=True)
    alpha6 = -beta**2/2
    sp_, sm_, s3 = sp.symbols("Splus Sminus S3")
    spx, smx, s3x = sp.symbols("Splus_x Sminus_x S3_x")
    spxx, smxx, s3xx = sp.symbols("Splus_xx Sminus_xx S3_xx")

    sigma = sp.Matrix([[s3, sp_], [sm_, -s3]])
    sigma_x = sp.Matrix([[s3x, spx], [smx, -s3x]])
    sigma_xx = sp.Matrix([[s3xx, spxx], [smxx, -s3xx]])
    z = beta*lam
    k = -4/beta
    q_time = I*beta**2/2  # = -i*alpha6

    u11 = u_eleven_vertex(z, sigma)
    # The eigenvalues of Sigma are +/-1 on the unit-spin constraint, so the
    # literature formula v=-k [Sigma,Sigma_x]/4 applies.
    v_aux = -k*comm(sigma, sigma_x)/4
    v11 = sp.Rational(1, 2)*(u11/z+2*m_top(z, sigma)-u_eleven_vertex(z, v_aux))

    u6_tr = sp.Matrix([
        [s3/(4*lam)+alpha6*lam*sp_/2,
         sm_/(4*lam)+alpha6*lam*s3-alpha6**2*lam**3*sp_],
        [sp_/(4*lam), -s3/(4*lam)-alpha6*lam*sp_/2],
    ])
    v6 = I/(4*lam**2)*sp.Matrix([
        [s3-2*alpha6*lam**2*sp_+lam*(sm_*spx-sp_*smx)
         +4*alpha6*lam**3*(sp_*s3x-s3*spx),
         sm_-4*alpha6*lam**2*s3+12*alpha6**2*lam**4*sp_
         +2*lam*(s3*smx-sm_*s3x)
         +4*alpha6*lam**3*(sm_*spx-sp_*smx)
         +8*alpha6**2*lam**5*(s3*spx-sp_*s3x)],
        [sp_+2*lam*(sp_*s3x-s3*spx),
         -s3+2*alpha6*lam**2*sp_-lam*(sm_*spx-sp_*smx)
         -4*alpha6*lam**3*(sp_*s3x-s3*spx)],
    ])

    checks: dict[str, str] = {}
    zero("unit_spin_characteristic",
         sigma*sigma-(s3**2+sp_*sm_)*ID2, checks)
    zero("spectral_and_coupling_scaling",
         sp.Matrix([beta**2+2*alpha6, k*beta+4,
                    q_time+I*alpha6, beta*k/4+1]), checks)
    zero("spatial_Lax_dictionary", u6_tr-beta*u11.T/4, checks)
    zero("time_Lax_dictionary", v6-q_time*v11.T, checks)

    # Rational eleven-vertex anisotropy matrix, Atalikov--Zotov Eq. (4.1).
    j_sigma = -sp.Matrix([[sp_, 0], [2*s3, -sp_]])
    alpha0 = k**2/8  # eigenvalue squared is one
    flow_11 = comm(j_sigma, sigma)-alpha0*comm(sigma, sigma_xx)
    flow_mapped = q_time*flow_11

    splus_t = 2*I*(sp_*s3xx-s3*spxx+alpha6*sp_**2)
    sminus_t = 2*I*(s3*smxx-sm_*s3xx+2*alpha6*s3**2-alpha6*sp_*sm_)
    s3_t = I*(sm_*spxx-sp_*smxx-2*alpha6*sp_*s3)
    flow6 = sp.Matrix([[s3_t, splus_t], [sminus_t, -s3_t]])
    zero("equation_of_motion_dictionary", flow_mapped-flow6, checks)

    # The two zero-curvature conventions match after transpose and the
    # coordinate/time rescalings: U6=c0 U11^T, V6=q V11^T,
    # tau=q t and c0*k=-1.
    c0 = beta/4
    zero("zero_curvature_convention_factors",
         sp.Matrix([c0*k+1, q_time-c0*(-k)*q_time]), checks)

    result = {
        "schema_version": 1,
        "base_commit": BASE_COMMIT,
        "status": "exact_eleven_vertex_dictionary_verified",
        "python_version": platform.python_version(),
        "sympy_version": sp.__version__,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "number_of_checks": len(checks),
        "checks": checks,
        "dictionary": {
            "coupling": "beta^2 = -2 alpha6",
            "spectral_parameter": "z = beta lambda",
            "spatial_Lax": "U6_traceless = (beta/4) U_11v(z,Sigma)^T",
            "LL_spatial_coefficient": "k = -4/beta",
            "time_coordinate": "t_11v = (i beta^2/2) t = -i alpha6 t",
            "time_Lax": "V6 = (i beta^2/2) V_11v(z,Sigma)^T = -i alpha6 V_11v^T",
            "spin_matrix": "Sigma = [[S3,Splus],[Sminus,-S3]], Sigma^2=I on Splus*Sminus+S3^2=1",
            "scalar_part": "manuscript U6 = U6_traceless + I/(4 lambda); the scalar has zero curvature"
        },
        "sources": [
            "Levin, Olshanetsky, Zotov, Nucl. Phys. B 887 (2014) 400-422, arXiv:1406.2995",
            "Atalikov, Zotov, J. Geom. Phys. 164 (2021) 104161, arXiv:2010.14297, Sec. 4, Eqs. (4.1)-(4.5)"
        ],
        "scope": "algebraic convention dictionary for the known rational eleven-vertex Landau-Lifshitz U-V pair and the manuscript Class 6 pair; no novelty claim for the eleven-vertex hierarchy"
    }
    out = Path(__file__).resolve().parents[2]/"results"/"class6"/"eleven_vertex_comparison.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"All {len(checks)} checks passed. Report: {out}")


if __name__ == "__main__":
    main()
