#!/usr/bin/env python3
"""Conditional shared-site correction for fundamental spin-1/2 chains.

Conventions: auxiliary factor 0, physical factors 1 and 2; h0=2P;
Lhat=I+eps*X+eps**2*Y+..., d_u Lhat=eps**2*Z+...;
h=h0+eps*h1+...; Ahat=-i Lhat^{-1}([h,Lhat]+d_u Lhat).

The symbolic generic calculation does NOT assume a Yang--Baxter solution.
It proves an identity with the order-eps**2 Sutherland residual retained.
The conditional consequence requires that residual to vanish and
P_plus*h1*P_plus=c*P_plus.  See research/shared_site_generalization.md.

Run from any directory.  Only SymPy and the Python standard library are
needed.  A failed exact check raises an exception and no success JSON is
written.  This program does not change any manuscript or prior results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
from typing import Callable

import sympy as sp

I = sp.I
ID2 = sp.eye(2)
ID4 = sp.eye(4)
P = sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0],
               [0, 1, 0, 0], [0, 0, 0, 1]])
PAULI = (sp.Matrix([[0, 1], [1, 0]]),
         sp.Matrix([[0, -I], [I, 0]]), sp.diag(1, -1))
BASE_COMMIT = "ac97ec22c21c5f0c8feb5ca3aaf73e1a72ea2325"


def comm(a: sp.MatrixBase, b: sp.MatrixBase) -> sp.Matrix:
    return a * b - b * a


def expand(matrix: sp.MatrixBase) -> sp.Matrix:
    return matrix.applyfunc(sp.expand)


def embed(operator: sp.MatrixBase, first: int, second: int) -> sp.Matrix:
    """Embed a 4x4 operator, respecting the ordered pair of tensor factors."""
    if operator.shape != (4, 4) or first == second or not {first, second} <= {0, 1, 2}:
        raise ValueError("Expected a 4x4 operator and two distinct factors in 0,1,2.")
    result = sp.zeros(8)
    for row in range(8):
        r = [(row >> k) & 1 for k in (2, 1, 0)]
        for col in range(8):
            c = [(col >> k) & 1 for k in (2, 1, 0)]
            if all(r[k] == c[k] for k in range(3) if k not in (first, second)):
                result[row, col] = operator[2*r[first]+r[second], 2*c[first]+c[second]]
    return result


def symbol(operator: sp.MatrixBase, rho: sp.MatrixBase) -> sp.Matrix:
    """Trace physical factor of auxiliary-first 4x4 operator; rho may be a tangent."""
    return sp.Matrix(2, 2, lambda a, b: sp.trace(rho * operator[2*a:2*a+2, 2*b:2*b+2]))


def two_site_symbol(operator: sp.MatrixBase, left: sp.MatrixBase,
                    right: sp.MatrixBase) -> sp.Matrix:
    """Bilinear contraction, also valid when either argument is a density tangent."""
    state = sp.kronecker_product(left, right)
    return sp.Matrix(2, 2, lambda a, b:
                     sp.expand(sp.trace(state * operator[4*a:4*a+4, 4*b:4*b+4])))


def overlap_data(x: sp.MatrixBase, y: sp.MatrixBase, z: sp.MatrixBase,
                 h1: sp.MatrixBase, rho: sp.MatrixBase,
                 drho: sp.MatrixBase) -> dict[str, sp.Matrix]:
    """Form operator products before taking symbols. No continuum V is an input."""
    xl, xr = embed(x, 0, 1), embed(x, 0, 2)
    yl, yr = embed(y, 0, 1), embed(y, 0, 2)
    h0, hfirst = 2 * embed(P, 1, 2), embed(h1, 1, 2)
    a1 = -I * comm(h0, xr)
    a2 = I*xr*comm(h0, xr) - I*comm(h0, yr) - I*comm(hfirst, xr) - I*embed(z, 0, 2)
    lower = lambda q: two_site_symbol(q, rho, rho)
    u, ux, ybar = symbol(x, rho), symbol(x, drho), symbol(y, rho)
    xi = expand(symbol(x*x, rho) - u*u)
    xi_x = expand(symbol(x*x, drho) - ux*u - u*ux)
    a1bar, a2bar = lower(a1), lower(a2)
    d2plus = lower(a1*xl) - a1bar*u
    d2minus = lower(xr*a1) - u*a1bar
    # Dplus varies at the right site, Dminus at the left site. Their
    # contributions ADD because the lattice equation subtracts Dminus.
    quadratic_taylor = (
        two_site_symbol(a1*xl, rho, drho)
        - two_site_symbol(a1, rho, drho)*u
        + two_site_symbol(xr*a1, drho, rho)
        - u*two_site_symbol(a1, drho, rho)
    )
    cubic = (lower(a2*xl + a1*yl - xr*a2 - yr*a1)
             - a2bar*u + u*a2bar - a1bar*ybar + ybar*a1bar)
    sutherland = (comm(h0, xr*xl + yr + yl) + comm(hfirst, xr+xl)
                  - embed(z, 0, 1) + embed(z, 0, 2))
    return {"U": u, "Xi": xi, "Xi_x": xi_x, "A1_symbol": a1bar,
            "A2_symbol": a2bar, "D2_plus": d2plus, "D2_minus": d2minus,
            "quadratic_taylor": quadratic_taylor, "cubic": cubic,
            "source": quadratic_taylor+cubic, "sutherland": sutherland,
            "sutherland_contraction": lower(xl*sutherland),
            "h1_commutator_symbol": lower(comm(hfirst, xr))}


class Checks:
    def __init__(self) -> None:
        self.results: dict[str, str] = {}

    def zero(self, name: str, matrix: sp.MatrixBase,
             reducer: Callable[[sp.Expr], sp.Expr] = sp.expand) -> None:
        result = matrix.applyfunc(reducer)
        if result != sp.zeros(*result.shape):
            raise AssertionError(f"{name}: nonzero residual {result}")
        self.results[name] = "exact_zero"
        print(f"{name}: exact_zero", flush=True)

    def nonzero(self, name: str, matrix: sp.MatrixBase) -> None:
        result = expand(matrix)
        if result == sp.zeros(*result.shape):
            raise AssertionError(f"{name}: negative control unexpectedly vanished")
        self.results[name] = "expected_nonzero"
        print(f"{name}: expected_nonzero", flush=True)


def generic_checks(checks: Checks) -> None:
    # Every entry is independent. This is a polynomial identity, not a
    # random numerical test or an evaluation at fitted coefficients.
    x = sp.Matrix(4, 4, sp.symbols("x0:16"))
    y = sp.Matrix(4, 4, sp.symbols("y0:16"))
    z = sp.Matrix(4, 4, sp.symbols("z0:16"))
    h = sp.Matrix(4, 4, sp.symbols("h0:16"))
    c, d, p, q = sp.symbols("c d p q")
    # This spans all h1 with scalar restriction to the spin-triplet sector.
    h1 = (h-P*h*P)/2 + c*ID4 + d*(ID4-P)
    plus = (ID4+P)/2
    rho = sp.diag(1, 0)
    drho = sp.Matrix([[0, p], [q, 0]])
    checks.zero("generic.triplet_condition", plus*h1*plus-c*plus)
    checks.zero("generic.projector_tangent", rho*drho+drho*rho-drho)
    data = overlap_data(x, y, z, h1, rho, drho)
    checks.zero("generic.coincident_A1", data["A1_symbol"])
    checks.zero("generic.quadratic_plus", data["D2_plus"]-2*I*data["Xi"])
    checks.zero("generic.quadratic_minus", data["D2_minus"]-2*I*data["Xi"])
    checks.zero("generic.quadratic_Taylor", data["quadratic_taylor"]-2*I*data["Xi_x"])
    checks.zero("generic.cubic_with_Sutherland_residual",
                data["cubic"]-2*I*comm(data["Xi"], data["U"])
                + I*data["sutherland_contraction"])
    checks.zero("generic.full_source_with_Sutherland_residual",
                data["source"]-2*I*(data["Xi_x"]+comm(data["Xi"], data["U"]))
                + I*data["sutherland_contraction"])
    checks.zero("generic.corrected_time_symbol",
                data["A2_symbol"]+2*I*data["Xi"]
                + I*data["h1_commutator_symbol"]+I*symbol(z, rho))
    # The cancellation of Y is an output, not an assumption Y=0.
    checks.zero("generic.cubic_independent_of_second_spatial_coefficient",
                expand(data["cubic"]) - expand(data["cubic"]).subs(dict.fromkeys(list(y), 0)))


def model_checks(checks: Checks) -> dict[str, str]:
    lam, alpha, a, b = sp.symbols("lambda alpha a b")
    spin = sp.Matrix(sp.symbols("s1:4"))
    spin_x = sp.Matrix(sp.symbols("v1:4"))
    rho = (ID2+sum((spin[k]*PAULI[k] for k in range(3)), sp.zeros(2)))/2
    drho = sum((spin_x[k]*PAULI[k] for k in range(3)), sp.zeros(2))/2
    ideal = sp.groebner([spin.dot(spin)-1, spin.dot(spin_x)],
                       *list(spin_x), *list(spin),
                       domain=sp.QQ_I.frac_field(lam, alpha, a, b))
    reduce_spin = lambda e: sp.factor(ideal.reduce(sp.expand(e))[1])
    k5 = sp.Matrix([[0, 1, -1, 0], [0, 0, 0, 0],
                    [0, 0, 0, 0], [0, 0, 0, 0]])
    k6 = sp.Matrix([[0, 1, 1, 0], [0, 0, 0, -1],
                    [0, 0, 0, -1], [0, 0, 0, 0]])
    x5 = P/(2*lam)+alpha*k5
    x6 = P/(2*lam)+alpha*lam*k6+alpha**2*lam**3*k6*k6
    xz = sp.Matrix([[a, 0, 0, 0], [0, 0, b, 0],
                    [0, b, 0, 0], [0, 0, 0, a]])
    # a=g*coth(2*g*lambda), b=g*csch(2*g*lambda);
    # hence a'=-2*b**2, b'=-2*a*b, a**2-b**2=g**2.
    zz = -2*b*b*xz.diff(a)-2*a*b*xz.diff(b)
    models = (
        ("class5", x5, sp.zeros(4), x5.diff(lam), 2*alpha*k5),
        ("class6", x6, alpha*k6/2+alpha**2*lam**2*k6*k6, x6.diff(lam), sp.zeros(4)),
        ("XXZ", xz, (a*a-b*b)*sp.diag(1, 0, 0, 1)/2, zz, sp.zeros(4)),
    )
    expressions = {}
    for name, x, y, z, h1 in models:
        data = overlap_data(x, y, z, h1, rho, drho)
        checks.zero(f"{name}.Sutherland_order_two", data["sutherland"])
        plus = (ID4+P)/2
        checks.zero(f"{name}.triplet_condition", plus*h1*plus)
        checks.zero(f"{name}.quadratic_cancellation", data["D2_plus"]-data["D2_minus"], reduce_spin)
        checks.zero(f"{name}.quadratic_value", data["D2_plus"]-2*I*data["Xi"], reduce_spin)
        checks.zero(f"{name}.quadratic_Taylor", data["quadratic_taylor"]-2*I*data["Xi_x"], reduce_spin)
        checks.zero(f"{name}.cubic_commutator", data["cubic"]-2*I*comm(data["Xi"], data["U"]), reduce_spin)
        checks.zero(f"{name}.complete_source", data["source"]-2*I*(data["Xi_x"]+comm(data["Xi"], data["U"])), reduce_spin)
        checks.zero(f"{name}.corrected_time_symbol", data["A2_symbol"]+2*I*data["Xi"]+I*data["h1_commutator_symbol"]+I*symbol(z, rho), reduce_spin)
        u = data["U"]
        gradient = -2*sum((spin.cross(spin_x)[k]*u.diff(spin[k]) for k in range(3)), sp.zeros(2))
        raw = gradient+data["A2_symbol"]
        if name == "class5":
            corrected = raw+2*I*data["Xi"]-I*ID2/(4*lam**2)
            expected = gradient+2*I*u*u-I*ID2/(4*lam**2)
        elif name == "class6":
            corrected = raw+2*I*data["Xi"]-I*ID2/(4*lam**2)
            u0 = u-ID2/(4*lam)
            expected = gradient-I*u0.diff(lam)
        else:
            zbar = symbol(z, rho)
            # A field-independent scalar fixes the traceless representative.
            corrected = raw+2*I*data["Xi"]+I*sp.trace(zbar)*ID2/2
            expected = gradient-I*(zbar-sp.trace(zbar)*ID2/2)
        checks.zero(f"{name}.final_time_matrix", corrected-expected, reduce_spin)
        if name == "XXZ":
            spin_xx = sp.Matrix(sp.symbols("w1:4"))
            spin_t = sp.Matrix(sp.symbols("t1:4"))
            flow = -2*spin.cross(spin_xx+(a*a-b*b)*spin[2]*sp.Matrix([0, 0, 1]))
            vx = sum((spin_x[k]*expected.diff(spin[k])
                      + spin_xx[k]*expected.diff(spin_x[k]) for k in range(3)), sp.zeros(2))
            ut = sum((spin_t[k]*u.diff(spin[k]) for k in range(3)), sp.zeros(2))
            eom_map = sum(((spin_t[k]-flow[k])*u.diff(spin[k]) for k in range(3)), sp.zeros(2))
            checks.zero("XXZ.off_shell_curvature_factorization", ut-vx+comm(u, expected)-eom_map, reduce_spin)
            expressions["XXZ_U"] = str(u.applyfunc(sp.factor))
            expressions["XXZ_V"] = str(expected.applyfunc(sp.factor))
            expressions["XXZ_flow"] = str(flow)
            expressions["XXZ_spectral_relations"] = "a=g*coth(2*g*lambda); b=g*csch(2*g*lambda); a**2-b**2=g**2"
    return expressions


def xxz_scaling_checks(checks: Checks) -> None:
    """Independent regular six-vertex R-matrix and weak-anisotropy scaling."""
    q, t, v = sp.symbols("q t v", nonzero=True)

    def r_exponential(w: sp.Expr) -> sp.Matrix:
        diagonal = (q*w-1/(q*w))/(q-1/q)
        middle = (w-1/w)/(q-1/q)
        return sp.Matrix([[diagonal, 0, 0, 0], [0, middle, 1, 0],
                          [0, 1, middle, 0], [0, 0, 0, diagonal]])

    checks.zero("XXZ.quantum_regularity", r_exponential(1)-P, sp.cancel)
    r12 = embed(r_exponential(t/v), 0, 1)
    r13 = embed(r_exponential(t), 0, 2)
    r23 = embed(r_exponential(v), 1, 2)
    checks.zero("XXZ.quantum_Yang_Baxter", r12*r13*r23-r23*r13*r12, sp.cancel)
    # q=exp(eta), t=exp(2*eta*u). For eta=eps*g, u=lambda/eps,
    # the normalization is sinh(eta)/sinh(2*eta*u).
    eps, g, lam = sp.symbols("epsilon g lambda", nonzero=True)
    theta = 2*g*lam
    diagonal = sp.sinh(theta+eps*g)/sp.sinh(theta)
    off = sp.sinh(eps*g)/sp.sinh(theta)
    checks.zero("XXZ.normalized_L_diagonal_first", sp.Matrix([
        diagonal.diff(eps).subs(eps, 0)-g*sp.coth(theta)]), sp.simplify)
    checks.zero("XXZ.normalized_L_diagonal_second", sp.Matrix([
        diagonal.diff(eps, 2).subs(eps, 0)/2-g*g/2]), sp.simplify)
    checks.zero("XXZ.normalized_L_off_diagonal_first", sp.Matrix([
        off.diff(eps).subs(eps, 0)-g/sp.sinh(theta)]), sp.simplify)
    a, b = g*sp.coth(theta), g/sp.sinh(theta)
    checks.zero("XXZ.spectral_derivatives", sp.Matrix([
        a.diff(lam)+2*b*b, b.diff(lam)+2*a*b, a*a-b*b-g*g]), sp.simplify)
    # h=P*d_u R(0): parallel entry 2*eta*coth(eta), exchange
    # entry 2*eta/sinh(eta). Subtract 2P before expanding.
    eta = sp.Symbol("eta")
    parallel_second = sp.series(2*eta*sp.coth(eta), eta, 0, 3).removeO().coeff(eta, 2)
    exchange_second = sp.series(2*eta/sp.sinh(eta), eta, 0, 3).removeO().coeff(eta, 2)
    checks.zero("XXZ.Hamiltonian_order_two", sp.Matrix([
        parallel_second-sp.Rational(2, 3), exchange_second+sp.Rational(1, 3)]))
    spin = sp.Matrix(sp.symbols("r1:4"))
    rho = (ID2+sum((spin[k]*PAULI[k] for k in range(3)), sp.zeros(2)))/2
    h2 = sp.Matrix([[sp.Rational(2,3), 0, 0, 0], [0, 0, -sp.Rational(1,3), 0],
                    [0, -sp.Rational(1,3), 0, 0], [0, 0, 0, sp.Rational(2,3)]])
    energy = sp.trace(sp.kronecker_product(rho, rho)*h2)
    ideal = sp.groebner([spin.dot(spin)-1], *list(spin), domain=sp.QQ)
    checks.zero("XXZ.Hamiltonian_anisotropy", sp.Matrix([
        energy-sp.Rational(1,6)-spin[2]**2/2]),
        lambda e: ideal.reduce(sp.expand(e))[1])


def negative_control(checks: Checks) -> dict[str, str]:
    # An arbitrary trial X that is NOT a spatial Lax operator of h0=2P.
    # This is a counterexample only to dropping the Sutherland condition.
    x = sp.kronecker_product(PAULI[0], PAULI[0])+sp.kronecker_product(PAULI[2], PAULI[2])
    data = overlap_data(x, sp.zeros(4), sp.zeros(4), sp.zeros(4),
                        sp.diag(1, 0), sp.zeros(2))
    difference = expand(data["source"]-2*I*comm(data["Xi"], data["U"]))
    checks.nonzero("negative.Sutherland_condition_violated", data["sutherland"])
    checks.nonzero("negative.unconditional_formula_fails", difference)
    checks.zero("negative.explicit_residual", difference+4*I*PAULI[2])
    checks.zero("negative.residual_explained_by_Sutherland", difference+I*data["sutherland_contraction"])
    return {"X": str(x), "Z": "0", "h1": "0", "rho": "diag(1,0)",
            "source_minus_covariant_second_moment": str(difference),
            "interpretation": "Arbitrary non-Lax trial operator; not a counterexample among integrable chains."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(__file__).resolve().parents[2]/"results"/"generalization"/"shared_site_identity.json"
    parser.add_argument("--out", type=Path, default=default)
    args = parser.parse_args()
    checks = Checks()
    generic_checks(checks)
    expressions = model_checks(checks)
    xxz_scaling_checks(checks)
    control = negative_control(checks)
    report = {
        "schema_version": 1, "base_commit": BASE_COMMIT,
        "python_version": platform.python_version(), "sympy_version": sp.__version__,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "status": "conditional_algebraic_identity_verified",
        "scope": "two-dimensional auxiliary and physical spaces; pure product symbols; h0=2P; scalar triplet restriction of h1; order-two Sutherland relation",
        "generic_proof": "Independent symbolic entries of X,Y,Z and an arbitrary admissible h1, at rho=diag(1,0) with the general projector tangent; constant simultaneous physical similarity gives any rank-one rho.",
        "number_of_checks": len(checks.results), "checks": checks.results,
        "negative_control": control, "expressions": expressions,
        "notation": {"I": "imaginary unit", "s1,s2,s3": "Cartesian unit-spin components",
                     "v1,v2,v3": "first spatial derivatives of s1,s2,s3",
                     "w1,w2,w3": "second spatial derivatives of s1,s2,s3",
                     "g": "continuum XXZ anisotropy scale; eta=epsilon*g",
                     "lambda": "continuum spectral parameter; u=lambda/epsilon",
                     "a,b": "spectral functions defined in XXZ_spectral_relations",
                     "X,Y,Z": "coefficients defined in the script module docstring"},
        "not_established": ["priority or absence from literature", "arbitrary spin representations",
                            "higher Hamiltonian flows", "finite-time quantum convergence",
                            "the repository's full pre-existing verification suite"],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"All {len(checks.results)} checks passed. Report: {args.out}")


if __name__ == "__main__":
    main()
