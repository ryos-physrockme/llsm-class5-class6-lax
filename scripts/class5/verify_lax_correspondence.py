#!/usr/bin/env python3
"""Exact checks of the Class-5 continuum Lax correspondence.

Requires SymPy. From the repository root run:
  python scripts/class5/verify_lax_correspondence.py
Writes results/class5/lax_correspondence/verification_results.json.
Use --output PATH to write a different report.

Sources (read separately; no network requests made by this script):
  de Leeuw, Fontanella, Nieto Garcia, arXiv:2506.13598v2,
    equations (1.1), (1.5), (1.8), (1.11), (1.12), (2.14), (2.23).
  Kameyama, Yoshida, arXiv:1405.4467v2,
    equations (4.1)--(4.9), (A.1), (A.4).
  https://github.com/ryos-physrockme/llsm-class5-class6-lax
    paper/sections/03a_spatial.tex, 03b_time.tex, 03c_eom.tex.

Conventions:
  alpha is the continuum deformation coupling; lam is the continuum
  spectral parameter. epsilon is the lattice spacing, u=lam/epsilon,
  theta=a_2/2=epsilon*alpha, with a_1=1 and a_3=0.
  S is the unit spin in the manuscript. T is the rotated spin,
  S=exp(-alpha*x*M/2)*T. All vector products are complex bilinear.
  AK=4*pi*L_K/sqrt(Lambda_K), where Lambda_K is the coupling (not a
  spectral parameter) in the Kameyama--Yoshida action.
  Their variables are n=diag(1,-i,-i)*T, C=alpha**2/(2*AK**2),
  z_K=2*i*AK*lam, and t_K=AK**2*t.

The reconstructed quantum L is explicitly distinguished from printed
Eq. (1.12). Agreement with its continuum leading coefficient does NOT
assert agreement with the uncorrected full printed quantum expression.
No claim concerning real forms, periodic gauge transformations, or
finite-time approximation to quantum evolution is made.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sympy as s

I = s.I
alpha, lam, x, AK, epsilon = s.symbols('alpha lam x AK epsilon', nonzero=True)
u, v, theta = s.symbols('u v theta')
a1, a2, a3 = s.symbols('a1 a2 a3')
Id = s.eye(2)
E = s.Matrix([[0, 1], [0, 0]])
sigma = [s.Matrix([[0, 1], [1, 0]]),
         s.Matrix([[0, -I], [I, 0]]), s.diag(1, -1)]
P = s.Matrix([[1, 0, 0, 0], [0, 0, 1, 0],
              [0, 1, 0, 0], [0, 0, 0, 1]])
K = s.Matrix([[0, 1, -1, 0], [0, 0, 0, 0],
              [0, 0, 0, 0], [0, 0, 0, 0]])
checks: list[dict[str, object]] = []


def check(name: str, expr: s.Expr | s.MatrixBase) -> None:
    values = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    residual = [s.simplify(s.expand(entry)) for entry in values]
    ok = all(entry == 0 for entry in residual)
    checks.append({'name': name, 'passed': ok})
    if not ok:
        raise AssertionError(f'{name}: {residual}')
    print(f'PASS  {name}')


def blocks(a, b, c, d):
    return s.BlockMatrix([[a, b], [c, d]]).as_explicit()


def quantum_L(uu, tt, e, f, h, k):
    """L = diag(1+theta*e,1) L_XXX [[1,-theta*h],[0,1]]."""
    identity, zero = s.eye(e.rows), s.zeros(e.rows)
    return (blocks(identity+tt*e, zero, zero, identity)
            * blocks(2*uu*identity+h, f, e, 2*uu*identity+k)
            * blocks(identity, -tt*h, zero, identity)).applyfunc(s.expand)


def printed_L(uu, tt, e, f, h, k):
    """Printed (1.12), with its a_2 coefficient replaced by tt.

    Keeping e, not f, in the inner upper-right parenthesis is intentional.
    """
    identity = s.eye(e.rows)
    a = 2*uu*identity+h
    return blocks((identity+tt*e)*a,
                  f-tt*a*h+tt*e*(e-tt*a*h),
                  e, 2*uu*identity+k-tt*e*h)


def pauli(vector):
    return sum((vector[j]*sigma[j] for j in range(3)), s.zeros(2))


def spatial_U(vector):
    p, q, z = vector[0]+I*vector[1], vector[0]-I*vector[1], vector[2]
    return s.Matrix([[1+z+2*alpha*lam*p, q-2*alpha*lam*(1+z)],
                     [p, 1-z]])/(4*lam)


def spin_linear_part(vector):
    return spatial_U(vector)-spatial_U(s.zeros(3, 1))


def time_V(vector, vector_x):
    return (-2*spin_linear_part(vector.cross(vector_x))
            +2*I*spatial_U(vector)**2-I*Id/(4*lam**2))


def traceless(matrix):
    return matrix-s.trace(matrix)*Id/2


def lower_symbol(matrix, rho):
    return s.Matrix(2, 2, lambda a, b:
                    s.trace(rho*matrix[2*a:2*a+2, 2*b:2*b+2]))


def spin_representation(two_j: int):
    jj, dim = s.Rational(two_j, 2), two_j+1
    sz = s.diag(*[jj-k for k in range(dim)])
    e = s.zeros(dim)
    for k in range(1, dim):
        m = jj-k
        e[k-1, k] = s.sqrt((jj-m)*(jj+m+1))
    return e, e.T, s.eye(dim)/2+sz, s.eye(dim)/2-sz


def embed_first_auxiliary(matrix, dim):
    """Embed in auxiliary 1 x auxiliary 2 x physical, identity on aux 2."""
    out = s.zeros(4*dim)
    for a in range(2):
        for b in range(2):
            for ap in range(2):
                for k in range(dim):
                    for l in range(dim):
                        out[(2*a+b)*dim+k, (2*ap+b)*dim+l] = matrix[a*dim+k, ap*dim+l]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output', type=Path,
        default=Path(__file__).resolve().parents[2]
                / 'results/class5/lax_correspondence/verification_results.json',
        help='Destination for the deterministic JSON verification report.',
    )
    args = parser.parse_args()

    # de Leeuw--Fontanella--Nieto Garcia (1.1) and (1.5), before fixing
    # a_1=1 and the a_3=0 local-basis representative.
    r_dfn = s.Matrix([
        [1+2*a1*u, a2*u, -a2*u, a2*a3*u**2],
        [0, 2*a1*u, 1, -a3*u],
        [0, 1, 2*a1*u, a3*u],
        [0, 0, 0, 1+2*a1*u],
    ])
    h_dfn = s.Matrix([
        [2, a2, -a2, 0],
        [0, 0, 2, a3],
        [0, 2, 0, -a3],
        [0, 0, 0, 2],
    ])

    e, f, h, k = E, E.T, s.diag(1, 0), s.diag(0, 1)
    l_direct = quantum_L(u, theta, e, f, h, k)
    r_class5 = 2*u*s.eye(4)+P+2*theta*u*K
    h_class5 = 2*P+2*theta*K
    check('DFN_general_R_reduces_to_manuscript_R',
          r_dfn.subs({a1: 1, a2: 2*theta, a3: 0})-r_class5)
    check('DFN_local_H_reduces_to_manuscript_H',
          h_dfn.subs({a2: 2*theta, a3: 0})-h_class5)
    check('twist_product_equals_fundamental_R', l_direct-r_class5)
    check('printed_upper_right_operator_difference_after_coefficient_correction',
          printed_L(u, theta, e, f, h, k)-l_direct
          +theta*s.kronecker_product(E, h))

    spin = s.Matrix(s.symbols('S1 S2 S3'))
    rho = (Id+pauli(spin))/2
    scaling = {u: lam/epsilon, theta: epsilon*alpha}
    direct_symbol = lower_symbol((l_direct/(2*u)).subs(scaling), rho)
    check('normalized_lower_symbol_exactly_I_plus_epsilon_U',
          direct_symbol-Id-epsilon*spatial_U(spin))
    half_symbol = lower_symbol((printed_L(u, theta, e, f, h, k)/(2*u)).subs(scaling), rho)
    check('coefficient_corrected_printed_L_differs_only_at_epsilon_squared',
          half_symbol-direct_symbol
          +epsilon**2*alpha*(1+spin[2])*E/(4*lam))
    literal = lower_symbol((printed_L(u, 2*theta, e, f, h, k)/(2*u)).subs(scaling), rho)
    literal_U = literal.applyfunc(lambda z: s.expand(z).coeff(epsilon, 1))
    k_symbol = s.Matrix([[spin[0]+I*spin[1], -(1+spin[2])], [0, 0]])/2
    check('uncorrected_printed_L_has_double_deformation_in_U',
          literal_U-spatial_U(spin)-alpha*k_symbol)

    for two_j in (1, 2):
        ee, ff, hh, kk = spin_representation(two_j)
        dim = two_j+1
        l1 = embed_first_auxiliary(quantum_L(u, theta, ee, ff, hh, kk), dim)
        l2 = s.kronecker_product(Id, quantum_L(v, theta, ee, ff, hh, kk))
        rr = s.kronecker_product(2*(u-v)*s.eye(4)+P+2*theta*(u-v)*K, s.eye(dim))
        check(f'exact_symbolic_RLL_spin_{s.Rational(two_j, 2)}', rr*l1*l2-l2*l1*rr)

    T = s.Matrix(s.symbols('T1 T2 T3'))
    Tx = s.Matrix(s.symbols('T1x T2x T3x'))
    Txx = s.Matrix(s.symbols('T1xx T2xx T3xx'))
    M = s.Matrix([[0, 0, -1], [0, 0, -I], [1, I, 0]])
    R = s.eye(3)-alpha*x*M/2+alpha**2*x**2*M**2/8
    check('nilpotence_M_cubed', M**3)
    check('field_rotation_complex_orthogonal', R.T*R-s.eye(3))
    S, Sx = R*T, R*(Tx-alpha*M*T/2)
    Sxx = R*(Txx-alpha*M*Tx+alpha**2*M**2*T/4)
    D = s.diag(1, -I, -I)
    n, nx, nxx = D*T, D*Tx, D*Txx
    check('sphere_to_hyperboloid_constraint',
          D.T*s.diag(-1, 1, 1)*D+s.eye(3))

    C, zK, beta = alpha**2/(2*AK**2), 2*I*AK*lam, -2/AK
    n0, n1, n2 = n
    np, nm = (n0+n1)/s.sqrt(2), (n0-n1)/s.sqrt(2)
    npx, nmx = (nx[0]+nx[1])/s.sqrt(2), (nx[0]-nx[1])/s.sqrt(2)
    T0, T1, T2 = I*sigma[1]/2, sigma[0]/2, sigma[2]/2
    Tp, Tm = (T0+T1)/s.sqrt(2), (T0-T1)/s.sqrt(2)
    # Direct transcription of Kameyama--Yoshida (4.7), with (4.8).
    UK = AK/zK*(-(np-C*zK**2*nm/2)*Tm+n2*T2-nm*Tp)
    VK = (beta/zK*(-((n2*npx-np*nx[2])-C*zK**2/2*(nm*nx[2]-n2*nmx))*Tm
                   +(nm*npx-np*nmx)*T2-(nm*nx[2]-n2*nmx)*Tp)
          +AK*beta/zK**2*(-(np+C*zK**2*nm/2)*Tm+n2*T2-nm*Tp))

    Q = s.Matrix([[0, I], [1, 0]])
    G = Id+alpha*(lam-x/2)*E
    J = G*Q.inv()
    expected_U = J*UK*J.inv()+s.diff(J, x)*J.inv()
    expected_V = AK**2*J*VK*J.inv()
    check('spatial_gauge_identity_all_entries', traceless(spatial_U(S))-expected_U)
    norm_minus_one = (T.T*T)[0]-1
    check('time_gauge_identity_modulo_unit_spin_all_entries',
          traceless(time_V(S, Sx))-expected_V-I*alpha*norm_minus_one*E/(2*lam))
    check('gauge_matrix_invertible', s.det(G)-1)

    Urot = pauli(T)/(4*lam)+alpha**2*lam*(T[0]+I*T[1])*E/4
    check('KY_spatial_matrix_in_rotated_spin_frame', Q.inv()*UK*Q-Urot)
    cross = T.cross(Tx)
    Urot_cross = pauli(cross)/(4*lam)+alpha**2*lam*(cross[0]+I*cross[1])*E/4
    Vrot = -2*Urot_cross-I*s.diff(Urot, lam)
    check('KY_time_matrix_in_rotated_spin_frame', AK**2*Q.inv()*VK*Q-Vrot)

    flow5 = -2*S.cross(Sxx+alpha*M*Sx)
    flowT = -2*T.cross(Txx-alpha**2*M**2*T/4)
    W = s.Matrix([n1*nxx[2]-n2*nxx[1], n0*nxx[2]-n2*nxx[0],
                  -n0*nxx[1]+n1*nxx[0]])
    Y = s.Matrix([(n0-n1)*n2, (n0-n1)*n2, (n0-n1)**2])
    check('field_rotation_of_Hamiltonian_equations', flow5-R*flowT)
    check('KY_equations_with_time_rescaling', D*flowT-AK**2*(2*W/AK**2-C*Y))

    # Trace checks are independent of the hand-written scalar curvature below.
    spin_x = s.Matrix(s.symbols('S1x S2x S3x'))
    spin_p, spin_px = spin[0]+I*spin[1], spin_x[0]+I*spin_x[1]
    scalar_v = (I*alpha*(spin_p*spin_x[2]-spin[2]*spin_px)/2
                +I*alpha**2*spin_p**2/4)
    check('spatial_scalar_part_from_trace',
          s.trace(spatial_U(spin))/2-1/(4*lam)-alpha*spin_p/4)
    check('time_scalar_part_from_trace_modulo_unit_spin',
          s.trace(time_V(spin, spin_x))/2-scalar_v
          -I*((spin.T*spin)[0]-1)/(8*lam**2))

    p, z, px, zx, pxx, zxx, pt = s.symbols('p z px zx pxx zxx pt')
    scalar_V_x = I*alpha*(p*zxx-z*pxx)/2+I*alpha**2*p*px/2
    Eplus = pt-2*I*(p*zxx-z*pxx+alpha*p*px)
    check('scalar_curvature_equals_alpha_Eplus_over_four',
          alpha*pt/4-scalar_V_x-alpha*Eplus/4)
    report = {'checks': checks, 'passed': len(checks),
              'scope': 'Exact algebra; unit-spin restriction is stated explicitly. '
                       'No equality of real or global boundary sectors asserted.'}
    destination = args.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'\n{len(checks)} exact checks passed. Results: {destination}')


if __name__ == '__main__':
    main()
