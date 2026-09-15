#!/usr/bin/env python3
"""Exact Hamiltonian/EOM checks for the Class-5/null-like correspondence."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sympy as s

I = s.I
alpha, x, AK = s.symbols('alpha x AK', nonzero=True)
checks: list[dict[str, object]] = []


def check(name: str, expr: s.Expr | s.MatrixBase) -> None:
    values = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    residual = [s.simplify(s.expand(v)) for v in values]
    passed = all(v == 0 for v in residual)
    checks.append({'name': name, 'passed': passed})
    if not passed:
        raise AssertionError(f'{name}: {residual}')
    print(f'PASS  {name}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output', type=Path,
        default=Path(__file__).resolve().parents[2]
                / 'results/class5/lax_correspondence/hamiltonian_verification.json',
    )
    args = parser.parse_args()

    T = s.Matrix(s.symbols('T1 T2 T3'))
    Tx = s.Matrix(s.symbols('T1x T2x T3x'))
    Txx = s.Matrix(s.symbols('T1xx T2xx T3xx'))

    M = s.Matrix([[0, 0, -1], [0, 0, -I], [1, I, 0]])
    R = s.eye(3)-alpha*x*M/2+alpha**2*x**2*M**2/8
    check('M_cubed_vanishes', M**3)
    check('field_rotation_is_complex_orthogonal', R.T*R-s.eye(3))

    S = R*T
    Sx = R*(Tx-alpha*M*T/2)
    h5 = -s.Rational(1, 2)*(Sx.T*Sx)[0] + alpha*(S.T*M*Sx)[0]/2
    hrot = (-s.Rational(1, 2)*(Tx.T*Tx)[0]
            -alpha**2*(T.T*M**2*T)[0]/8)
    check('Class5_Hamiltonian_density_after_field_rotation', h5-hrot)

    D = s.diag(1, -I, -I)
    eta = s.diag(1, -1, -1)
    n, nx, nxx = D*T, D*Tx, D*Txx
    check('sphere_to_hyperboloid_metric_map', D.T*eta*D-s.eye(3))

    n0, n1, n2 = n
    nm = (n0-n1)/s.sqrt(2)
    check('Class5_anisotropy_to_null_component', (T.T*M**2*T)[0]+2*nm**2)

    C = alpha**2/(2*AK**2)
    hK = C*nm**2-(nx.T*eta*nx)[0]/AK**2
    check('Hamiltonian_density_dictionary', hrot-AK**2*hK/2)

    flowT = -2*T.cross(Txx-alpha**2*M**2*T/4)
    W = s.Matrix([
        n1*nxx[2]-n2*nxx[1],
        n0*nxx[2]-n2*nxx[0],
        -n0*nxx[1]+n1*nxx[0],
    ])
    Y = s.Matrix([(n0-n1)*n2, (n0-n1)*n2, (n0-n1)**2])
    check('KY_equations_with_time_rescaling', D*flowT-AK**2*(2*W/AK**2-C*Y))

    report = {
        'checks': checks,
        'passed': len(checks),
        'scope': (
            'Exact local complexified Hamiltonian-density and EOM mapping. '
            'Real sections and periodic global phase spaces are not identified.'
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'\n{len(checks)} exact checks passed. Results: {args.output}')


if __name__ == '__main__':
    main()
