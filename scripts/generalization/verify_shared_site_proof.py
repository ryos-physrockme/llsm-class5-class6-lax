#!/usr/bin/env python3
"""Independent noncommutative verification of the shared-site correction.

The physical local space is C^2. Auxiliary matrix entries are free,
noncommuting symbols, so no auxiliary dimension or matrix polynomial
identity is used. The physical basis is (00,01,10,11). X,Y,Z are the
coefficients in Lhat=1+epsilon*X+epsilon**2*Y+... and
partial_u Lhat=epsilon**2*Z+.... The bond is h=2P+epsilon*h1+....

Run: python scripts/generalization/verify_shared_site_proof.py
A report is written only after every exact assertion has passed.
This script is independent of verify_shared_site_identity.py.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp

I = sp.I
P = sp.Matrix([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]])
RHO = sp.diag(1,0)
PAULI = (sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-I],[I,0]]), sp.diag(1,-1))
BASE_COMMIT = '722b5bf14bb64ca895375ab05569e11dd69f9712'


def comm(a, b):
    return a*b-b*a


def physical_embeddings(x: sp.MatrixBase):
    """Physical-site blocks; entries multiply in their auxiliary order."""
    return sp.kronecker_product(x,sp.eye(2)), sp.kronecker_product(sp.eye(2),x)


def lower(q: sp.MatrixBase, left=RHO, right=RHO):
    return sp.expand(sp.trace(sp.kronecker_product(left,right)*q))


class Checks:
    def __init__(self):
        self.results = {}

    def zero(self, name, value):
        if isinstance(value,sp.MatrixBase):
            result = value.applyfunc(sp.expand)
            ok = result == sp.zeros(*result.shape)
        else:
            result = sp.expand(value)
            ok = result == 0
        if not ok:
            raise AssertionError(f'{name}: {result}')
        self.results[name] = 'exact_zero'
        print(name+': exact_zero', flush=True)

    def nonzero(self,name,value):
        result = sp.expand(value)
        if result == 0:
            raise AssertionError(name+': expected a nonzero control')
        self.results[name] = 'expected_nonzero'
        print(name+': expected_nonzero',flush=True)


def verify(checks: Checks):
    x00,x01,x10,x11 = sp.symbols('X00 X01 X10 X11',commutative=False)
    x = sp.Matrix([[x00,x01],[x10,x11]])
    y = sp.Matrix(2,2,sp.symbols('Y00 Y01 Y10 Y11',commutative=False))
    z = sp.Matrix(2,2,sp.symbols('Z00 Z01 Z10 Z11',commutative=False))
    x1,x2 = physical_embeddings(x)
    y1,y2 = physical_embeddings(y)
    z1,z2 = physical_embeddings(z)
    h = sp.Matrix(4,4,sp.symbols('h0:16'))
    p,q,c,d = sp.symbols('p q c d')
    drho = sp.Matrix([[0,p],[q,0]])
    plus = (sp.eye(4)+P)/2
    allowed = (h-P*h*P)/2+c*sp.eye(4)+d*(sp.eye(4)-P)
    condition = {h[i,j]:allowed[i,j] for i in range(4) for j in range(4)}
    # Simultaneous substitution is important: the general parameter symbols
    # also occur in allowed. Sequential substitution would change the family.
    impose = lambda e: e.subs(condition, simultaneous=True)
    u=x00
    xi=x01*x10
    ux = sp.trace(drho*x)
    xi_x = sp.expand(sp.trace(drho*x*x)-ux*u-u*ux)
    a1=-I*comm(2*P,x2)
    a2_exchange=I*x2*comm(2*P,x2)
    a2_y=-I*comm(2*P,y2)
    a2_h=-I*comm(h,x2)
    a2_z=-I*z2
    a2=a2_exchange+a2_y+a2_h+a2_z
    sutherland=(comm(2*P,x2*x1+y2+y1)+comm(h,x1+x2)-z1+z2)
    kh=lower(comm(h,x2))
    checks.zero('permutation.intertwining', P*x1-x2*P)
    checks.zero('projector.tangent',RHO*drho+drho*RHO-drho)
    checks.zero('projector.second_moment',lower(x2*x2)-u*u-xi)
    checks.zero('triplet.general_family',plus*allowed*plus-c*plus)
    checks.zero('time.first_coincident_symbol',lower(a1))
    checks.zero('quadratic.plus',lower(a1*x1)-2*I*xi)
    checks.zero('quadratic.minus',lower(x2*a1)-2*I*xi)
    plus_taylor=lower(a1*x1,RHO,drho)-lower(a1,RHO,drho)*u
    minus_taylor=lower(x2*a1,drho,RHO)-u*lower(a1,drho,RHO)
    checks.zero('quadratic.right_tangent',plus_taylor-2*I*p*(x11-x00)*x10)
    checks.zero('quadratic.left_tangent',minus_taylor-2*I*q*x01*(x11-x00))
    checks.zero('quadratic.total_tangent',plus_taylor+minus_taylor-2*I*xi_x)
    checks.zero('time.second_coincident_symbol',lower(a2)+2*I*xi+I*kh+I*z[0,0])
    cubic=lambda a: lower(a*x1-x2*a)-comm(lower(a),u)
    bx=comm(2*P,x2*x1)
    by=comm(2*P,y1+y2)
    bh=comm(h,x1+x2)
    bz=-z1+z2
    checks.zero('cubic.exchange',cubic(a2_exchange)-2*I*comm(xi,u)+I*lower(x1*bx))
    checks.zero('cubic.second_spatial_cancellation',cubic(a2_y)+lower(a1*y1-y2*a1))
    checks.zero('Sutherland.second_spatial_cancellation',by)
    checks.zero('cubic.spectral',cubic(a2_z)+I*lower(x1*bz))
    # Before imposing the triplet restriction, retain the complete defect.
    h_defect=I*(lower((x1+x2)*comm(h,x1+x2))-lower(comm(h,x2*x1))+comm(kh,u))
    checks.zero('cubic.arbitrary_h_defect',cubic(a2_h)+I*lower(x1*bh)-h_defect)
    checks.zero('triplet.collective_contraction',impose(lower((x1+x2)*comm(h,x1+x2))))
    checks.zero('triplet.product_commutator',impose(lower(comm(h,x2*x1))-comm(kh,u)))
    checks.zero('triplet.defect_vanishes',impose(h_defect))
    cubic_full=cubic(a2)+lower(a1*y1-y2*a1)
    checks.zero('source.arbitrary_h_identity',plus_taylor+minus_taylor+cubic_full
                -2*I*(xi_x+comm(xi,u))+I*lower(x1*sutherland)-h_defect)
    checks.zero('source.conditional_identity',impose(plus_taylor+minus_taylor+cubic_full
                -2*I*(xi_x+comm(xi,u))+I*lower(x1*sutherland)))
    h_even=(h+P*h*P)/2
    checks.zero('rank.exchange_parity',(bx+P*bx*P)/2)
    checks.zero('rank.spectral_parity',(bz+P*bz*P)/2)
    checks.zero('rank.even_Sutherland',(sutherland+P*sutherland*P)/2-comm(h_even,x1+x2))
    # Three independent auxiliary coefficients of the physical Pauli matrices
    # make even_Sutherland=0 imply [h_even, sigma1^i+sigma2^i]=0.
    hvars=list(h)
    collective=[sum(physical_embeddings(sigma),sp.zeros(4)) for sigma in PAULI]
    equations=[entry for generator in collective for entry in comm(h,generator)]
    coefficient,_=sp.linear_eq_to_matrix(equations,hvars)
    checks.zero('rank.collective_commutant_dimension',len(hvars)-coefficient.rank()-2)
    for i,generator in enumerate(collective,1):
        checks.zero(f'rank.collective_commutant_basis_{i}',comm(c*sp.eye(4)+d*P,generator))
    # Exact ranks of the spin-to-auxiliary maps in the three model examples.
    lam,alpha,a,b=sp.symbols('lambda alpha a b', nonzero=True)
    perm=P
    k5=sp.Matrix([[0,1,-1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    k6=sp.Matrix([[0,1,1,0],[0,0,0,-1],[0,0,0,-1],[0,0,0,0]])
    models={
        'class5':perm/(2*lam)+alpha*k5,
        'class6':perm/(2*lam)+alpha*lam*k6+alpha**2*lam**3*k6*k6,
        'XXZ':sp.Matrix([[a,0,0,0],[0,0,b,0],[0,b,0,0],[0,0,0,a]])}
    rank_minors={}
    for name,full in models.items():
        # These model matrices have the conventional auxiliary-first ordering.
        spin_coeffs=[sp.Matrix(2,2,lambda r,c:sp.trace(sigma*full[2*r:2*r+2,2*c:2*c+2])/2)
                     for sigma in PAULI]
        # Rows U12,U21,U11 give a nonzero 3x3 minor.
        minor=sp.Matrix([[m[r,c] for m in spin_coeffs] for r,c in [(0,1),(1,0),(0,0)]]).det()
        expected=I*a*b*b/4 if name=='XXZ' else I/(32*lam**3)
        checks.zero(name+'.rank_three_minor',minor-expected)
        rank_minors[name]=str(sp.factor(minor))
    # A rank-deficient example shows that rank must not be silently dropped.
    x_rank_one=sp.Matrix([[x00,0],[0,-x00]])
    r1,r2=physical_embeddings(x_rank_one)
    h_bad=sp.kronecker_product(PAULI[2],PAULI[2])
    b_rank_one=comm(2*P,r2*r1)+comm(h_bad,r1+r2)
    checks.zero('rank_deficient.Sutherland',b_rank_one)
    nonscalar=plus*h_bad*plus-plus
    if nonscalar==sp.zeros(4):
        raise AssertionError('rank-deficient control must violate scalar triplet restriction')
    checks.results['rank_deficient.triplet_not_scalar']='expected_nonzero'
    # This is NOT a counterexample to the correction: only to deriving the
    # triplet condition from Sutherland without rank-three spin dependence.
    # An arbitrary first-order h need not have vanishing defect either.
    bad_sub={h[i,j]:h_bad[i,j] for i in range(4) for j in range(4)}
    bad_defect=sp.expand(h_defect.subs(bad_sub, simultaneous=True))
    checks.nonzero('arbitrary_h.defect_can_be_nonzero',bad_defect)
    return {
        'rank_minors':rank_minors,
        'generic_h_defect':str(h_defect),
        'nonzero_h_defect_example':str(bad_defect),
        'definition_h_defect':'i*(<(X1+X2)[h1,X1+X2]>-<[h1,X2 X1]>+[<[h1,X2]>,U])',
        'scope':'Physical dimension two; arbitrary finite auxiliary dimension; rank-one product lower symbols; h0=2P.',
        'proof':'Free noncommuting auxiliary blocks; exact polynomial cancellation. Scalar triplet restriction follows from order-two Sutherland when the three auxiliary Pauli coefficients of X are linearly independent.',
        'not_claimed':['literature priority','higher physical spins','higher Hamiltonian flows','finite-time quantum convergence']
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    default=Path(__file__).resolve().parents[2]/'results/generalization/shared_site_proof.json'
    parser.add_argument('--out',type=Path,default=default)
    args=parser.parse_args()
    checks=Checks()
    content=verify(checks)
    report={'schema_version':1,'base_commit':BASE_COMMIT,'python_version':platform.python_version(),
            'sympy_version':sp.__version__,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'number_of_checks':len(checks.results),'checks':checks.results,**content}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f'All {len(checks.results)} checks passed. Report: {args.out}')


if __name__=='__main__':
    main()
