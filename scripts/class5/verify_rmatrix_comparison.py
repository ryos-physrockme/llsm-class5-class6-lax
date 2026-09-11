import sympy as sp

I = sp.I
p,q,z,px,qx,zx,alpha,lam,mu = sp.symbols('p q z p_x q_x z_x alpha lambda mu', nonzero=True)

# Physical definitions
rho = sp.Matrix([[1+z, q],[p, 1-z]])/2
rho_x = sp.Matrix([[zx, qx],[px, -zx]])/2
Pup = sp.Matrix([[1,0],[0,0]])
Eplus = sp.Matrix([[0,1],[0,0]])
K = sp.Matrix([[p, -(1+z)],[0,0]])/2
Id2 = sp.eye(2)
Pperm = sp.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
Cj = sp.Matrix([[0,1,-1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])


def U(spec):
    return rho/(2*spec) + alpha*K

U_lam = U(lam)
U_mu = U(mu)

# Classical Poisson brackets in the normalization of the previous Class-5 report.
vars3 = [p,q,z]
pb_table = {
    (p,q):-4*I*z, (q,p):4*I*z,
    (p,z):2*I*p, (z,p):-2*I*p,
    (q,z):-2*I*q, (z,q):2*I*q,
}
def PB(f,g):
    out = 0
    for a in vars3:
        for b in vars3:
            out += sp.diff(f,a)*sp.diff(g,b)*pb_table.get((a,b),0)
    return sp.expand(out)

# Continuum classical r-matrix
r12 = I*Pperm/(2*(lam-mu)) + I*alpha*Cj
U1pU2 = sp.kronecker_product(U_lam,Id2) + sp.kronecker_product(Id2,U_mu)
r_rhs = sp.expand(r12*U1pU2-U1pU2*r12)
r_lhs = sp.zeros(4)
for a in range(2):
    for b in range(2):
        for c in range(2):
            for d in range(2):
                r_lhs[2*a+c,2*b+d] = PB(U_lam[a,b],U_mu[c,d])

# Formal monodromy eigenprojector Pi = rho + lambda Pi1 + ...
def comm(A,B):
    return A*B-B*A
Pi1 = sp.expand(2*comm(rho,rho_x)-2*alpha*comm(rho,comm(K,rho)))

# Groebner reduction by spin-length and tangent constraints.
G = sp.groebner(
    [p*q+z**2-1, px*q+p*qx+2*z*zx],
    px,qx,zx,p,q,z,alpha,lam,mu,
    order='lex', domain=sp.QQ_I
)
def reduce_constraints(expr):
    expr = sp.together(sp.expand(expr))
    num, den = sp.fraction(expr)
    rem = G.reduce(sp.expand(num))[1]
    return sp.factor(rem/den)

def mat_reduce(M):
    return M.applyfunc(reduce_constraints)

# Projector recursion checks
proj_linear = sp.expand(rho*Pi1 + Pi1*rho - Pi1)
transport_order0 = sp.expand(rho_x - (sp.Rational(1,2)*comm(rho,Pi1) + alpha*comm(K,rho)))

# STS generating time matrix: coefficient of lambda^1 in the local hierarchy.
# tr_a(Pi_a P_ab)=Pi_b and Cj=Pup\otimes Eplus-Eplus\otimes Pup.
tr_Pi1_Pup = sp.trace(Pi1*Pup)
tr_Pi1_Eplus = sp.trace(Pi1*Eplus)
Vgen1 = -I/(2*mu)*(Pi1 + rho/mu) + I*alpha*(tr_Pi1_Pup*Eplus - tr_Pi1_Eplus*Pup)
# The previous report uses U_T={U,H}, whereas the STS formula is written as {H,U}=... .
V_rmatrix = -Vgen1

# Quantum-finite-lattice continuum time matrix from the previous Class-5 report.
W = sp.Matrix([
    [-I/(4*mu)*(p*qx-px*q+4*alpha*mu*(z*px-p*zx)),
      I/(2*mu)*(z*qx-q*zx+alpha*mu*(p*qx-px*q))],
    [ I/(2*mu)*(p*zx-z*px),
      I/(4*mu)*(p*qx-px*q)]
])
Umu = U_mu
Z = 2*I*(Umu*Umu) - I*Id2/(4*mu**2)
V_quantum_limit = sp.expand(W+Z)
V_difference = sp.expand(V_rmatrix - V_quantum_limit)
expected_scalar = I*Id2/(4*mu**2)

checks = {
    'linear r-matrix Poisson algebra': r_lhs-r_rhs,
    'projector constraint at O(lambda)': proj_linear,
    'projector transport at O(lambda^0)': transport_order0,
    'STS vs quantum-limit V up to scalar': V_difference-expected_scalar,
}

all_ok = True
for name, M in checks.items():
    Mr = mat_reduce(M)
    ok = all(x == 0 for x in Mr)
    all_ok &= ok
    print(f'{name}:', 'PASS' if ok else 'FAIL')
    if not ok:
        print(Mr)

print('\nCentral result after constraints:')
print('V_rmatrix - V_quantum_limit =')
sp.pprint(mat_reduce(V_difference))
print('\nALL CHECKS:', 'PASS' if all_ok else 'FAIL')
if not all_ok:
    raise SystemExit(1)
