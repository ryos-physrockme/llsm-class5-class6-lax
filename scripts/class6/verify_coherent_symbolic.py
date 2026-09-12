"""Exact Class-6 quantum-to-classical Lax calculation (spin 1/2).

Input: R(u,a) of arXiv:2506.13598v2, Eq. (A.1), with a1=1.
The connected source is computed from operator products BEFORE a classical
V is constructed. No saved matrices, previous Class-6 programs or fitted
coefficients are imported. The tensor-contraction method extends the
Class-5 report; epsilon**2 terms of the Class-6 spatial operator are retained.

For short code expressions only, p=S^1+i*S^2, q=S^1-i*S^2, z=S^3.
They are spin components, not momenta or spectral parameters.
Run: python code/verify_symbolic.py
"""
from pathlib import Path
import json
import sys
import sympy as s

OUT = Path(__file__).resolve().parents[2] / 'results' / 'class6'
OUT.mkdir(exist_ok=True)
I = s.I
P = s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
K = s.Matrix([[0,1,1,0],[0,0,0,-1],[0,0,0,-1],[0,0,0,0]])
la = s.Symbol('lambda', nonzero=True)
al = s.Symbol('alpha')
u, aq, eps = s.symbols('u a epsilon')
checks = {}

def embed(X: s.Matrix, first: int, second: int) -> s.Matrix:
    """Ordered two-qubit operator, tensor order (auxiliary, left, right)."""
    if first == second or any(t not in (0,1,2) for t in (first,second)):
        raise ValueError('Two distinct tensor positions are required.')
    Y = s.zeros(8)
    for row in range(8):
        r = [(row >> k) & 1 for k in (2,1,0)]
        for col in range(8):
            c = [(col >> k) & 1 for k in (2,1,0)]
            if all(r[k] == c[k] for k in range(3) if k not in (first,second)):
                Y[row,col] = X[2*r[first]+r[second],2*c[first]+c[second]]
    return Y

def comm(A, B):
    return A*B-B*A

def raw_check(name, M):
    value = M.applyfunc(lambda x: s.factor(s.expand(x)))
    if value != s.zeros(*value.shape):
        raise AssertionError((name, value))
    checks[name] = 'exact_zero'
    print(name + ': exact_zero', flush=True)

# All quantum coefficients are checked against the full R, not assumed.
R = 2*u*s.eye(4)+P+aq*u*(1+2*u)*K+aq**2*u**2*(1+2*u)**2*K*K/2
h = 2*P+aq*K
ell = P/(2*la)+al*la*K+al**2*la**3*K*K
m = al*K/2+al**2*la**2*K*K
n = al**2*la*K*K/4
raw_check('regularity_R_zero_is_P', R.subs(u,0)-P)
raw_check('Hamiltonian_from_R_derivative', P*R.diff(u).subs(u,0)-h)
raw_check('P_squared', P*P-s.eye(4))
raw_check('K_cubed', K**3)
raw_check('P_K_equals_K', P*K-K)
raw_check('K_P_equals_K', K*P-K)
raw_check('inverse_relation_R_u_R_minus_u', R*R.subs(u,-u)-(1-4*u**2)*s.eye(4))
scaled_R = (R/(2*u)).subs({u:la/eps,aq:eps**2*al})
raw_check('exact_spatial_epsilon_expansion', scaled_R-s.eye(4)-eps*ell-eps**2*m-eps**3*n)
raw_check('ell_squared_contains_second_order', ell**2-s.eye(4)/(4*la**2)-2*m)
raw_check('fixed_coupling_spectral_derivative',
    (R/(2*u)).diff(u).subs({u:la/eps,aq:eps**2*al})
    -eps**2*ell.diff(la)-eps**3*m.diff(la)-eps**4*n.diff(la))
RL, RR = embed(R,0,1), embed(R,0,2)
hbond = embed(h,1,2)
raw_check('Sutherland_relation', comm(hbond,RR*RL)-RR*RL.diff(u)+RR.diff(u)*RL)

# Expand the quantum time component to the order required for the source.
ll, lr = embed(ell,0,1), embed(ell,0,2)
ml, mr = embed(m,0,1), embed(m,0,2)
h0 = 2*embed(P,1,2)
A1 = -I*comm(h0,lr)
A2 = I*lr*comm(h0,lr)-I*comm(h0,mr)-I*lr.diff(la)
# Check the expansion by multiplying back; no matrix inverse is needed here.
raw_check('time_operator_expansion_to_second_order',
    A2+lr*A1+I*comm(h0,mr)+I*lr.diff(la))

left = s.Matrix(s.symbols('left_plus left_minus left_z'))
right = s.Matrix(s.symbols('right_plus right_minus right_z'))
S = s.Matrix(s.symbols('p q z'))
X = s.Matrix(s.symbols('p_x q_x z_x'))
Y = s.Matrix(s.symbols('p_xx q_xx z_xx'))

def rho(v):
    return s.Matrix([[1+v[2],v[1]],[v[0],1-v[2]]])/2

weight = s.kronecker_product(rho(left),rho(right))
def low(Q):
    """Physical partial trace, leaving the auxiliary matrix order intact."""
    Q = Q.applyfunc(s.expand)
    answer = s.zeros(2)
    for i in range(2):
        for j in range(2):
            answer[i,j] = s.expand(s.trace(weight*Q[4*i:4*i+4,4*j:4*j+4]))
    return answer

def d(M, variables, direction):
    return sum((M.diff(variables[k])*direction[k] for k in range(3)), s.zeros(*M.shape))

def dx(M):
    return d(M,S,X)+d(M,X,Y)

same = dict(zip(list(left)+list(right),list(S)+list(S)))
def at(M):
    return M.subs(same, simultaneous=True).applyfunc(s.expand)

f1, f2 = low(A1), low(A2)
Ul, Ur, Ml, Mr = low(ll), low(lr), low(ml), low(mr)
R2, L2 = low(A1*ll), low(lr*A1)
# Class 6 requires the products involving the second spatial coefficient m.
R3, L3 = low(A2*ll+A1*ml), low(lr*A2+mr*A1)
CR2, CL2 = R2-f1*Ul, L2-Ur*f1
CR3, CL3 = R3-f2*Ul-f1*Ml, L3-Ur*f2-Mr*f1
source = at(d(CR2,right,X)+d(CL2,left,X)+CR3-CL3)
Vdown = at(f2-d(f1,left,X))
# Direct quantum RHS, including both the ordinary and connected products.
quantum_flow = at((d(f1,right,Y)-d(f1,left,Y)
    +d(d(f1,right,X),right,X)-d(d(f1,left,X),left,X))/2
    +d(f2,right,X)+d(f2,left,X)+d(R2,right,X)+d(L2,left,X)+R3-L3)
U, M = at(Ur), at(Mr)
raw_check('permutation_lower_symbol', at(low(embed(P,0,2)))-rho(S))
raw_check('deformation_lower_symbol', at(low(embed(K,0,2)))
    -s.Matrix([[S[0]/2,S[2]],[0,-S[0]/2]]))
raw_check('deformation_square_lower_symbol', at(low(embed(K*K,0,2)))
    -s.Matrix([[0,-S[0]],[0,0]]))
bond_symbol = s.trace(s.kronecker_product(rho(left),rho(right))*h)
bond_formula = (1+(left[0]*right[1]+left[1]*right[0])/2
    +left[2]*right[2]+aq*(left[0]*right[2]+left[2]*right[0])/2)
raw_check('two_site_Hamiltonian_lower_symbol',s.Matrix([bond_symbol-bond_formula]))
raw_check('quantum_rhs_decomposition', quantum_flow-dx(Vdown)-comm(Vdown,U)-source)
raw_check('coincident_A1', at(f1))

# Only kinematic unit-spin constraints are used, not equations of motion.
p,q,z = S
px,qx,zx = X
pxx,qxx,zxx = Y
constraints = [p*q+z**2-1,
               q*px+p*qx+2*z*zx,
               q*pxx+p*qxx+2*px*qx+2*z*zxx+2*zx**2]
variables = list(Y)+list(X)+list(S)
G = s.groebner(constraints,*variables,domain=s.QQ_I.frac_field(la,al))
def reduce(M):
    return M.applyfunc(lambda e:s.factor(G.reduce(s.expand(e))[1]))

def check(name, M):
    value = reduce(M)
    if value != s.zeros(*value.shape):
        raise AssertionError((name,value))
    checks[name] = 'exact_zero'
    print(name+': exact_zero',flush=True)

check('connected_order_epsilon_squared_cancels',at(CR2-CL2))

def cart(v):
    return s.Matrix([(v[0]+v[1])/2,(v[0]-v[1])/(2*I),v[2]])

def components(v):
    return s.Matrix([v[0]+I*v[1],v[0]-I*v[1],v[2]])

def Q(v):
    """v uses (S-plus,S-minus,S-z) components, not Cartesian components."""
    return s.Matrix([[v[2]+2*al*la**2*v[0],
                     v[1]+4*al*la**2*v[2]-4*al**2*la**4*v[0]],
                     [v[0],-v[2]-2*al*la**2*v[0]]])/(4*la)

raw_check('spatial_matrix_formula', U-s.eye(2)/(4*la)-Q(S))
raw_check('A1_two_site_formula', f1+2*Q(components(cart(left).cross(cart(right)))))
W = -2*Q(components(cart(S).cross(cart(X))))
U0 = U-s.eye(2)/(4*la)
# Intermediate formulas and component entries printed in manuscript section 4.2
# are checked against the independently contracted quantum coefficients.
raw_check('second_spatial_lower_symbol_formula',
    M-s.Matrix([[al*p/4,al*z/2-al**2*la**2*p],[0,-al*p/4]]))
check('coincident_second_time_coefficient',
    at(f2)+2*I*U0.diff(la)-I*al*p*z*s.eye(2))
check('spatial_matrix_determinant',s.Matrix([U.det()+al*p*z/2]))
raw_check('spectral_derivative_from_second_spatial_symbol',
    U0.diff(la)-4*M+U0/la)
Xi = at(low(lr*lr))-U*U
# Check the derivative and algebraic pieces separately, before combining them.
check('right_quadratic_overlap_at_equal_spins',at(CR2)-2*I*Xi)
check('left_quadratic_overlap_at_equal_spins',at(CL2)-2*I*Xi)
check('quadratic_overlap_Taylor_contribution',
      at(d(CR2,right,X)+d(CL2,left,X))-2*I*dx(Xi))
check('coincident_cubic_overlap_commutator',at(CR3-CL3)-2*I*comm(Xi,U))
# The following identity is verified against the independently contracted source.
check('source_equals_covariant_second_moment', source-2*I*(dx(Xi)+comm(Xi,U)))
# Construct Delta V only AFTER the previous identity has been verified.
D = 2*I*Xi-I*s.eye(2)/(4*la**2)
V = Vdown+D
check('source_absorbed_by_local_correction', source-dx(D)-comm(D,U))
check('local_correction_from_m', D-(4*I*M-2*I*U*U+I*s.eye(2)/(4*la**2)))
check('local_correction_spectral_formula', D-I*U0.diff(la)+I*al*p*z*s.eye(2))
check('raw_time_symbol_spectral_formula', Vdown-W+2*I*U0.diff(la)-I*al*p*z*s.eye(2))
check('corrected_time_component_spectral_formula', V-W+I*U0.diff(la))
check('corrected_time_component_trace_zero',s.Matrix([s.trace(V)]))
v11 = I*(z-2*al*la**2*p+la*(q*px-p*qx)
         +4*al*la**3*(p*zx-z*px))/(4*la**2)
v12 = I*(q-4*al*la**2*z+12*al**2*la**4*p
         +2*la*(z*qx-q*zx)+4*al*la**3*(q*px-p*qx)
         +8*al**2*la**5*(z*px-p*zx))/(4*la**2)
v21 = I*(p+2*la*(p*zx-z*px))/(4*la**2)
check('corrected_time_component_explicit_entries',
    V-s.Matrix([[v11,v12],[v21,-v11]]))
check('source_explicit_spatial_and_commutator_terms',
      source+2*I*dx(U*U)-4*I*dx(M)-4*I*comm(M,U))

# Independent Hamiltonian variation and equation-of-motion comparison.
H = -(px*qx+zx*zx)/2+al*p*z
variation = s.Matrix([s.diff(H,S[k])-dx(s.Matrix([s.diff(H,X[k])]))[0] for k in range(3)])
Pi = s.Matrix([[0,-4*I*z,2*I*p],[4*I*z,0,-2*I*q],[-2*I*p,2*I*q,0]])
flow = s.Matrix([2*I*(p*zxx-z*pxx+al*p*p),
                 2*I*(z*qxx-q*zxx+2*al*z*z-al*p*q),
                 I*(q*pxx-p*qxx-2*al*p*z)])
raw_check('Hamiltonian_variation',Pi*variation-flow)
check('unit_spin_is_preserved',s.Matrix([q*flow[0]+p*flow[1]+2*z*flow[2]]))
check('quantum_rhs_equals_Hamiltonian_flow',quantum_flow-Q(flow))
check('corrected_zero_curvature_on_flow',Q(flow)-dx(V)+comm(U,V))
T = s.Matrix(s.symbols('p_T q_T z_T'))
F = Q(T)-dx(V)+comm(U,V)
check('off_shell_factorization',F-Q(T-flow))
# Inverse map on arbitrary residual components.
e = s.Matrix(s.symbols('E_plus E_minus E_z'))
Ftest=Q(e)
recovered = s.Matrix([4*la*Ftest[1,0],
    4*la*(Ftest[0,1]-4*al*la**2*Ftest[0,0]+12*al**2*la**4*Ftest[1,0]),
    4*la*(Ftest[0,0]-2*al*la**2*Ftest[1,0])])
raw_check('inverse_EOM_map',recovered-e)
# Matrix of the three independent components (21,11,12), determinant nonzero.
qcoeff=s.Matrix([Ftest[1,0],Ftest[0,0],Ftest[0,1]]).jacobian(e)
raw_check('EOM_map_determinant',s.Matrix([qcoeff.det()+1/(64*la**3)]))

# Direct Hamiltonian symbol expansion.
state = s.kronecker_product(rho(S),rho(S+eps*X+eps**2*Y/2))
hsymbol=s.expand(s.trace(state*h.subs(aq,eps**2*al)))
check('energy_density_order_zero',s.Matrix([hsymbol.coeff(eps,0)-2]))
check('energy_density_order_one',s.Matrix([hsymbol.coeff(eps,1)]))
check('energy_density_order_two',s.Matrix([hsymbol.coeff(eps,2)-H]))
N=s.Matrix([[0,0,-1],[0,0,-I],[-1,-I,0]])
raw_check('anisotropy_nilpotent',N**3)
raw_check('anisotropy_Hamiltonian_density',
    s.Matrix([H+cart(X).dot(cart(X))/2+al*(cart(S).T*N*cart(S))[0]/2]))
raw_check('vector_Hamiltonian_flow',cart(flow)+2*cart(S).cross(cart(Y)-al*N*cart(S)))

# The instantaneous spin derivative printed in manuscript section 4.3 is
# checked directly from the two adjacent quantum bonds, before taking a limit.
pauli = [s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-I],[I,0]]), s.diag(1,-1)]
spin_hamiltonian = embed(h,0,1)+embed(h,1,2)
spin_state = s.kronecker_product(rho(left),rho(S),rho(right))
quantum_spin_flow = s.Matrix([
    s.trace(spin_state*I*comm(spin_hamiltonian,
        s.kronecker_product(s.eye(2),sigma,s.eye(2)))) for sigma in pauli
]).applyfunc(s.expand)
neighbors = cart(left)+cart(right)
raw_check('finite_lattice_spin_Heisenberg_symbol',
    quantum_spin_flow+2*cart(S).cross(neighbors-aq*N*neighbors/2))
spin_taylor = dict(zip(left,S-eps*X+eps**2*Y/2))
spin_taylor.update(zip(right,S+eps*X+eps**2*Y/2))
spin_taylor[aq] = eps**2*al
continuum_spin_flow = quantum_spin_flow.subs(spin_taylor,simultaneous=True)
raw_check('finite_lattice_spin_continuum_limit',
    continuum_spin_flow.applyfunc(lambda value:s.expand(value).coeff(eps,2))
    -cart(flow))

reflection=s.diag(-1,1,1)
Nref=s.Matrix([[0,0,1],[0,0,-I],[1,-I,0]])
Sref,Yref=reflection*cart(S),reflection*cart(Y)
raw_check('source_EOM_dictionary',reflection*cart(flow)-2*Sref.cross(Yref-al*Nref*Sref))

# XXX limit and a negative control: a spatially uniform spin is not stationary.
check('XXX_time_component',V.subs(al,0)-W.subs(al,0)-I*(2*rho(S)-s.eye(2))/(4*la**2))
isotropic_spin = sum((cart(S)[k]*pauli[k] for k in range(3)),s.zeros(2))
isotropic_spin_x = sum((cart(X)[k]*pauli[k] for k in range(3)),s.zeros(2))
check('XXX_local_correction',D.subs(al,0)+I*isotropic_spin/(4*la**2))
check('XXX_connected_source',source.subs(al,0)+I*isotropic_spin_x/(4*la**2))
uniform=dict(zip(list(S)+list(X)+list(Y),[0,0,1]+[0]*6))
uniform_source=s.Matrix([[0,-I*al/la],[0,0]])
raw_check('uniform_spin_nonzero_connected_source',source.subs(uniform)-uniform_source)
D_class5=-2*I*U*U+I*s.eye(2)/(4*la**2)
wrong_residual=reduce(source-dx(D_class5)-comm(D_class5,U))
if wrong_residual.subs(uniform) == s.zeros(2):
    raise AssertionError('The Class-5-only negative control unexpectedly passed.')
negative={'Class5_only_correction_residual_at_uniform_spin':str(wrong_residual.subs(uniform)),
          'result':'nonzero_for_alpha_nonzero_and_lambda_nonzero'}

summary={'python_version':sys.version.split()[0], 'sympy_version':s.__version__,
         'number_of_exact_checks':len(checks),'checks':checks,'negative_control':negative,
         'scope':'spin 1/2; Class 6 R of 2506.13598v2 (A.1), a1=1; a=epsilon^2 alpha, u=lambda/epsilon, T=epsilon^2 t_lat',
         'constraint_polynomials':[str(c) for c in constraints],
         'groebner_variable_order':[str(v) for v in variables],
         'groebner_coefficient_field':'QQ_I(lambda,alpha)'}
(OUT/'symbolic_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
expressions={name:str(reduce(M)) for name,M in {
    'U':U,'U_traceless':U0,'second_spatial_symbol':M,'second_moment_difference':Xi,
    'A1_symbol':f1,'A2_symbol':f2,'time_symbol':Vdown,'connected_source':source,
    'local_correction':D,'corrected_time_component':V,'flow_plus_minus_z':flow,
    'negative_control_residual':wrong_residual}.items() if name not in ('A1_symbol','A2_symbol')}
expressions.update({'A1_symbol':str(f1),'A2_symbol':str(f2)})
(OUT/'symbolic_expressions.json').write_text(json.dumps(expressions,indent=2)+'\n')
print('All %d exact checks passed; negative control rejected as expected.' % len(checks))
