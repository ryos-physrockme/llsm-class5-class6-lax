"""Class-5 quantum-to-classical time Lax correction: exact independent derivation.

Starts from the fundamental R matrix, not from an input continuum V.
All asserted continuum equalities are polynomial identities modulo the unit-spin
constraint and its first two x derivatives. No numerical fitted coefficients,
previous project helper functions, cached pickles, or network access are used.
"""
import sympy as s
import json, sys
from pathlib import Path
out=Path(__file__).resolve().parents[2]/'results'/'class5'
out.mkdir(exist_ok=True)
checks={}
I=s.I; Id=s.eye(2)
sig=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-I],[I,0]]),s.diag(1,-1)]
P=s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
K=s.Matrix([[0,1,-1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
def embed(X,a,b):
    Y=s.zeros(8)
    for row in range(8):
      r=[(row>>k)&1 for k in (2,1,0)]
      for col in range(8):
        c=[(col>>k)&1 for k in (2,1,0)]
        if all(r[k]==c[k] for k in range(3) if k not in (a,b)):
          Y[row,col]=X[2*r[a]+r[b],2*c[a]+c[b]]
    return Y
la=s.Symbol('lambda', nonzero=True)
al=s.Symbol('alpha')
ll=embed(P,0,1)/(2*la)+al*embed(K,0,1)
lr=embed(P,0,2)/(2*la)+al*embed(K,0,2)
h0=2*embed(P,1,2); h1=2*al*embed(K,1,2)
comm=lambda A,B:A*B-B*A
A1=-I*comm(h0,lr)
A2=I*lr*comm(h0,lr)-I*comm(h1,lr)+I*embed(P,0,2)/(2*la**2)
a=s.Matrix(s.symbols('a1:4')); b=s.Matrix(s.symbols('b1:4'))
S=s.Matrix(s.symbols('p q z')); X=s.Matrix(s.symbols('p_x q_x z_x')); Y=s.Matrix(s.symbols('p_xx q_xx z_xx'))
def rho(v): return s.Matrix([[1+v[2],v[1]],[v[0],1-v[2]]])/2
def low(Q):
    Q=Q.applyfunc(s.expand)
    rr=s.kronecker_product(rho(a),rho(b)); ret=s.zeros(2)
    for i in range(2):
      for j in range(2): ret[i,j]=s.expand(s.trace(rr*Q[4*i:4*i+4,4*j:4*j+4]))
    return ret
print('operators built',flush=True)
f1=low(A1); f2=low(A2)
print('A symbols done',flush=True)
Lleft=low(ll); Lright=low(lr)
R2=low(A1*ll); L2=low(lr*A1)
R3=low(A2*ll); L3=low(lr*A2)
print('product symbols done',flush=True)
CR2=(R2-f1*Lleft).applyfunc(s.expand); CL2=(L2-Lright*f1).applyfunc(s.expand)
CR3=(R3-f2*Lleft).applyfunc(s.expand); CL3=(L3-Lright*f2).applyfunc(s.expand)
subs=dict(zip(list(a)+list(b),list(S)+list(S)))
def at(M): return M.subs(subs, simultaneous=True).applyfunc(s.expand)
def d(M,vs,v): return sum([M.diff(vs[k])*v[k] for k in range(3)],s.zeros(2))
Vdown=at(f2-d(f1,a,X))
C2=at(CR2-CL2)
C=at(d(CR2,b,X)+d(CL2,a,X)+CR3-CL3)
U=rho(S)/(2*la)+al*s.Matrix([[S[0]/2,-(1+S[2])/2],[0,0]])
# Obtain exact continuum evolution of U directly from expansion of the quantum RHS.
E3=at((d(f1,b,Y)-d(f1,a,Y)+d(d(f1,b,X),b,X)-d(d(f1,a,X),a,X))/2+d(f2,b,X)+d(f2,a,X)+d(R2,b,X)+d(L2,a,X)+R3-L3)
Fnaive=d(Vdown,S,X)+d(Vdown,X,Y)+comm(Vdown,U)
res=(E3-Fnaive-C).applyfunc(s.expand)
print('exact decomposition:',res,flush=True)
# polynomial constraints with a specified lexicographic leading order
polys=[S[0]*S[1]+S[2]**2-1,S[1]*X[0]+S[0]*X[1]+2*S[2]*X[2],S[1]*Y[0]+S[0]*Y[1]+2*X[0]*X[1]+2*S[2]*Y[2]+2*X[2]**2]
vars=list(Y)+list(X)+list(S)
G=s.groebner(polys,*vars,domain=s.QQ_I.frac_field(la,al))
def red(e): return s.factor(G.reduce(s.expand(e))[1])
def reduce(M): return M.applyfunc(red)

def check(name, M, constrained=True):
    value=reduce(M) if constrained else M.applyfunc(s.simplify)
    assert value == s.zeros(*value.shape), (name,value)
    checks[name] = "exact_zero"
    print(name+": exact_zero",flush=True)

# These algebraic identities use the quantum operators before taking symbols.
check('P_squared',P*P-s.eye(4),False)
check('K_squared',K*K,False)
check('PK_minus_K',P*K-K,False)
check('KP_plus_K',K*P+K,False)
check('ell_squared',lr*lr-s.eye(8)/(4*la**2),False)
u,a_q=s.symbols('u a_q')
R_L=2*u*s.eye(8)+embed(P,0,1)+2*a_q*u*embed(K,0,1)
R_R=2*u*s.eye(8)+embed(P,0,2)+2*a_q*u*embed(K,0,2)
h=2*embed(P,1,2)+2*a_q*embed(K,1,2)
check('Sutherland_relation',comm(h,R_R*R_L)-R_R*R_L.diff(u)+R_R.diff(u)*R_L,False)

check('exact_rhs_decomposition',res,False)
check('order_epsilon_squared_connected_cancellation',C2)
check('coincident_A1',at(f1),False)

def dx(M):
    return sum([M.diff(S[k])*X[k]+M.diff(X[k])*Y[k] for k in range(3)],s.zeros(*M.shape))

def cart(v):
    return s.Matrix([(v[0]+v[1])/2,(v[0]-v[1])/(2*I),v[2]])
def components(v):
    return s.Matrix([v[0]+I*v[1],v[0]-I*v[1],v[2]])

# Section 5 compares the common-site product with its coherent-state symbol.
spin_cart = cart(S)
pauli_product_difference = s.Matrix(3,3,lambda i,j:
    s.trace(rho(S)*sig[i]*sig[j])
    -s.trace(rho(S)*sig[i])*s.trace(rho(S)*sig[j])
    -s.KroneckerDelta(i,j)+spin_cart[i]*spin_cart[j]
    -I*sum(s.LeviCivita(i,j,k)*spin_cart[k] for k in range(3)))
check('single_site_Pauli_product_difference',pauli_product_difference,False)

def Q(v):
    # v is expressed in p,q,z components here.
    return s.Matrix([[v[2]+2*al*la*v[0],v[1]-2*al*la*v[2]],[v[0],-v[2]]])/(4*la)
check('two_site_A1',f1+2*Q(components(cart(a).cross(cart(b)))),False)

W=-2*Q(components(cart(S).cross(cart(X))))
Z=2*I*U*U-I*s.eye(2)/(4*la**2)
D=-Z
V=W+Z
check('coincident_A2',at(f2)-2*Z)
check('time_symbol_decomposition',Vdown-W-2*Z)
check('quantum_derived_connected_source',C+2*I*dx(U*U))
check('local_correction_commutes_with_U',comm(D,U),False)
check('local_correction_absorbs_source',C-dx(D)-comm(D,U))
Xi = at(low(lr*lr))-U*U
# Resolve the two contributions in the manuscript's Taylor extraction.
# Each left-hand side comes from the independent quantum product contractions.
check('right_quadratic_overlap_at_equal_spins',at(CR2)-2*I*Xi)
check('left_quadratic_overlap_at_equal_spins',at(CL2)-2*I*Xi)
check('quadratic_overlap_Taylor_contribution',
      at(d(CR2,b,X)+d(CL2,a,X))-2*I*dx(Xi))
check('coincident_cubic_overlap_difference',at(CR3-CL3))
check('source_equals_covariant_second_moment',C-2*I*(dx(Xi)+comm(Xi,U)))
check('local_correction_from_second_moment',D-2*I*Xi+I*Id/(4*la**2))

p,q,z=S; px,qx,zx=X; pxx,qxx,zxx=Y
flow=s.Matrix([2*I*(p*zxx-z*pxx+al*p*px),
               2*I*(z*qxx-q*zxx+al*p*qx),
               I*(q*pxx-p*qxx+2*al*p*zx)])
check('quantum_rhs_equals_continuum_flow',E3-Q(flow))
check('unit_spin_flow',s.Matrix([q*flow[0]+p*flow[1]+2*z*flow[2]]))
check('corrected_curvature',Q(flow)-dx(V)+comm(U,V))
pt,qt,zt=s.symbols('p_T q_T z_T')
F=Q(s.Matrix([pt,qt,zt]))-dx(V)+comm(U,V)
check('off_shell_EOM_factorization',F-Q(s.Matrix([pt,qt,zt])-flow))

# Hamiltonian verification using the explicitly declared continuum Poisson tensor.
H=-(px*qx+zx*zx)/2+al*((1+z)*px-p*zx)/2
variations=s.Matrix([H.diff(S[k])-sum(H.diff(X[k]).diff(S[m])*X[m]
                         +H.diff(X[k]).diff(X[m])*Y[m] for m in range(3)) for k in range(3)])
Pi=s.Matrix([[0,-4*I*z,2*I*p],[4*I*z,0,-2*I*q],[-2*I*p,2*I*q,0]])
check('Hamiltonian_flow',Pi*variations-flow)
# Coefficient of the energy-density Taylor expansion, independently obtained from h.
h_01=2*P+2*s.Symbol('eps')*al*K
eps=s.Symbol('eps')
state=s.kronecker_product(rho(S),rho(S+eps*X+eps**2*Y/2))
h_symbol=s.expand(s.trace(state*h_01))
check('Hamiltonian_density_expansion',s.Matrix([s.expand(h_symbol).coeff(eps,2)-H]))
check('Hamiltonian_density_order_zero',s.Matrix([s.expand(h_symbol).coeff(eps,0)-2]))
check('Hamiltonian_density_order_one',s.Matrix([s.expand(h_symbol).coeff(eps,1)]))

# Explicit EOM dictionary with the conventions of arXiv:2506.13598v2 (2.22).
reflect=s.diag(1,-1,1)
Sref=reflect*cart(S); Xref=reflect*cart(X); Yref=reflect*cart(Y)
Mref=s.Matrix([[0,0,-1],[0,0,I],[1,-I,0]])
reference_rhs=Sref.cross(Yref+al*Mref*Xref) # alpha_reference=-alpha
check('reference_EOM_dictionary',reflect*cart(flow)-2*reference_rhs)
check('XXX_reduction',V.subs(al,0)-(-2*Q(components(cart(S).cross(cart(X)))).subs(al,0)+I*(2*rho(S)-s.eye(2))/(4*la**2)))
isotropic_spin = sum((cart(S)[k]*sig[k] for k in range(3)),s.zeros(2))
isotropic_spin_x = sum((cart(X)[k]*sig[k] for k in range(3)),s.zeros(2))
check('XXX_local_correction',D.subs(al,0)+I*isotropic_spin/(4*la**2))
check('XXX_connected_source',C.subs(al,0)+I*isotropic_spin_x/(4*la**2))

summary={"sympy_version":s.__version__,"python_version":sys.version.split()[0],
         "number_of_checks":len(checks),"checks":checks,
         "scope":"spin-1/2 Class 5, a1=1, a2=2a, a3=0, a=epsilon*alpha, u=lambda/epsilon",
         "constraint_polynomials":[str(t) for t in polys],
         "groebner_variable_order":[str(t) for t in vars],
         "groebner_domain":"QQ_I(lambda,alpha)"}
(out/'symbolic_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
expressions={name:str(M.applyfunc(s.factor)) for name,M in {
    "spatial_U":U,"A1_symbol":f1,"A2_symbol":f2,"raw_time_symbol":Vdown,
    "quantum_connected_source":C,"local_correction":D,"corrected_time_component":V,
    "flow_p_q_z":flow}.items()}
(out/'symbolic_expressions.json').write_text(json.dumps(expressions,indent=2)+'\n')
print("All %d exact checks passed." % len(checks))
