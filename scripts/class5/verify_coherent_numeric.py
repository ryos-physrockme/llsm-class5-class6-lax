"""Independent finite-spacing checks for the Class-5 spin-1/2 Lax limit.
Run with Python 3, NumPy and SymPy. Results are written next to this directory.
No previously computed continuum time Lax component is imported.
"""
from pathlib import Path
import csv, json
import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parents[2] / 'results' / 'class5'
OUT.mkdir(exist_ok=True)
P = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], complex)
K = np.array([[0,1,-1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]], complex)
SIGMA = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)

def embed_pair(mat, first, second, count=3):
    """Embed an ordered two-qubit operator into a count-qubit tensor product."""
    size = 2**count
    out = np.zeros((size,size), complex)
    for row in range(size):
        r = [(row >> k)&1 for k in reversed(range(count))]
        for col in range(size):
            c = [(col >> k)&1 for k in reversed(range(count))]
            if all(r[k] == c[k] for k in range(count) if k not in (first,second)):
                out[row,col] = mat[2*r[first]+r[second],2*c[first]+c[second]]
    return out

PL = embed_pair(P,0,1); PR = embed_pair(P,0,2)
KL = embed_pair(K,0,1); KR = embed_pair(K,0,2)
PH = embed_pair(P,1,2); KH = embed_pair(K,1,2)
I8 = np.eye(8, dtype=complex)

def rho(s):
    return (np.eye(2)+np.einsum('a,aij->ij',s,SIGMA))/2

def lower(op, left, right):
    weight = np.kron(rho(left),rho(right))
    return np.array([[np.trace(weight @ op[4*i:4*i+4,4*j:4*j+4])
                      for j in range(2)] for i in range(2)])

def comm(a,b):
    return a@b-b@a

def spatial(s, lam, alpha):
    p = s[0]+1j*s[1]; q = s[0]-1j*s[1]; z=s[2]
    return np.array([[1+z+2*alpha*lam*p, q-2*alpha*lam*(1+z)],
                     [p,1-z]], complex)/(4*lam)

def linear_map(v, lam, alpha):
    p=v[0]+1j*v[1]; q=v[0]-1j*v[1]; z=v[2]
    return np.array([[z+2*alpha*lam*p, q-2*alpha*lam*z], [p,-z]], complex)/(4*lam)

# Smooth unit-vector profile, with its derivatives generated independently.
x = sp.symbols('x', real=True)
theta = sp.Rational(9,10)+sp.sin(x)/5
phi = sp.Rational(3,10)+2*x/5+sp.cos(2*x)/10
profile_expr = sp.Matrix([sp.sin(theta)*sp.cos(phi),sp.sin(theta)*sp.sin(phi),sp.cos(theta)])
profile = [sp.lambdify(x, profile_expr.diff(x,k), 'numpy') for k in range(3)]
def spin(x0, order=0):
    return np.asarray(profile[order](x0), dtype=complex).reshape(3)

def continuum(s, sx, sxx, lam, alpha):
    U=spatial(s,lam,alpha); Ux=linear_map(sx,lam,alpha)
    W=-2*linear_map(np.cross(s,sx),lam,alpha)
    Wx=-2*linear_map(np.cross(s,sxx),lam,alpha)
    Z=2j*U@U-1j*np.eye(2)/(4*lam**2)
    Zx=2j*(Ux@U+U@Ux)
    Vdown=W+2*Z
    naive=Wx+2*Zx+comm(Vdown,U)
    C=-Zx
    p,q,z=s[0]+1j*s[1],s[0]-1j*s[1],s[2]
    px,qx,zx=sx[0]+1j*sx[1],sx[0]-1j*sx[1],sx[2]
    pxx,qxx,zxx=sxx[0]+1j*sxx[1],sxx[0]-1j*sxx[1],sxx[2]
    pt=2j*(p*zxx-z*pxx+alpha*p*px)
    qt=2j*(z*qxx-q*zxx+alpha*p*qx)
    zt=1j*(q*pxx-p*qxx+2*alpha*p*zx)
    st=np.array([(pt+qt)/2,(pt-qt)/(2j),zt])
    Ut=linear_map(st,lam,alpha)
    assert np.linalg.norm(Ut-naive-C)<1e-12
    return U,Vdown,C,naive,Ut

def finite(eps, alpha, lam, x0):
    h=2*PH+2*eps*alpha*KH
    LL=I8+eps*(PL/(2*lam)+alpha*KL)
    LR=I8+eps*(PR/(2*lam)+alpha*KR)
    # u derivative at fixed lattice coupling a, followed by u=lambda/eps.
    derivative=-eps**2*PR/(2*lam**2)
    A=-1j*np.linalg.solve(LR,comm(h,LR)+derivative)
    sm,s0,sp_=spin(x0-eps),spin(x0),spin(x0+eps)
    Ar=lower(A,s0,sp_); Al=lower(A,sm,s0)
    l0=lower(LL,s0,sp_)
    exact=lower(A@LL,s0,sp_)-lower(LR@A,sm,s0)
    naive=Ar@l0-l0@Al
    return Al/eps**2, exact/eps**3, naive/eps**3, (exact-naive)/eps**3

x0=0.37; lam=0.8+0.2j
rows=[]; limits={}
for name,alpha in [('XXX',0j),('Class5',0.7-0.2j)]:
    U,Vdown,C,naive,Ut=continuum(spin(x0),spin(x0,1),spin(x0,2),lam,alpha)
    limits[name]={'correction_norm':float(np.linalg.norm(C)),
                  'naive_continuum_error':float(np.linalg.norm(naive-Ut))}
    for eps in [0.16,0.08,0.04,0.02,0.01,0.005]:
        vd,ex,na,co=finite(eps,alpha,lam,x0)
        rows.append(dict(model=name,epsilon=eps,
            time_symbol_error=float(np.linalg.norm(vd-Vdown)),
            exact_flow_error=float(np.linalg.norm(ex-Ut)),
            connected_error=float(np.linalg.norm(co-C)),
            naive_flow_error=float(np.linalg.norm(na-Ut))))
    last=[r for r in rows if r['model']==name][-3:]
    limits[name]['exact_flow_order']=float(np.polyfit(np.log([r['epsilon'] for r in last]),
                                            np.log([r['exact_flow_error'] for r in last]),1)[0])
    limits[name]['connected_order']=float(np.polyfit(np.log([r['epsilon'] for r in last]),
                                            np.log([r['connected_error'] for r in last]),1)[0])

# Exact discrete zero curvature, on auxiliary + left + central + right quantum spaces.
# Only the two Hamiltonian bonds not commuting with the central Lax operator enter.
rng=np.random.default_rng(20260909)
residuals=[]
for _ in range(24):
    a=(rng.normal()+1j*rng.normal())/4
    u=0.8+0.3j+(rng.normal()+1j*rng.normal())/5
    R=2*u*np.eye(4)+P+2*a*u*K
    Rp=2*np.eye(4)+2*a*K
    hh=2*P+2*a*K
    Lc=embed_pair(R,0,2,4); Lr=embed_pair(R,0,3,4)
    Lpc=embed_pair(Rp,0,2,4); Lpr=embed_pair(Rp,0,3,4)
    hl=embed_pair(hh,1,2,4); hr=embed_pair(hh,2,3,4)
    Ac=1j*hl-1j*np.linalg.solve(Lc,hl@Lc+Lpc)
    Ar=1j*hr-1j*np.linalg.solve(Lr,hr@Lr+Lpr)
    lhs=1j*comm(hl+hr,Lc); rhs=Ar@Lc-Lc@Ac
    residuals.append(float(np.linalg.norm(lhs-rhs)/max(1,np.linalg.norm(lhs),np.linalg.norm(rhs))))
assert max(residuals)<1e-12

# Noncommuting auxiliary coefficients are retained in the shared-site contraction.
star_res=[]
for _ in range(32):
    v=rng.normal(size=3); v/=np.linalg.norm(v)
    aa=rng.normal(size=(4,2,2))+1j*rng.normal(size=(4,2,2))
    bb=rng.normal(size=(4,2,2))+1j*rng.normal(size=(4,2,2))
    basis=np.concatenate([np.eye(2,dtype=complex)[None],SIGMA],axis=0)
    A=sum(np.kron(aa[k],basis[k]) for k in range(4))
    B=sum(np.kron(bb[k],basis[k]) for k in range(4))
    l=lambda M:np.array([[np.trace(rho(v)@M[2*i:2*i+2,2*j:2*j+2]) for j in range(2)] for i in range(2)])
    expected=np.zeros((2,2),complex)
    for k in range(3):
        for m in range(3):
            epsilon_term=np.dot(np.cross(np.eye(3)[k],np.eye(3)[m]),v)
            expected+=aa[k+1]@bb[m+1]*((1 if k==m else 0)-v[k]*v[m]+1j*epsilon_term)
    star_res.append(float(np.linalg.norm(l(A@B)-l(A)@l(B)-expected)))
assert max(star_res)<1e-12
with (OUT/'finite_spacing.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'seed':20260909,'local_zero_curvature_cases':len(residuals),
         'max_relative_local_zero_curvature_residual':max(residuals),
         'shared_site_star_cases':len(star_res),'max_absolute_star_residual':max(star_res),
         'profile_x':x0,'lambda':[lam.real,lam.imag], 'limits':limits,
         'numpy_version':np.__version__,'sympy_version':sp.__version__}
(OUT/'numeric_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
for r in rows: print(r)
