"""Independent finite-matrix/finite-spacing tests for the Class-6 Lax limit.

Builds R and its inverse at finite lattice spacing; does NOT load the symbolic
program or its output. Tensor spaces and continuum time use the conventions
stated in the accompanying report. NumPy and SymPy are the only dependencies.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parents[2] / 'results' / 'class6'
OUT.mkdir(exist_ok=True)
P = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
K = np.array([[0,1,1,0],[0,0,0,-1],[0,0,0,-1],[0,0,0,0]],complex)
K2 = K@K
SIGMA = np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
N = np.array([[0,0,-1],[0,0,-1j],[-1,-1j,0]],complex)
I2 = np.eye(2,dtype=complex)
I4 = np.eye(4,dtype=complex)

def comm(A,B):
    return A@B-B@A

def embed_pair(mat, first, second, count=3):
    """Ordered two-qubit operator in a tensor product of count qubits."""
    if first == second or min(first,second)<0 or max(first,second)>=count:
        raise ValueError('Invalid tensor positions.')
    size=2**count
    out=np.zeros((size,size),complex)
    for row in range(size):
        r=[(row>>k)&1 for k in reversed(range(count))]
        for col in range(size):
            c=[(col>>k)&1 for k in reversed(range(count))]
            if all(r[k]==c[k] for k in range(count) if k not in(first,second)):
                out[row,col]=mat[2*r[first]+r[second],2*c[first]+c[second]]
    return out

def R(u,a):
    return 2*u*I4+P+a*u*(1+2*u)*K+a*a*u*u*(1+2*u)**2*K2/2

def Rprime(u,a):
    return 2*I4+a*(1+4*u)*K+a*a*u*(1+2*u)*(1+4*u)*K2

def rho(s):
    return (I2+np.einsum('a,aij->ij',s,SIGMA))/2

def lower(op,left,right):
    weight=np.kron(rho(left),rho(right))
    return np.array([[np.trace(weight@op[4*i:4*i+4,4*j:4*j+4])
        for j in range(2)] for i in range(2)])

def linear_map(v,lam,alpha):
    p,q,z=v[0]+1j*v[1],v[0]-1j*v[1],v[2]
    return np.array([[z+2*alpha*lam**2*p,q+4*alpha*lam**2*z-4*alpha**2*lam**4*p],
                     [p,-z-2*alpha*lam**2*p]],complex)/(4*lam)

def linear_map_derivative(v,lam,alpha):
    p,q,z=v[0]+1j*v[1],v[0]-1j*v[1],v[2]
    return np.array([[alpha*p/2-z/(4*lam**2),
                     -3*alpha**2*lam**2*p+alpha*z-q/(4*lam**2)],
                     [-p/(4*lam**2),-alpha*p/2+z/(4*lam**2)]],complex)

def second_spatial_map(v,lam,alpha):
    p,z=v[0]+1j*v[1],v[2]
    return np.array([[alpha*p/4,-alpha**2*lam**2*p+alpha*z/2],[0,-alpha*p/4]],complex)

def third_spatial_map(v,lam,alpha):
    p=v[0]+1j*v[1]
    return np.array([[0,-alpha**2*lam*p/4],[0,0]],complex)

def continuum(s,sx,sxx,lam,alpha):
    U=I2/(4*lam)+linear_map(s,lam,alpha)
    Ux=linear_map(sx,lam,alpha)
    M=second_spatial_map(s,lam,alpha)
    Mx=second_spatial_map(sx,lam,alpha)
    W=-2*linear_map(np.cross(s,sx),lam,alpha)
    Wx=-2*linear_map(np.cross(s,sxx),lam,alpha)
    qprime=linear_map_derivative(s,lam,alpha)
    qprimex=linear_map_derivative(sx,lam,alpha)
    p,z=s[0]+1j*s[1],s[2]
    px,zx=sx[0]+1j*sx[1],sx[2]
    Vdown=W-2j*qprime+1j*alpha*p*z*I2
    Vdownx=Wx-2j*qprimex+1j*alpha*(px*z+p*zx)*I2
    D=4j*M-2j*U@U+1j*I2/(4*lam**2)
    Dx=4j*Mx-2j*(Ux@U+U@Ux)
    source=Dx+comm(D,U)
    flow=-2*np.cross(s,sxx-alpha*N@s)
    Ut=linear_map(flow,lam,alpha)
    naive=Vdownx+comm(Vdown,U)
    V=W-1j*qprime
    Vx=Wx-1j*qprimex
    residual=np.linalg.norm(Ut-naive-source)
    scale=max(1,np.linalg.norm(Ut),np.linalg.norm(naive),np.linalg.norm(source))
    if residual/scale > 1e-11:
        raise AssertionError(('continuum formula consistency',residual/scale))
    return dict(U=U,V=V,Vx=Vx,Vdown=Vdown,source=source,Ut=Ut,flow=flow)

# Smooth unit-vector profile; no numerical differentiation of spatial data.
x=sp.symbols('x',real=True)
theta=sp.Rational(9,10)+sp.sin(x)/5
phi=sp.Rational(3,10)+2*x/5+sp.cos(2*x)/10
profile=sp.Matrix([sp.sin(theta)*sp.cos(phi),sp.sin(theta)*sp.sin(phi),sp.cos(theta)])
profile_functions=[sp.lambdify(x,profile.diff(x,k),'numpy') for k in range(3)]
def spin(x0,order=0):
    return np.asarray(profile_functions[order](x0),complex).reshape(3)

def finite(eps,alpha,lam,x0):
    a,u=eps**2*alpha,lam/eps
    raw=R(u,a)
    normalized=raw/(2*u)
    # Derivative at FIXED quantum coupling a, not fixed a*u or a*u^2.
    derivative=-P/(2*u*u)+a*K+a*a*(1+8*u+12*u*u)*K2/4
    LL=embed_pair(normalized,0,1)
    LR=embed_pair(normalized,0,2)
    dLR=embed_pair(derivative,0,2)
    h=embed_pair(2*P+a*K,1,2)
    A=-1j*np.linalg.solve(LR,comm(h,LR)+dLR)
    sm,s0,sp_=spin(x0-eps),spin(x0),spin(x0+eps)
    Ar,Al=lower(A,s0,sp_),lower(A,sm,s0)
    l0=lower(LL,s0,sp_)
    exact=lower(A@LL,s0,sp_)-lower(LR@A,sm,s0)
    naive=Ar@l0-l0@Al
    return dict(time_symbol=Al/eps**2,exact=exact/eps**3,
                naive=naive/eps**3,connected=(exact-naive)/eps**3)

x0=0.37
lam=0.8+0.2j
rows=[]
limits={}
finite_hamiltonian_residuals=[]
for model,alpha in [('XXX',0j),('Class6',0.7-0.2j)]:
    cv=continuum(spin(x0),spin(x0,1),spin(x0,2),lam,alpha)
    limits[model]={'alpha':[alpha.real,alpha.imag],
                   'connected_source_norm':float(np.linalg.norm(cv['source']))}
    for eps in [0.16,0.08,0.04,0.02,0.01,0.005]:
        fv=finite(eps,alpha,lam,x0)
        s0=spin(x0); neighbors=spin(x0-eps)+spin(x0+eps)
        lattice_flow=-2*np.cross(s0,(neighbors-2*s0)/eps**2-alpha*(N@neighbors)/2)
        lattice_target=linear_map(lattice_flow,lam,alpha)+eps*second_spatial_map(lattice_flow,lam,alpha)+eps**2*third_spatial_map(lattice_flow,lam,alpha)
        finite_hamiltonian_residuals.append(float(np.linalg.norm(fv['exact']-lattice_target)/max(1,np.linalg.norm(lattice_target))))
        # A finite-epsilon spatial symbol also has m_T and n_T contributions.
        finite_target=cv['Ut']+eps*second_spatial_map(cv['flow'],lam,alpha)\
                              +eps**2*third_spatial_map(cv['flow'],lam,alpha)
        rows.append(dict(model=model,epsilon=eps,
            time_symbol_error=float(np.linalg.norm(fv['time_symbol']-cv['Vdown'])),
            exact_leading_error=float(np.linalg.norm(fv['exact']-cv['Ut'])),
            exact_full_spatial_error=float(np.linalg.norm(fv['exact']-finite_target)),
            connected_error=float(np.linalg.norm(fv['connected']-cv['source'])),
            factorized_error=float(np.linalg.norm(fv['naive']-cv['Ut']))))
    last=[r for r in rows if r['model']==model][-3:]
    for label in ['exact_leading_error','exact_full_spatial_error','connected_error']:
        limits[model][label+'_order']=float(np.polyfit(np.log([r['epsilon'] for r in last]),
                                np.log([r[label] for r in last]),1)[0])

if max(finite_hamiltonian_residuals)>1e-8:
    raise AssertionError('Finite Hamiltonian-symbol identity failed.')

rng=np.random.default_rng(20260909)
zc_res=[]
ybe_res=[]
for _ in range(32):
    a=(rng.normal()+1j*rng.normal())/4
    u=0.8+0.3j+(rng.normal()+1j*rng.normal())/5
    raw,rp=R(u,a),Rprime(u,a)
    h=2*P+a*K
    Lc,Lr=embed_pair(raw,0,2,4),embed_pair(raw,0,3,4)
    dLc,dLr=embed_pair(rp,0,2,4),embed_pair(rp,0,3,4)
    hl,hr=embed_pair(h,1,2,4),embed_pair(h,2,3,4)
    Ac=1j*hl-1j*np.linalg.solve(Lc,hl@Lc+dLc)
    Ar=1j*hr-1j*np.linalg.solve(Lr,hr@Lr+dLr)
    lhs=1j*comm(hl+hr,Lc)
    rhs=Ar@Lc-Lc@Ac
    zc_res.append(float(np.linalg.norm(lhs-rhs)/max(1,np.linalg.norm(lhs),np.linalg.norm(rhs))))
    v=(rng.normal()+1j*rng.normal())/3
    r12,r13,r23=embed_pair(R(u-v,a),0,1),embed_pair(R(u,a),0,2),embed_pair(R(v,a),1,2)
    lhs,rhs=r12@r13@r23,r23@r13@r12
    ybe_res.append(float(np.linalg.norm(lhs-rhs)/max(1,np.linalg.norm(lhs),np.linalg.norm(rhs))))
if max(zc_res)>1e-11 or max(ybe_res)>1e-11:
    raise AssertionError('Finite operator identity failed.')

# Off-shell tests: time derivatives are sampled independently of the EOM.
off_res=[]
for _ in range(128):
    s=rng.normal(size=3);s/=np.linalg.norm(s)
    sx=rng.normal(size=3);sx-=s*np.dot(s,sx)
    sxx=rng.normal(size=3);sxx+=s*(-np.dot(sx,sx)-np.dot(s,sxx))
    st=rng.normal(size=3);st-=s*np.dot(s,st)
    alpha=(rng.normal()+1j*rng.normal())/3
    spectral=0.7+0.2j+(rng.normal()+1j*rng.normal())/6
    cv=continuum(s,sx,sxx,spectral,alpha)
    curvature=linear_map(st,spectral,alpha)-cv['Vx']+comm(cv['U'],cv['V'])
    predicted=linear_map(st-cv['flow'],spectral,alpha)
    off_res.append(float(np.linalg.norm(curvature-predicted)/max(1,np.linalg.norm(curvature),np.linalg.norm(predicted))))
if max(off_res)>1e-11:
    raise AssertionError('Off-shell factorization failed.')

uniform=continuum(np.array([0,0,1]),np.zeros(3),np.zeros(3),lam,0.7-0.2j)
uniform_norm=float(np.linalg.norm(uniform['source']))
summary={'python_version':sys.version.split()[0],'numpy_version':np.__version__,
         'sympy_version':sp.__version__,'seed':20260909,
         'finite_zero_curvature_cases':len(zc_res),'finite_zero_curvature_max_relative_error':max(zc_res),
         'Yang_Baxter_cases':len(ybe_res),'Yang_Baxter_max_relative_error':max(ybe_res),
         'off_shell_cases':len(off_res),'off_shell_max_relative_error':max(off_res),
         'finite_Hamiltonian_match_max_relative_error':max(finite_hamiltonian_residuals),
         'profile_x':x0,'spectral_parameter':[lam.real,lam.imag],
         'uniform_source_norm':uniform_norm,'limits':limits,
         'norm':'Frobenius; relative operator error divided by max(1,norm(lhs),norm(rhs)); finite-spacing errors absolute'}
(OUT/'numeric_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
with (OUT/'finite_spacing.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(json.dumps(summary,indent=2))
for row in rows:print(row)
