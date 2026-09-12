#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd,torch
from class5_lax import RDTYPE,class6_curvature,class6_eom_rhs,class6_exact_coefficients,rms,sample_complex_jets

def train_one(a,alpha,seed,ntrain,nvalid,steps):
    torch.manual_seed(seed)
    s,sx,sxx=sample_complex_jets(ntrain,700000+seed); st=class6_eom_rhs(s,sx,sxx,alpha)
    sv,sxv,sxxv=sample_complex_jets(nvalid,800000+seed); stv=class6_eom_rhs(sv,sxv,sxxv,alpha)
    # Start near the undeformed LL point but do not supply deformation dependence.
    init=np.array([0,0,-a*a,0,0,a,0,0],float)+0.08*np.random.default_rng(seed).normal(size=8)
    p=torch.nn.Parameter(torch.tensor(init,dtype=RDTYPE)); opt=torch.optim.Adam([p],lr=2e-2)
    t0=time.perf_counter(); best=None; bestloss=float('inf')
    for _ in range(steps):
        opt.zero_grad(set_to_none=True); f=class6_curvature(s,sx,sxx,st,a,p)
        loss=torch.mean(torch.abs(f)**2); loss.backward(); torch.nn.utils.clip_grad_norm_([p],10); opt.step()
        val=float(loss.detach())
        if val<bestloss: bestloss=val; best=p.detach().clone()
    if best is not None: p.data.copy_(best)
    lb=torch.optim.LBFGS([p],lr=.7,max_iter=120,tolerance_grad=1e-13,tolerance_change=1e-15,line_search_fn='strong_wolfe')
    def closure():
        lb.zero_grad(set_to_none=True); f=class6_curvature(s,sx,sxx,st,a,p); loss=torch.mean(torch.abs(f)**2); loss.backward(); return loss
    lb.step(closure)
    with torch.no_grad(): tr=class6_curvature(s,sx,sxx,st,a,p); va=class6_curvature(sv,sxv,sxxv,stv,a,p)
    learned=p.detach().numpy(); exact=class6_exact_coefficients(a,alpha); names=['b','c','d','e','f','g','h','j']
    row={'a':a,'alpha':alpha,'seed':seed,'elapsed_s':time.perf_counter()-t0,'train_rms':rms(tr),'valid_rms':rms(va),'max_coefficient_error':float(np.max(np.abs(learned-exact)))}
    for name,v,t in zip(names,learned,exact,strict=True): row[f'{name}_learned']=float(v); row[f'{name}_exact']=float(t); row[f'{name}_error']=float(v-t)
    return row

def main():
    q=argparse.ArgumentParser(); q.add_argument('--out',type=Path,required=True); q.add_argument('--n-train',type=int,default=768); q.add_argument('--n-valid',type=int,default=3072); q.add_argument('--steps',type=int,default=1300); args=q.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for alpha,a in [(0.35,.65),(.7,1.0),(1.1,1.45)]:
        for seed in [0,1]:
            row=train_one(a,alpha,100+seed+int(alpha*100),args.n_train,args.n_valid,args.steps); rows.append(row)
            print(f"alpha={alpha} a={a} seed={seed} valid={row['valid_rms']:.3e} err={row['max_coefficient_error']:.3e}",flush=True)
    df=pd.DataFrame(rows); df.to_csv(args.out/'class6_ml_coefficients.csv',index=False)
    summary={'runs':len(df),'max_validation_rms':float(df.valid_rms.max()),'median_validation_rms':float(df.valid_rms.median()),'max_coefficient_error':float(df.max_coefficient_error.max()),'median_coefficient_error':float(df.max_coefficient_error.median())}
    (args.out/'class6_ml_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
