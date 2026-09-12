#!/usr/bin/env python3
"""Original eight-coefficient search in the manuscript conventions.

Adapted from archive/ml_original/run_class6_ml_discovery.py, preserving
the original spatial sampling, initialization, Adam iterate selection and L-BFGS.
The default off-shell loss is exactly equivalent for this ansatz; select
--training-mode on-shell to run the original curvature objective directly.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np,pandas as pd,torch
from paper_conventions import (RDTYPE, class6_training_residual, class6_eom_rhs,
    class6_exact_coefficients, rms, sample_complex_jets, spectral_coefficient,
    sample_class6_time_derivative)

torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)

def train_one(lambda_,alpha6,seed,ntrain,nvalid,steps,training_mode="off-shell"):
    torch.manual_seed(seed)
    s,sx,sxx=sample_complex_jets(ntrain,700000+seed)
    sv,sxv,sxxv=sample_complex_jets(nvalid,800000+seed)
    if training_mode == "off-shell":
        st=sample_class6_time_derivative(s,900000+seed)
        stv=sample_class6_time_derivative(sv,1000000+seed)
    else:
        st=class6_eom_rhs(s,sx,sxx,alpha6)
        stv=class6_eom_rhs(sv,sxv,sxxv,alpha6)
    def residual(spin,spin_x,spin_xx,spin_t,coefficients):
        return class6_training_residual(spin,spin_x,spin_xx,spin_t,lambda_,coefficients,alpha6,training_mode)
    # Start near the undeformed LL point but do not supply deformation dependence.
    leading=spectral_coefficient(lambda_)
    init=np.array([0,0,-leading**2,0,0,leading,0,0],float)+0.08*np.random.default_rng(seed).normal(size=8)
    p=torch.nn.Parameter(torch.tensor(init,dtype=RDTYPE)); opt=torch.optim.Adam([p],lr=2e-2)
    t0=time.perf_counter(); best=None; bestloss=float('inf')
    for _ in range(steps):
        opt.zero_grad(set_to_none=True); f=residual(s,sx,sxx,st,p)
        loss=torch.mean(torch.abs(f)**2); loss.backward(); torch.nn.utils.clip_grad_norm_([p],10); opt.step()
        val=float(loss.detach())
        if val<bestloss: bestloss=val; best=p.detach().clone()
    if best is not None: p.data.copy_(best)
    lb=torch.optim.LBFGS([p],lr=.7,max_iter=120,tolerance_grad=1e-13,tolerance_change=1e-15,line_search_fn='strong_wolfe')
    def closure():
        lb.zero_grad(set_to_none=True); f=residual(s,sx,sxx,st,p); loss=torch.mean(torch.abs(f)**2); loss.backward(); return loss
    lb.step(closure)
    with torch.no_grad(): tr=residual(s,sx,sxx,st,p); va=residual(sv,sxv,sxxv,stv,p)
    learned=p.detach().numpy(); exact=class6_exact_coefficients(lambda_,alpha6); names=['a1','a2','b0','b1','b2','c0','c1','c2']
    row={'training_mode':training_mode,'lambda_real':complex(lambda_).real,'lambda_imag':complex(lambda_).imag,'alpha6':alpha6,'seed':seed,'elapsed_s':time.perf_counter()-t0,'train_rms':rms(tr),'valid_rms':rms(va),'max_coefficient_error':float(np.max(np.abs(learned-exact)))}
    for name,v,t in zip(names,learned,exact,strict=True): row[f'{name}_learned']=float(v); row[f'{name}_exact']=float(t); row[f'{name}_error']=float(v-t)
    return row

def main():
    q=argparse.ArgumentParser(); q.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[2]/'results/ml/class6'); q.add_argument('--n-train',type=int,default=768); q.add_argument('--n-valid',type=int,default=3072); q.add_argument('--steps',type=int,default=1300)
    q.add_argument('--training-mode',choices=['off-shell','on-shell'],default='off-shell')
    args=q.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for alpha6,lambda_ in [(0.35,10j/13),(.7,0.5j),(1.1,10j/29)]:
        for seed in [0,1]:
            row=train_one(lambda_,alpha6,100+seed+int(alpha6*100),args.n_train,args.n_valid,args.steps,args.training_mode); rows.append(row)
            print(f"alpha6={alpha6} lambda={lambda_} seed={seed} valid={row['valid_rms']:.3e} err={row['max_coefficient_error']:.3e}",flush=True)
    df=pd.DataFrame(rows); df.to_csv(args.out/'class6_ml_coefficients.csv',index=False)
    summary={'protocol':'original-derived search in manuscript conventions; equivalent off-shell extension',
        'training_mode':args.training_mode,
        'torch_version':torch.__version__,'train_samples':args.n_train,'validation_samples':args.n_valid,
        'adam_steps':args.steps,'lbfgs_max_iter':120,'learned_real_parameters':8,'runs':len(df),'max_validation_rms':float(df.valid_rms.max()),'median_validation_rms':float(df.valid_rms.median()),'max_coefficient_error':float(df.max_coefficient_error.max()),'median_coefficient_error':float(df.max_coefficient_error.median())}
    (args.out/'class6_ml_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8'); print(json.dumps(summary,indent=2))
    if summary['max_validation_rms'] > 1e-5 or summary['max_coefficient_error'] > 1e-4:
        raise SystemExit('Class 6 original-protocol recovery failed.')
if __name__=='__main__': main()
