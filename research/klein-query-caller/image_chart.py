"""Search the two-shortest-character image chart.

The chart covers the common-image branch only. It is not a general exclusion.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
import sympy as sp
from scipy.optimize import least_squares


def chart(labels):
    m=len(labels);n=m-2
    ka=sp.zeros(m,1);kb=sp.zeros(m,1)
    for k,missing in [(ka,0),(kb,1)]:
        for pos,i in enumerate(j for j,c in enumerate(labels) if c!=missing): k[i]=(-1)**pos
    dep=[max(i for i,c in enumerate(labels) if c==j) for j in [0,1]]
    basis=[i for i in range(m) if i not in dep]
    u=sp.zeros(n,m)
    for a,i in enumerate(basis):u[a,i]=1
    for k,d in [(kb,dep[0]),(ka,dep[1])]:
        u[:,d]=-sum((k[i]*u[:,i] for i in basis),sp.zeros(n,1))/k[d]
    v0=[];directions=[]
    for i,c in enumerate(labels):
        known=[j for j,t in enumerate(labels) if c==t or c==2 or t==2]
        a=u[:,known].T;b=sp.Matrix([sp.sign(i-j) for j in known]);sol=a.gauss_jordan_solve(b)
        v0.append(sol[0].subs(dict.fromkeys(sol[1],0)).T)
        directions.append(a.nullspace())
    base=sp.Matrix.vstack(*v0);deriv=[]
    for i,rows in enumerate(directions):
        for row in rows:
            dv=sp.zeros(m,n);dv[i,:]=row.T;deriv.append(dv)
    return np.array(u,float),np.array(base,float),np.array([d.tolist() for d in deriv],float)


def model(labels):
    u,v0,dv=chart(labels);n=len(u);m=len(labels);count=len(dv);ident=np.eye(n)
    masks=[[c!=missing for c in labels] for missing in [0,1,2]]
    def word(v,mask,derivative):
        matrices=[ident-2*np.outer(u[:,i],v[i]) if mask[i] else ident for i in range(m)]
        prefix=[ident]
        for a in matrices:prefix.append(a@prefix[-1])
        if not derivative:return prefix[-1]
        ds=np.zeros((count,n,n));after=ident
        for i in range(m-1,-1,-1):
            if mask[i]:ds-=2*(after@u[:,i])[None,:,None]*(dv[:,i,:]@prefix[i])[:,None,:]
            after=after@matrices[i]
        return prefix[-1],ds
    def evaluate(x,derivative=True):
        v=v0+np.einsum('k,kij->ij',x,dv)
        if not derivative:
            a,b,c=[word(v,k,False) for k in masks];return np.r_[(a@b-b@a).ravel(),(c-a@b).ravel()]
        (a,da),(b,db),(c,dc)=[word(v,k,True) for k in masks]
        f=np.r_[(a@b-b@a).ravel(),(c-a@b).ravel()]
        j=np.concatenate([(da@b+a@db-db@a-b@da).reshape(count,-1),(dc-da@b-a@db).reshape(count,-1)],axis=1).T
        return f,j
    return evaluate,count


def necklaces():
    all_words=set(itertools.permutations([0]*3+[1]*3+[2]*2));representatives=set()
    for word in all_words:
        variants=[]
        for x in [word,word[::-1],tuple(1-c if c<2 else c for c in word),tuple(1-c if c<2 else c for c in word[::-1])]:
            variants.extend(x[i:]+x[:i] for i in range(8))
        representatives.add(min(variants))
    return sorted(representatives)


def run(args):
    rng=np.random.default_rng(290642);start=time.monotonic();records=[];best=None
    for index,labels in enumerate(necklaces()):
        calc,count=model(labels);x=rng.normal(0,.1,count);_,j=calc(x)
        for col in range(count):
            dx=np.zeros(count);dx[col]=1e-6
            if np.max(abs((calc(x+dx,False)-calc(x-dx,False))/2e-6-j[:,col]))>1e-7:raise ValueError('Jacobian mismatch')
        for trial in range(args.trials):
            x=np.zeros(count) if trial==0 else rng.normal(0,.2,count)
            result=least_squares(lambda x:calc(x,False),x,jac=lambda x:calc(x)[1],max_nfev=args.steps,
                                bounds=(-8,8),ftol=1e-12,xtol=1e-12,gtol=1e-12)
            row={'labels':labels,'pattern':index,'trial':trial,'variables':count,'norm':float(np.linalg.norm(result.fun)),
                 'max_residual':float(max(abs(result.fun))),'evaluations':result.nfev,'elapsed_seconds':time.monotonic()-start}
            records.append(row);print(json.dumps(row),flush=True)
            if best is None or row['norm']<best['norm']:best=row|{'parameters':result.x.tolist()}
            args.output.write_text(json.dumps({'scope':'common-image chart, real coefficients','patterns':len(necklaces()),
                    'records':records,'best':best,'new_kappa':False},indent=2)+'\n')
            if row['max_residual']<1e-10:return

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--trials',type=int,default=2)
    p.add_argument('--steps',type=int,default=400);p.add_argument('--output',type=Path,required=True);run(p.parse_args())
