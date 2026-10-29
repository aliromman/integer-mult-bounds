"""Search a mixed three-involution query contract in all four characters.

Numerical discovery only. Every accepted construction needs exact replay.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares

CHARACTERS=np.array([[1,1,1],[-1,1,-1],[1,-1,-1],[-1,-1,1]])


def model(n, labels, targets):
    length=len(labels); ident=np.eye(n)
    factors=CHARACTERS[:, labels].T-1
    target=np.array([np.diag(c[targets]) for c in CHARACTERS])
    def evaluate(x, derivative=True):
        values=x.reshape(length,2,n); us,vs=values[:,0],values[:,1]
        matrices=ident+factors[:,:,None,None]*us[:,None,:,None]*vs[:,None,None,:]
        prefix=[np.broadcast_to(ident,(4,n,n))]
        for matrix in matrices: prefix.append(matrix@prefix[-1])
        f=np.r_[(prefix[-1]-target).ravel(),np.sum(us*vs,axis=1)]
        if not derivative: return f
        jac=np.zeros((len(f),len(x)));after=np.broadcast_to(ident,(4,n,n))
        for site in range(length-1,-1,-1):
            u,v=us[site],vs[site]
            right=v@prefix[site];left=after@u
            for a in range(n):
                jac[:4*n*n,site*2*n+a]=(factors[site,:,None,None]*after[:,:,a,None]*right[:,None,:]).ravel()
                jac[:4*n*n,site*2*n+n+a]=(factors[site,:,None,None]*left[:,:,None]*prefix[site][:,a,None,:]).ravel()
                jac[4*n*n+site,site*2*n+a]=v[a]
                jac[4*n*n+site,site*2*n+n+a]=u[a]
            after=after@matrices[site]
        return f,jac
    return evaluate


def canonical_triangle(offset,n):
    us=[[1,0],[0,1],[-1,1]];vs=[[0,-1],[1,0],[1,1]]
    rows=[]
    for u,v in zip(us,vs):
        a=np.zeros((2,n));a[0,offset:offset+2]=u;a[1,offset:offset+2]=v;rows.append(a)
    return rows


def controls():
    n=6;labels=np.repeat(np.arange(3),3);targets=np.repeat(np.arange(3),2)
    x=np.array(sum((canonical_triangle(2*i,n) for i in range(3)),[])).ravel()
    if np.max(abs(model(n,labels,targets)(x,False)))!=0: raise ValueError('The exact triangle control fails')
    if np.max(abs(model(n,labels,np.roll(targets,1))(x,False)))==0: raise ValueError('Changed targets pass')
    rng=np.random.default_rng(290641);x=x+rng.normal(0,.01,len(x));calc=model(n,labels,targets);_,j=calc(x)
    for col in np.linspace(0,len(x)-1,24,dtype=int):
        step=np.zeros_like(x);step[col]=1e-6
        if np.max(abs((calc(x+step,False)-calc(x-step,False))/2e-6-j[:,col]))>1e-7: raise ValueError('The Jacobian control fails')
    return {'three_triangle_control':'PASS','changed_targets_rejected':True,'finite_difference_columns':24}


def run(args):
    checks=controls(); n=args.n; length=n+args.extra
    targets=np.repeat(np.arange(3),n//3)
    counts=[length//3+(i<length%3) for i in range(3)]
    labels=np.repeat(np.arange(3),counts)
    rng=np.random.default_rng(args.seed); start=time.monotonic();records=[];best=None
    for trial in range(args.trials):
        if args.resume:
            previous=json.loads(args.resume.read_text())['best'];labels=np.array(previous['labels']);x=np.array(previous['parameters'])
            if trial: x+=rng.normal(0,.03,len(x))
        else:
            labels=rng.permutation(labels)
            x=rng.normal(0,.3,(length,2,n)).ravel()
        calc=model(n,labels,targets)
        result=least_squares(lambda x:calc(x,False),x,jac=lambda x:calc(x)[1],max_nfev=args.steps,
                bounds=(-args.bound,args.bound),ftol=1e-11,xtol=1e-11,gtol=1e-11)
        row={'trial':trial,'labels':labels.tolist(),'norm':float(np.linalg.norm(result.fun)),
             'max_residual':float(max(abs(result.fun))),'max_parameter':float(max(abs(result.x))),
             'evaluations':result.nfev,'elapsed_seconds':time.monotonic()-start}
        print(json.dumps(row),flush=True);records.append(row)
        if best is None or row['norm']<best['norm']: best=row|{'parameters':result.x.tolist()}
        output={'fields':n,'queries':length,'target_labels':targets.tolist(),'seed':args.seed,'bound':args.bound,
                'steps':args.steps,'controls':checks,'records':records,'best':best,'status':'NUMERICAL_DISCOVERY_ONLY','new_kappa':False}
        args.output.write_text(json.dumps(output,indent=2)+'\n')
        if row['max_residual']<1e-10: break

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=6)
    p.add_argument('--extra',type=int,default=2);p.add_argument('--trials',type=int,default=12)
    p.add_argument('--steps',type=int,default=500);p.add_argument('--bound',type=float,default=8)
    p.add_argument('--seed',type=int,default=290641);p.add_argument('--resume',type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.n%6 or a.n<6: p.error('Use a positive multiple of six fields')
    run(a)
