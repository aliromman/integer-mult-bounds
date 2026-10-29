"""Check every cyclic stage character of the retained mixed-query chart.

This is numerical discovery. A zero needs exact reconstruction and pricing.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares


def model(data,period):
    pc=data['projectors'];rs=data['queries'];labels=data['labels'];block=int(data['block'])
    length,terms,n,_=pc.shape;count=terms-1;ident=np.eye(n);dh=pc[:,1:].sum(axis=0)
    chars=np.array([[1,1,1],[-1,1,-1],[1,-1,-1],[-1,-1,1]])
    z=np.tile(np.exp(2j*np.pi*np.arange(period)/period),4);character=np.repeat(chars,period,axis=0)
    signs=np.repeat(character,block,axis=1)
    ts=[ident+(character[:,label]-1)[:,None,None]*r for label,r in zip(labels,rs)]
    spots=[tuple(np.argwhere(r!=0)[0]) for r in rs]
    def calculate(x,derivative=True):
        ps=pc[:,0]+np.einsum('mkij,k->mij',pc[:,1:],x)
        h=ps.sum(axis=0)
        matrices=[t+(z-1)[:,None,None]*(p@t) for p,t in zip(ps,ts)]
        prefix=[np.broadcast_to(ident,(len(z),n,n))]
        for a in matrices:prefix.append(a@prefix[-1])
        target=signs[:,:,None]*(z*z)[:,None,None]*((3-2*z)[:,None,None]*ident+(z-1)[:,None,None]*h)
        projector=np.array([(p@p-p)[i,j] for p,(i,j) in zip(ps,spots)])
        f=np.r_[((prefix[-1]-target)/np.sqrt(len(z))).ravel(),(h@h-5*h+6*ident).ravel(),projector]
        if not derivative:return np.r_[f.real,f.imag]
        j=np.zeros((len(f),count),complex);after=np.broadcast_to(ident,(len(z),n,n));rows=len(z)*n*n
        for site in range(length-1,-1,-1):
            dp=pc[site,1:]
            delta=after[:,None,:,:]@dp[None,:,:,:]@(ts[site]@prefix[site])[:,None,:,:]
            delta*=(z-1)[:,None,None,None]/np.sqrt(len(z))
            j[:rows]+=delta.transpose(0,2,3,1).reshape(rows,count)
            p=ps[site];a,b=spots[site];j[rows+n*n+site]=(p@dp+dp@p-dp)[:,a,b]
            after=after@matrices[site]
        dt=signs[:,None,:,None]*(z*z*(z-1))[:,None,None,None]*dh[None,:,:,:]/np.sqrt(len(z))
        j[:rows]-=dt.transpose(0,2,3,1).reshape(rows,count)
        j[rows:rows+n*n]=(h@dh+dh@h-5*dh).reshape(count,n*n).T
        return np.r_[f.real,f.imag],np.vstack([j.real,j.imag])
    return calculate,count


def positive_controls():
    n=6;ps=[];rs=[];labels=[]
    us=[[1,0],[0,1],[-1,1]];vs=[[0,-1],[1,0],[1,1]]
    for site in range(3):
        for block in range(3):
            p=np.zeros((1,n,n));p[0,2*block:2*block+2,2*block:2*block+2]=np.eye(2)
            r=np.zeros((n,n));r[2*block:2*block+2,2*block:2*block+2]=np.outer(us[site],vs[site])
            ps.append(p);rs.append(r);labels.append(block)
    data={'projectors':np.array(ps),'queries':np.array(rs),'labels':np.array(labels),'block':2}
    calc,_=model(data,6)
    residual=float(max(abs(calc(np.array([]),False))))
    if residual>1e-11:raise ValueError('The three-triangle stage control fails')
    data['queries'][0]*=-1
    changed,_=model(data,6)
    if max(abs(changed(np.array([]),False)))<.1:raise ValueError('A changed query passes')
    return {'three_triangle_max_residual':residual,'changed_query_rejected':True}


def run(args):
    control=positive_controls()
    data=np.load(args.chart);calc,count=model(data,args.period);rng=np.random.default_rng(290643);start=time.monotonic()
    # This check compares a derived endpoint Jacobian with direct perturbations.
    x=rng.normal(0,.005,count);_,j=calc(x)
    for col in np.linspace(0,count-1,min(24,count),dtype=int):
        dx=np.zeros(count);dx[col]=1e-6
        diff=(calc(x+dx,False)-calc(x-dx,False))/2e-6
        if np.max(abs(diff-j[:,col]))>1e-6:raise ValueError('The full endpoint Jacobian differs')
    print(json.dumps({'controls':control,'Jacobian_columns':'PASS','variables':count}),flush=True)
    records=[];best=None
    for trial in range(args.trials):
        x=np.zeros(count) if trial==0 else rng.normal(0,.05,count)
        if args.resume:x=np.array(json.loads(args.resume.read_text())['best']['parameters'])+(rng.normal(0,.01,count) if trial else 0)
        result=least_squares(lambda x:calc(x,False),x,jac=lambda x:calc(x)[1],max_nfev=args.steps,
                bounds=(-args.bound,args.bound),ftol=1e-11,xtol=1e-11,gtol=1e-11)
        row={'trial':trial,'norm':float(np.linalg.norm(result.fun)),'max_residual':float(max(abs(result.fun))),
             'max_parameter':float(max(abs(result.x))),'evaluations':result.nfev,'elapsed_seconds':time.monotonic()-start}
        records.append(row);print(json.dumps(row),flush=True)
        if best is None or row['norm']<best['norm']:best=row|{'parameters':result.x.tolist()}
        args.output.write_text(json.dumps({'block':int(data['block']),'variables':count,'period':args.period,'records':records,
                'best':best,'bound':args.bound,'steps':args.steps,'status':'NUMERICAL_DISCOVERY_ONLY','new_kappa':False},indent=2)+'\n')
        if row['max_residual']<1e-10:break

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--chart',type=Path,required=True)
    p.add_argument('--period',type=int,default=6);p.add_argument('--trials',type=int,default=2)
    p.add_argument('--steps',type=int,default=200);p.add_argument('--bound',type=float,default=8)
    p.add_argument('--resume',type=Path);p.add_argument('--output',type=Path,required=True);run(p.parse_args())
