"""Search inverse boundary stages for a complete monomial stage endpoint.

The interior and query word are fixed. Boundary projector charts are complete.
Numerical residuals are discovery only. The whole signed invoice is still due.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import least_squares
from history_discovery import projectors
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'free-query-work'))
from search import frames


def model(n, target):
    m=n//2; count=4*(n-2); local=2*(n-2)
    basis=[frames(n,0)[i] for i in (0,n)]
    inverse=[np.linalg.inv(b) for b in basis]
    ps=[np.array(p.tolist(),dtype=complex) for p in projectors(n)]
    initial=[]
    for p,b in zip((ps[0],ps[-1]),inverse):
        core=b@p@np.linalg.inv(b)
        initial.extend(core[2:,1]);initial.extend(core[0,2:])
    initial=np.array(initial)
    da=np.zeros((local,n-2),dtype=complex);db=da.copy()
    da[:n-2]=np.eye(n-2);db[n-2:]=np.eye(n-2)
    roots=np.exp(2j*np.pi*np.arange(m+3)/(m+3))
    degrees=list(range(-2,m+1)); ident=np.eye(n,dtype=complex)
    interior=[]
    for z in roots:
        w=ident.copy()
        for p in ps[1:-1]:w=(ident+(z-1)*p)@w
        interior.append(w)

    def calculate(real):
        values=real[:count]+1j*real[count:]
        ends=[]
        for site in range(2):
            a,b=values[site*local:(site+1)*local].reshape(2,n-2)
            core=np.zeros((n,n),dtype=complex);core[0,0]=core[1,1]=1
            core[0,1]=-b@a;core[0,2:]=b;core[2:,1]=a
            derivative=np.zeros((local,n,n),dtype=complex)
            derivative[:,0,1]=-db@a-da@b
            derivative[:,0,2:]=db;derivative[:,2:,1]=da
            ends.append((basis[site]@core@inverse[site],basis[site]@derivative@inverse[site]))
        words=[];jacs=[]
        for z,inside in zip(roots,interior):
            left=ident+(1/z-1)*ends[1][0];right=ident+(1/z-1)*ends[0][0]
            words.append(left@inside@right)
            jacs.append(np.concatenate((left@inside@((1/z-1)*ends[0][1]),((1/z-1)*ends[1][1])@inside@right)))
        coefficients={};derivatives={}
        for degree in degrees:
            coefficients[degree]=sum(w*z**(-degree) for w,z in zip(words,roots))/len(roots)
            derivatives[degree]=sum(j*z**(-degree) for j,z in zip(jacs,roots))/len(roots)
        residuals=[];jac=[]
        for degree in degrees:
            if degree not in target:
                residuals.append(coefficients[degree].ravel());jac.append(derivatives[degree].reshape(count,-1).T)
        for degree,rank in target.items():
            c,d=coefficients[degree],derivatives[degree]
            residuals.append(np.array([np.trace(c)-rank]));jac.append(np.trace(d,axis1=1,axis2=2).reshape(1,count))
            for other in target:
                c2,d2=coefficients[other],derivatives[other]
                residuals.append((c@c2-(c if degree==other else 0)).ravel())
                jac.append((d@c2+c@d2-(d if degree==other else 0)).reshape(count,-1).T)
        f=np.concatenate(residuals);j=np.vstack(jac)
        return np.r_[f.real,f.imag],np.block([[j.real,-j.imag],[j.imag,j.real]])
    return np.r_[initial.real,initial.imag],calculate


def search(args):
    start=time.monotonic();target={int(k):v for k,v in json.loads(args.target).items()}
    x0,calculate=model(args.n,target);rng=np.random.default_rng(args.seed)
    test=x0+rng.normal(0,.01,len(x0));_,j=calculate(test)
    for column in rng.choice(len(x0),12,replace=False):
        d=np.zeros(len(x0));d[column]=1e-6
        error=(calculate(test+d)[0]-calculate(test-d)[0])/2e-6-j[:,column]
        if np.max(abs(error))>2e-5:raise ValueError('The boundary Jacobian control fails')
    for n,t in ((4,{0:4}),(6,{0:4,1:2})):
        x,c=model(n,t)
        if np.max(abs(c(x)[0]))>1e-10:raise ValueError('The inverse cancellation control fails')
    records=[];best=None;cache_x=None;cache_value=None
    def cached(x):
        nonlocal cache_x,cache_value
        if cache_x is None or not np.array_equal(x,cache_x):cache_x=x.copy();cache_value=calculate(x)
        return cache_value
    for trial in range(args.trials):
        initial=x0 if trial==0 else x0+rng.normal(0,args.scale,len(x0))
        try:
            fit=least_squares(lambda x:cached(x)[0],initial,jac=lambda x:cached(x)[1],max_nfev=args.steps,ftol=1e-11,xtol=1e-11,gtol=1e-11)
            row={'trial':trial,'norm':float(np.linalg.norm(fit.fun)),'max_residual':float(max(abs(fit.fun))),'evaluations':fit.nfev}
            if best is None or row['norm']<best['norm']:best=row|{'parameters':fit.x.tolist()}
        except np.linalg.LinAlgError as error:row={'trial':trial,'status':'NUMERICAL_FAILURE','error':str(error)}
        records.append(row);print(json.dumps(row),flush=True)
        receipt={'n':args.n,'target':target,'boundary_signs':'--','seed':args.seed,'scale':args.scale,'steps':args.steps,'records':records,'best':best,'seconds':time.monotonic()-start,'jacobian_control':'PASS','inverse_controls':'PASS','new_kappa':False}
        args.output.write_text(json.dumps(receipt,indent=2)+'\n')
        if row.get('max_residual',1)>1e-9:continue
        break


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=10)
    p.add_argument('--target',default='{"-2":2,"0":1,"1":5,"2":1,"3":1}')
    p.add_argument('--trials',type=int,default=4);p.add_argument('--steps',type=int,default=400)
    p.add_argument('--seed',type=int,default=7331);p.add_argument('--scale',type=float,default=.2)
    p.add_argument('--output',type=Path,required=True)
    search(p.parse_args())
