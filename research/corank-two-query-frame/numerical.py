"""Search full Laurent endpoints in a saved fixed-query projector chart.

The residual checks every query coefficient, every projector identity,
and the complete one-parameter stage identity. Results need exact replay.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
from exact_replay import canonical


def model(data,ends,profile):
    n=data['n'];vs=data['chart_variables'];lookup={v:i for i,v in enumerate(vs)}
    count=len(vs)+4*(n-2);ident=np.eye(n,dtype=complex)
    interiors=len(data['chart'])
    constant=np.zeros((interiors,n,n),dtype=complex);linear=np.zeros((interiors,count,n,n),dtype=complex)
    for k,mat in enumerate(data['chart']):
        for i,row in enumerate(mat):
            for j,terms in enumerate(row):
                for word,c in terms:
                    from fractions import Fraction
                    c=float(Fraction(c))
                    if not word:constant[k,i,j]+=c
                    elif len(word)==1:linear[k,lookup[word[0]],i,j]+=c
                    else:raise ValueError('The saved chart is not affine')
    us,rows,queries=canonical(n)
    if 'raw_queries' in data:
        from fractions import Fraction
        number=lambda x:float(Fraction(x))
        us=[[number(x) for x in row] for row in data['query_u']]
        rows=[[number(x) for x in row] for row in data['query_v']]
        queries=[[[number(x) for x in row] for row in mat] for mat in data['raw_queries']]
    bases=[]
    for site in (0,len(queries)-1):
        u,v=np.array(us[site],float),np.array(rows[site],float)
        pivot=next(j for j,x in enumerate(v) if x);drop=next(j for j,x in enumerate(u) if x and j!=pivot)
        h=np.eye(n)[:,pivot]/v[pivot]
        basis=np.column_stack([u,h,*[np.eye(n)[:,j]-v[j]*h for j in range(n) if j not in (pivot,drop)]])
        bases.append((basis,np.linalg.inv(basis)))
    signs=[ends[0]]+[-1 if i in data['inverse_interior_sites'] else 1 for i in range(1,interiors+1)]+[ends[1]]
    ts=[np.asarray(t,complex) for t in queries]
    zero=np.zeros((n,n),complex);dzero=np.zeros((count,n,n),complex)
    def boundary(params,offset,bi):
        d=n-2;a=params[offset:offset+d];b=params[offset+d:offset+2*d]
        p=np.zeros((n,n),complex);p[0,0]=p[1,1]=1;p[0,1]=-b@a;p[2:,1]=a;p[0,2:]=b
        dp=np.zeros((count,n,n),complex)
        for j in range(d):
            dp[offset+j,j+2,1]=dp[offset+d+j,0,j+2]=1
            dp[offset+j,0,1]=-b[j];dp[offset+d+j,0,1]=-a[j]
        basis,inv=bases[bi];return basis@p@inv,basis@dp@inv
    def calculate(x):
        params=x[:count]+1j*x[count:]
        p0,dp0=boundary(params,len(vs),0);pn,dpn=boundary(params,len(vs)+2*(n-2),1)
        ps=[p0]+list(constant+np.einsum('kvij,v->kij',linear,params))+[pn]
        dps=[dp0]+list(linear)+[dpn]
        channels=[]
        for minus in (False,True):
            word={0:ident};dw={0:dzero}
            for p,dp,t,sign in zip(ps,dps,ts,signs):
                query=t if minus else ident;old={k:query@a for k,a in word.items()};od={k:query@a for k,a in dw.items()};word={};dw={}
                for k,a in old.items():
                    for shift,op,dop in [(0,ident-p,-dp),(sign,p,dp)]:
                        key=k+shift;word[key]=word.get(key,zero)+op@a;dw[key]=dw.get(key,dzero)+op@od[k]+dop@a
            channels.append((word,dw))
        plus,der=channels[0];minus,md=channels[1]
        h=sum(s*p for s,p in zip(signs,ps));dh=sum(s*p for s,p in zip(signs,dps))
        val=[];jac=[]
        def append(a,da):val.append(a.ravel());jac.append(da.reshape(count,a.size).T)
        for k in plus:
            a,da=plus[k],der[k]
            append((a+minus[k])/10,(da+md[k])/10)
            append((h@a-k*a)/10,(dh@a+h@da-k*da)/10)
            append(np.array([np.trace(a)-profile.get(k,0)]),np.trace(da,axis1=1,axis2=2)[:,None])
        for p,dp in zip(ps,dps):append(p@p-p,dp@p+p@dp-dp)
        f=np.concatenate(val);j=np.vstack(jac)
        return np.r_[f.real,f.imag],np.block([[j.real,-j.imag],[j.imag,j.real]])
    return count,calculate


def search(args):
    data=json.loads(args.chart.read_text());profile={int(k):int(v) for k,v in json.loads(args.profile).items()};ends=tuple(map(int,args.ends.split(',')))
    if sum(profile.values())!=data['n']:raise ValueError('The target multiplicities do not match the input dimension')
    count,calculate=model(data,ends,profile);rng=np.random.default_rng(args.seed)
    control={'n':2,'chart_variables':[],'inverse_interior_sites':[],
             'chart':[[[[[[], '1']] if i==j else [] for j in range(2)] for i in range(2)]]}
    _,positive=model(control,(1,1),{3:2});_,negative=model(control,(1,1),{2:2})
    if np.any(positive(np.zeros(0))[0]) or not np.any(negative(np.zeros(0))[0]):
        raise ValueError('The exact three-query control or its changed target fails')
    probe=rng.normal(0,.02,2*count);_,jac=calculate(probe)
    for i in np.linspace(0,2*count-1,min(24,2*count),dtype=int):
        delta=np.zeros(2*count);delta[i]=1e-6;fd=(calculate(probe+delta)[0]-calculate(probe-delta)[0])/2e-6
        if np.max(abs(fd-jac[:,i]))>1e-4:raise ValueError('The analytic derivative control fails')
    records=[];best=None;start=time.monotonic()
    for trial in range(args.trials):
        x=np.zeros(2*count) if trial==0 else rng.normal(0,.1,2*count)
        result=least_squares(lambda x:calculate(x)[0],x,jac=lambda x:calculate(x)[1],max_nfev=args.steps,ftol=1e-10,xtol=1e-10,gtol=1e-10)
        row={'trial':trial,'norm':float(np.linalg.norm(result.fun)),'max_residual':float(max(abs(result.fun))),'evaluations':result.nfev};records.append(row);print(json.dumps(row),flush=True)
        if best is None or row['norm']<best['norm']:best=row|{'parameters':result.x.tolist()}
        if row['max_residual']<1e-10:break
    return {'n':data['n'],'source_chart':args.chart.name,'boundary_signs':ends,'target_profile':profile,'records':records,'best':best,'seconds':time.monotonic()-start,'derivative_control':'PASS_24_COLUMNS','exact_triangle_control':'PASS','changed_target_rejected':True,'status':'NUMERICAL_DISCOVERY_ONLY','new_kappa':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--chart',type=Path,required=True);p.add_argument('--ends',default='1,1');p.add_argument('--profile',default='{"2":3,"-2":1,"-3":2}');p.add_argument('--trials',type=int,default=6);p.add_argument('--steps',type=int,default=800);p.add_argument('--seed',type=int,default=290510);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=search(a);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ('records','best')}))
