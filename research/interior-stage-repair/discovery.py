"""Test one internal stage projector before boundary completion.

The query identity and the stage degree bound are linear in the inserted
matrix. The full rank and idempotence obligations remain separate.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'alternating-history-caller'))
from boundary_curvature import e,fixed,nc


def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def add(a,b,sign=1):return [[x+sign*y for x,y in zip(row,other)] for row,other in zip(a,b)]
def product(a,b):
    n=len(a[0]);out=[[[0]*n for _ in range(n)] for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=add(out[i+j],mm(x,y))
    return out


def setup(n):
    ident=[[int(i==j) for j in range(n)] for i in range(n)]
    zero=[[0]*n for _ in range(n)]
    ps=[zero]+[[[int(x.get((),0)) for x in row] for row in p] for p in fixed(n)]+[zero]
    ts=e.canonical(n)[2];channels=[]
    for minus in (False,True):
        factors=[]
        for p,t in zip(ps,ts):
            query=t if minus else ident
            factor=[mm(add(ident,p,-1),query)]
            if any(v for row in p for v in row):factor.append(mm(p,query))
            factors.append(factor)
        prefixes=[[ident]]
        for f in factors:prefixes.append(product(f,prefixes[-1]))
        suffixes=[None]*(len(factors)+1);suffixes[-1]=[ident]
        for i in range(len(factors)-1,-1,-1):suffixes[i]=product(suffixes[i+1],factors[i])
        channels.append((prefixes,suffixes))
    if any(add(a,b)!=zero for a,b in zip(channels[0][0][-1],channels[1][0][-1])):raise ValueError('The original full query identity fails')
    return channels


def linear_product(left,right,n):
    out=[[[{} for _ in range(n)] for _ in range(n)] for _ in range(len(left)+len(right)-1)]
    for a,l in enumerate(left):
        ls=[(i,k,x) for i,row in enumerate(l) for k,x in enumerate(row) if x]
        for b,r in enumerate(right):
            rs=[(h,j,y) for h,row in enumerate(r) for j,y in enumerate(row) if y]
            for i,k,x in ls:
                for h,j,y in rs:
                    key=(k*n+h,);p=out[a+b][i][j];p[key]=p.get(key,F(0))+x*y
    return [[[{m:v for m,v in p.items() if v} for p in row] for row in c] for c in out]


def check(n,cut,cap,seconds,inverse=False):
    start=time.monotonic();channels=setup(n);linear=[]
    for prefix,suffix in channels:linear.append(linear_product(suffix[cut+1],prefix[cut+1],n))
    equations=[e.add(a,b) for ca,cb in zip(*linear) for ra,rb in zip(ca,cb) for a,b in zip(ra,rb)]
    stage=channels[0][0][-1];zero=e.zero(n)
    for degree in range(cap+1+int(inverse),len(linear[0])+1):
        now=linear[0][degree] if degree<len(linear[0]) else zero
        before=linear[0][degree-1]
        old_degree=degree-int(inverse)
        old=stage[old_degree] if old_degree<len(stage) else [[0]*n for _ in range(n)]
        sign=-1 if inverse else 1
        equations += [e.add(e.const(c),e.scale(b,sign),e.scale(a,-sign)) for rr,rb,ra in zip(old,before,now) for c,b,a in zip(rr,rb,ra)]
    def expired(signum,frame):raise TimeoutError('The exact linear solve time limit expired')
    signal.signal(signal.SIGALRM,expired);signal.alarm(seconds)
    try:
        values=e.linear_solution(equations)
        if values is None:result={'status':'EXACT_LINEAR_CONTRADICTION_AT_THIS_CUT'}
        else:
            candidate=[[values.get(i*n+j,{(i*n+j,):F(1)}).get((),F(0)) for j in range(n)] for i in range(n)]
            remaining=sorted(set(range(n*n))-set(values))
            result={'status':'LINEAR_SYSTEM_CONSISTENT','free_variables':len(remaining),'constant_choice_idempotent':mm(candidate,candidate)==candidate,'constant_choice':[[str(x) for x in row] for row in candidate]}
            matrix=[[values.get(i*n+j,{(i*n+j,):F(1)}) for j in range(n)] for i in range(n)]
            system=[p for row in e.ma(nc.mm(matrix,matrix),matrix,-1) for p in row]
            steps=[]
            while True:
                reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
                steps.append({'rank':len(reduced),'affine':len(affine)})
                if any(set(p)=={()} for p in affine):
                    result['status']='EXACT_PROJECTOR_CONTRADICTION_AT_THIS_CUT';break
                if not affine:
                    result['status']='PROJECTOR_SYSTEM_UNRESOLVED' if any(system) else 'EXACT_PROJECTOR_FAMILY';break
                sub=e.linear_solution(affine);system=[nc.substitute(p,sub) for p in system];matrix=nc.ms(matrix,sub)
            result['projector_steps']=steps
            if result['status']=='EXACT_PROJECTOR_FAMILY':
                candidate=[[p.get((),F(0)) for p in row] for row in matrix]
                result['constant_choice']=[[str(x) for x in row] for row in candidate]
                result['constant_choice_idempotent']=mm(candidate,candidate)==candidate
    except TimeoutError as error:result={'status':'UNKNOWN_TIMEOUT','error':str(error)}
    finally:signal.alarm(0)
    return {'n':n,'cut_after_query':cut,'degree_cap':cap,'inverse_stage':inverse,'equations':len(equations),'limit_seconds':seconds,'seconds':time.monotonic()-start,'new_kappa':False,**result}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=8);p.add_argument('--cut',type=int,default=4);p.add_argument('--cap',type=int,default=2);p.add_argument('--seconds',type=int,default=120);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inverse',action='store_true')
    a=p.parse_args();result=check(a.n,a.cut,a.cap,a.seconds,a.inverse);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='constant_choice'}))
