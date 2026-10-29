"""Test signed full-stage callers with the shortest query frame.

Multiplication by one scalar z per inverse stage replaces its projector
by the complementary projector. This gives exact polynomial equations.
The complete original rank-two projector chart is retained.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
import exact_replay as e
import free_replay as nc


def stage(ps,inverse):
    return [e.ma(e.eye(len(p)),p,-1) if i+1 in inverse else p for i,p in enumerate(ps)]


def search(n,inverse,maxdegree,seconds):
    start=time.monotonic();names,ps,_,ts=e.chart(n,())
    steps=[];status='UNRESOLVED';system=[];highest=0
    def expire(signum,frame):raise TimeoutError('The exact elimination time limit expired')
    signal.signal(signal.SIGALRM,expire);signal.alarm(seconds)
    try:
        values=e.linear_solution(nc.coefficient(stage(ps,inverse),ts,1))
        if values is None:status='EXACT_FIRST_COEFFICIENT_CONTRADICTION'
        else:
            ps=[nc.ms(p,values) for p in ps]
            system=[p for mat in ps for row in e.ma(nc.mm(mat,mat),mat,-1) for p in row]
            for degree in range(2,maxdegree+1):
                highest=degree
                system+=nc.coefficient(stage(ps,inverse),ts,degree)
                while True:
                    reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
                    row={'degree':degree,'variables':len({v for p in system for w in p for v in w}),'rank':len(reduced),'affine':len(affine)}
                    steps.append(row);print(json.dumps(row),flush=True)
                    if any(set(p)=={()} for p in affine):status='EXACT_FREE_ALGEBRA_CONTRADICTION';break
                    if not affine:
                        status='UNRESOLVED' if any(system) else 'COEFFICIENT_FAMILY';break
                    values=e.linear_solution(affine);system=[nc.substitute(p,values) for p in system];ps=[nc.ms(p,values) for p in ps]
                if 'CONTRADICTION' in status:break
    except TimeoutError as error:status='UNKNOWN_TIMEOUT'
    finally:signal.alarm(0)
    result={'n':n,'inverse_interior_sites':list(inverse),'requested_degree':maxdegree,'last_degree':highest,'status':status,'steps':steps,'limit_seconds':seconds,'seconds':time.monotonic()-start,'new_kappa':False}
    if 'CONTRADICTION' not in status and status!='UNKNOWN_TIMEOUT':
        raw=[[[p.get((),F(0)) for p in row] for row in mat] for mat in ps];numeric=[e.numeric(mat) for mat in raw]
        result['constant_choice_projectors']=not any(p for mat in numeric for row in e.ma(nc.mm(mat,mat),mat,-1) for p in row)
        result['constant_choice_query']=not any(p for degree in range(1,n) for p in nc.coefficient(stage(numeric,inverse),ts,degree))
        result['constant_choice']=[[[str(x) for x in row] for row in mat] for mat in raw]
        variables=sorted({v for mat in ps for row in mat for p in row for w in p for v in w})
        result['chart_variables']=variables
        result['chart']=[[[[[list(w),str(v)] for w,v in sorted(p.items())] for p in row] for row in mat] for mat in ps]
        result['equations']=[[[list(w),str(v)] for w,v in sorted(p.items())] for p in system if p]
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=6);p.add_argument('--inverse',type=int,nargs='*',default=[1,3,5]);p.add_argument('--degree',type=int);p.add_argument('--seconds',type=int,default=180);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.n<2 or a.n%2 or any(not 1<=i<a.n for i in a.inverse):raise ValueError('Use even size and interior inverse sites')
    r=search(a.n,tuple(a.inverse),a.degree or a.n-1,a.seconds);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ('constant_choice','chart','equations')}))
