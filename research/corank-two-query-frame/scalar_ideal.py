"""Test the exact scalar second-order ideal of a named query frame."""
import argparse
import importlib.util
import json
from pathlib import Path
import signal
import time
import sympy as s
spec=importlib.util.spec_from_file_location('corank_frame',Path(__file__).with_name('discovery.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);e=m.e


def check(n,mode,seconds):
    start=time.monotonic();_,names,ps,ts=m.construct(n,mode)
    values=e.linear_solution(m.coefficient(ps,ts,1))
    if values is None:return {'n':n,'mode':mode,'status':'EXACT_FIRST_COEFFICIENT_CONTRADICTION','new_kappa':False}
    ps=[e.ms(p,values) for p in ps]
    system=[p for mat in ps for row in e.ma(e.mm(mat,mat),mat,-1) for p in row]+m.coefficient(ps,ts,2)
    while True:
        reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
        if any(set(p)=={()} for p in affine):return {'n':n,'mode':mode,'status':'EXACT_AFFINE_CONTRADICTION','new_kappa':False}
        if not affine:break
        values=e.linear_solution(affine);system=[e.substitute(p,values) for p in system]
    indices=sorted({v for p in system for w in p for v in w});variables=s.symbols('x:'+str(len(indices)));lookup=dict(zip(indices,variables))
    expressions=list(dict.fromkeys(sum(s.Rational(c.numerator,c.denominator)*s.prod(lookup[v] for v in w) for w,c in p.items()) for p in system if p))
    def expire(signum,frame):raise TimeoutError('The Groebner time limit expired')
    signal.signal(signal.SIGALRM,expire);signal.alarm(seconds)
    try:
        basis=s.groebner(expressions,*variables,order='grevlex',method='f5b')
        result={'status':'EXACT_SCALAR_UNIT_IDEAL' if list(basis)==[s.Integer(1)] else 'NONUNIT_GROEBNER_BASIS','basis':[str(p) for p in basis]}
    except TimeoutError:result={'status':'UNKNOWN_GROEBNER_TIMEOUT'}
    finally:signal.alarm(0)
    return {'n':n,'mode':mode,'variables':len(variables),'equations':len(expressions),'limit_seconds':seconds,'seconds':time.monotonic()-start,'sympy':s.__version__,'scope':'scalar commutative second-order equations','new_kappa':False,**result}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=6);p.add_argument('--mode',default='integer');p.add_argument('--seconds',type=int,default=90);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=check(a.n,a.mode,a.seconds);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
