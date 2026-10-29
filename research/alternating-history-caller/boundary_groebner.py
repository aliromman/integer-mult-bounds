"""Check scalar boundary curvature with a bounded exact Groebner calculation.

An unresolved basis or timeout is not a no-caller proof.
"""
import argparse
import json
from pathlib import Path
import signal
import time
import sympy as s
from boundary_curvature import e, construct, fixed


def check(n, limit):
    start = time.monotonic()
    names, ps, _, _ = construct(n, 0, tuple(range(1,n)))
    ps = [ps[0]]+fixed(n)+[ps[-1]]
    h, second = e.zero(n), e.zero(n)
    for p in ps:
        second = e.ma(second,e.mm(p,h))
        h=e.ma(h,p)
    curvature=e.ma(e.ma(e.mm(h,h),h,-1),[[e.scale(p,2) for p in row] for row in second],-1)
    system=[p for row in curvature for p in row]
    system += [p for end in (ps[0],ps[-1]) for row in e.ma(e.mm(end,end),end,-1) for p in row]
    steps=[]
    while True:
        reduced=e.row_reduce(system)
        affine=[p for p in reduced if all(len(m)<=1 for m in p)]
        steps.append({'variables':len({v for p in system for m in p for v in m}), 'rank':len(reduced), 'affine':len(affine)})
        if any(set(p)=={()} for p in affine):
            return {'n':n,'status':'EXACT_AFFINE_CURVATURE_CONTRADICTION','steps':steps,'seconds':time.monotonic()-start,'new_kappa':False}
        if not affine: break
        values=e.linear_solution(affine)
        system=[e.substitute(p,values) for p in system]
    used=sorted({v for p in system for m in p for v in m})
    symbols={i:s.Symbol(names[i]) for i in used}
    polys=[sum(s.Rational(c.numerator,c.denominator)*s.prod(symbols[i] for i in word) for word,c in p.items()) for p in system if p]
    def expired(signum,frame): raise TimeoutError('The exact Groebner time limit expired')
    signal.signal(signal.SIGALRM,expired)
    signal.alarm(limit)
    try:
        basis=s.groebner(polys,*symbols.values(),order='grevlex')
        answer={'status':'EXACT_GROEBNER_CURVATURE_CONTRADICTION' if list(basis)==[1] else 'UNRESOLVED_CURVATURE_BASIS', 'basis':[str(p) for p in basis]}
    except TimeoutError as error:
        answer={'status':'UNKNOWN_TIMEOUT','error':str(error)}
    finally: signal.alarm(0)
    return {'n':n,'steps':steps,'limit_seconds':limit,'seconds':time.monotonic()-start,'new_kappa':False,**answer}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=6)
    p.add_argument('--seconds',type=int,default=30)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); result=check(a.n,a.seconds)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='basis'}))
