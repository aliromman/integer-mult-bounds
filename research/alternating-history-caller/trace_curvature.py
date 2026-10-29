"""Search exact cyclic trace certificates for boundary curvature.

The trace of a commutator is zero at every finite matrix width. A unit
in the cyclic row span therefore excludes every positive finite width.
A failed finite search makes no impossibility claim.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import time
from boundary_curvature import e, nc, construct, fixed


def equations(n, signs='++'):
    names, ps, _, _ = construct(n, 0, tuple(range(1,n)))
    ps = [ps[0]]+fixed(n)+[ps[-1]]
    h, second = e.zero(n), e.zero(n)
    for site,p in enumerate(ps):
        sign=-1 if (site==0 and signs[0]=='-') or (site==n and signs[1]=='-') else 1
        weighted=[[e.scale(x,sign) for x in row] for row in p]
        second=e.ma(second,nc.mm(weighted,h));h=e.ma(h,weighted)
    c=e.ma(e.ma(nc.mm(h,h),h,-1),[[e.scale(p,2) for p in row] for row in second],-1)
    for site,sign in ((0,signs[0]),(-1,signs[1])):
        if sign=='-': c=e.ma(c,[[e.scale(x,2) for x in row] for row in ps[site]],-1)
    system=[p for row in c for p in row]
    system += [p for end in (ps[0],ps[-1]) for row in e.ma(nc.mm(end,end),end,-1) for p in row]
    steps=[]
    while True:
        reduced=e.row_reduce(system)
        affine=[p for p in reduced if all(len(m)<=1 for m in p)]
        steps.append({'variables':len({v for p in system for m in p for v in m}), 'rank':len(reduced), 'affine':len(affine)})
        if any(set(p)=={()} for p in affine): return names,[{():F(1)}],steps
        if not affine: return names,[p for p in system if p],steps
        values=e.linear_solution(affine)
        system=[nc.substitute(p,values) for p in system]


def cyclic(poly):
    out={}
    for word,value in poly.items():
        key=min(word[i:]+word[:i] for i in range(len(word))) if word else ()
        out[key]=out.get(key,F(0))+value
    return {m:v for m,v in out.items() if v}


def subtract(left,right,scalar):
    for word,value in right.items():
        left[word]=left.get(word,F(0))-scalar*value
        if not left[word]: del left[word]


def search(n, degree, signs='++'):
    start=time.monotonic()
    names,system,steps=equations(n,signs)
    variables=sorted({v for p in system for m in p for v in m})
    pivots={}; tested=0
    for d in range(degree+1):
        for word in itertools.product(variables,repeat=d):
            for index,p in enumerate(system):
                row=cyclic(nc.mul({word:F(1)},p)); witness={(index,word):F(1)}
                tested+=1
                while row:
                    pivot=max(row,key=lambda m:(len(m),m)); value=row[pivot]
                    if pivot not in pivots:
                        row={m:v/value for m,v in row.items()}
                        witness={m:v/value for m,v in witness.items()}
                        if pivot==():
                            check={}
                            for (i,w),coefficient in witness.items():
                                term=cyclic(nc.mul({w:coefficient},system[i]))
                                subtract(check,term,F(-1))
                            if check!={():F(1)}: raise ValueError('The trace certificate does not equal one')
                            certificate=[{'equation':i,'left_word':list(w),'coefficient':str(c)} for (i,w),c in sorted(witness.items())]
                            return {'n':n,'boundary_signs':signs,'status':'EXACT_FINITE_WIDTH_TRACE_CONTRADICTION','multiplier_degree':d,'rows_tested':tested,'rank':len(pivots)+1,'steps':steps,'certificate':certificate,'seconds':time.monotonic()-start,'new_kappa':False}
                        pivots[pivot]=(row,witness);break
                    old,proof=pivots[pivot]
                    subtract(row,old,value);subtract(witness,proof,value)
        print(json.dumps({'degree':d,'rows':tested,'rank':len(pivots)}),flush=True)
    return {'n':n,'boundary_signs':signs,'status':'NO_CERTIFICATE_IN_TESTED_TRACE_SPAN','multiplier_degree':degree,'rows_tested':tested,'rank':len(pivots),'steps':steps,'seconds':time.monotonic()-start,'new_kappa':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n',type=int,default=4)
    p.add_argument('--degree',type=int,default=1)
    p.add_argument('--signs',choices=('++','+-','-+','--'),default='++')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=search(a.n,a.degree,a.signs)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='certificate'}))
