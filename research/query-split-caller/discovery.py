"""Split one canonical query into two full calls.

The split adds one query to the shortest word. Exact necessary equations
are tested over the free associative algebra. A survivor is not a caller.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
import exact_replay as e
import free_replay as nc


def construct(n,site,weight):
    us,vs,_=e.canonical(n)
    pairs=[(u,v,F(1)) for u,v in zip(us,vs)]
    u,v,_=pairs[site];pairs[site:site+1]=[(u,v,weight),(u,v,1-weight)]
    names=[];ps=[];ts=[]
    def variable(name):
        i=len(names);names.append(name);return {(i,):F(1)}
    for call,(u,v,scale) in enumerate(pairs):
        ts.append(e.numeric([[int(i==j)-2*scale*u[i]*v[j] for j in range(n)] for i in range(n)]))
        if call in (0,len(pairs)-1):ps.append(e.zero(n));continue
        pivot=next(i for i in range(n) if v[i]);drop=next(i for i in range(n) if u[i] and i!=pivot)
        h=[F(int(i==pivot),v[pivot]) for i in range(n)]
        columns=[u,h]+[[int(i==j)-v[j]*h[i] for i in range(n)] for j in range(n) if j not in (pivot,drop)]
        basis=list(map(list,zip(*columns)))
        a=[variable(f'a{call}_{j}') for j in range(n-2)]
        b=[variable(f'b{call}_{j}') for j in range(n-2)]
        c=variable(f'c{call}');core=e.zero(n);core[0][0]=core[1][1]=e.const(1);core[0][1]=c
        for j in range(n-2):core[j+2][1]=a[j];core[0][j+2]=b[j]
        ps.append(nc.mm(nc.mm(e.numeric(basis),core),e.numeric(e.inverse(basis))))
    return names,ps,ts


def coefficient(ps,ts,degree):
    n=len(ps[0]);plus=[e.eye(n)];minus=[e.eye(n)]
    for p,t in zip(ps,ts):
        for word,query in ((plus,e.eye(n)),(minus,t)):
            old=[nc.mm(query,x) for x in word]
            word[:]=[old[0]]+[e.ma(old[j],nc.mm(p,old[j-1])) for j in range(1,len(old))]
            if len(old)<=degree:word.append(nc.mm(p,old[-1]))
    return [x for row in e.ma(plus[degree],minus[degree]) for x in row]


def check(n,site,weight,full=False):
    start=time.monotonic();names,ps,ts=construct(n,site,weight)
    if any(coefficient(ps,ts,0)):raise ValueError('The split query word is not minus identity')
    values=e.linear_solution(coefficient(ps,ts,1));steps=[]
    if values is None:status='EXACT_FIRST_COEFFICIENT_CONTRADICTION'
    else:
        ps=[nc.ms(p,values) for p in ps]
        system=[x for p in ps for row in e.ma(nc.mm(p,p),p,-1) for x in row]+coefficient(ps,ts,2)
        while True:
            reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
            row={'variables':len({v for p in system for w in p for v in w}),'equations':sum(bool(p) for p in system),'rank':len(reduced),'affine':len(affine)}
            steps.append(row);print(json.dumps(row),flush=True)
            if any(set(p)=={()} for p in affine):status='EXACT_FREE_ALGEBRA_CONTRADICTION';break
            if not affine:status='UNRESOLVED_SYSTEM' if any(system) else 'SECOND_COEFFICIENT_CONSISTENT';break
            values=e.linear_solution(affine);system=[nc.substitute(p,values) for p in system];ps=[nc.ms(p,values) for p in ps]
    if full and ps is not None and 'CONTRADICTION' not in status:
        for degree in range(3,n+1):
            system+=coefficient(ps,ts,degree)
            while True:
                reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
                row={'degree':degree,'variables':len({v for p in system for w in p for v in w}),'rank':len(reduced),'affine':len(affine)}
                steps.append(row);print(json.dumps(row),flush=True)
                if any(set(p)=={()} for p in affine):status='EXACT_HIGHER_COEFFICIENT_CONTRADICTION';break
                if not affine:
                    status='FULL_SYSTEM_UNRESOLVED' if any(system) else 'FULL_QUERY_FAMILY';break
                values=e.linear_solution(affine);system=[nc.substitute(p,values) for p in system];ps=[nc.ms(p,values) for p in ps]
            if 'CONTRADICTION' in status:break
    result={'n':n,'split_site':site,'split_weight':str(weight),'queries':n+2,'full_calls':n+2,'boundary_stages':'cancelled from necessary query equations','initial_variables':len(names),'steps':steps,'status':status,'seconds':time.monotonic()-start,'new_kappa':False}
    if full and ps is not None and 'CONTRADICTION' not in status:
        candidate=[[[p.get((),F(0)) for p in row] for row in mat] for mat in ps]
        numeric=[e.numeric(mat) for mat in candidate]
        projector_ok=not any(x for p in numeric for row in e.ma(nc.mm(p,p),p,-1) for x in row)
        query_ok=not any(x for degree in range(n+1) for x in coefficient(numeric,ts,degree))
        result.update(constant_choice_projectors=projector_ok,constant_choice_full_query=query_ok,constant_choice=[[[str(x) for x in row] for row in mat] for mat in candidate])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=4);p.add_argument('--site',type=int,default=2);p.add_argument('--weight',type=F,default=F(1,2));p.add_argument('--output',type=Path,required=True)
    p.add_argument('--full',action='store_true');a=p.parse_args()
    if not 1<=a.site<a.n or a.weight in (0,1):raise ValueError('The split must be a nonzero interior pair')
    result=check(a.n,a.site,a.weight,a.full);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='constant_choice'}))
