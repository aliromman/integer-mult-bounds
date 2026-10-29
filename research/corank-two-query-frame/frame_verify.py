"""Check the exact Gram and asymmetric query constructions."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
spec=importlib.util.spec_from_file_location('corank_frame',Path(__file__).with_name('discovery.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def product(a,b):
    out=[[F(0) for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for k in range(len(b)):
            for j in range(len(b[0])):out[i][j]+=a[i][k]*b[k][j]
    return out


def rank(a):
    a=[list(map(F,row)) for row in a];r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];d=a[r][col];a[r]=[x/d for x in a[r]]
        for i in range(r+1,len(a)):
            c=a[i][col]
            if c:a[i]=[x-c*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


def pfaffian(a):
    @lru_cache(None)
    def rec(indices):
        if not indices:return F(1)
        return sum((-1)**(j+1)*a[indices[0]][indices[j]]*rec(indices[1:j]+indices[j+1:]) for j in range(1,len(indices)))
    return rec(tuple(range(len(a))))


def check_word(n,pairs,queries):
    ident=[[F(i==j) for j in range(n)] for i in range(n)];out=ident
    for (u,v),query in zip(pairs,queries):
        r=[[u[i]*v[j] for j in range(n)] for i in range(n)]
        if rank(r)!=1 or any(x for row in product(r,r) for x in row):raise ValueError('A query is not rank-one square zero')
        out=product(query,out)
    if out!=[[-x for x in row] for row in ident]:raise ValueError('The literal query endpoint differs')
    u=[list(row) for row in zip(*(a for a,b in pairs))];v=[b for a,b in pairs];g=product(v,u);length=len(pairs)
    z=[[F(i==j)+(F(-1,2) if i>j else F(1,2))*g[i][j] for j in range(length)] for i in range(length)]
    if rank(g)!=n or rank(z)!=2:raise ValueError('The Gram or complement rank differs')


def verify():
    rational=[];expanded=0
    for n in (4,6,8,12,16):
        for mode in ('integer','dyadic'):
            t,pairs,queries=m.frame(n,mode);check_word(n,pairs,queries)
            if n<=8:
                h=[[F(0) if i==j else (1 if i<j else -1)*(t[i]+t[j]) for j in range(n+2)] for i in range(n+2)]
                if pfaffian(h):raise ValueError('The prescribed Pfaffian does not vanish')
                t[-1]*=2
                h=[[F(0) if i==j else (1 if i<j else -1)*(t[i]+t[j]) for j in range(n+2)] for i in range(n+2)]
                expected=F(2)**(n//2)*(m.math.prod(t[::2])+m.math.prod(t[1::2]))
                if pfaffian(h)!=expected or rank(h)!=n+2:raise ValueError('The changed-product control fails')
            rational.append({'n':n,'mode':mode,'query_endpoint':'PASS','complement_rank':2})
    for n in (4,6,8):
        for site in range(n):
            for kernel in range(n-2):
                _,pairs,queries=m.frame(n,f'expand-{site}-{kernel}');check_word(n,pairs,queries)
                triple=pairs[site:site+3]
                if rank([list(row) for row in zip(*(u for u,v in triple))])!=3 or rank([v for u,v in triple])!=2:raise ValueError('The asymmetric image and covector ranks differ')
                expanded+=1
    _,pairs,queries=m.frame(4,'triangles');check_word(4,pairs,queries)
    result={'rational_frames':rational,'asymmetric_frames':expanded,'two_triangle_control':'PASS','changed_parameter_rejected':True,'new_kappa':False}
    print(json.dumps(result,indent=2));print('PASS_CORANK_TWO_QUERY_FRAMES_ONLY')


if __name__=='__main__':verify()
