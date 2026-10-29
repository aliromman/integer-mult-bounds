"""Check saved rational projector words by literal Laurent multiplication."""
from fractions import Fraction as F
import json
from pathlib import Path


def add(a,b):return [[x+y for x,y in zip(row,other)] for row,other in zip(a,b)]
def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def scale(a,c):return [[c*x for x in row] for row in a]


def word(ps,ts,signs):
    n=len(ps[0]);ident=[[F(i==j) for j in range(n)] for i in range(n)];zero=scale(ident,0)
    out={0:ident}
    for p,t,sign in zip(ps,ts,signs):
        new={}
        for degree,a in out.items():
            value=mul(t,a)
            for key,term in ((degree,mul(add(ident,scale(p,-1)),value)),(degree+sign,mul(p,value))):new[key]=add(new.get(key,zero),term)
        out={k:v for k,v in new.items() if v!=zero}
    return out


def verify():
    results=[]
    for n,file in ((6,'n6-alternating.json'),(8,'n8-alternating-chart.json')):
        data=json.loads((Path(__file__).with_name('discovery-results')/file).read_text())
        ident=[[F(i==j) for j in range(n)] for i in range(n)];zero=scale(ident,0)
        us=[[F(i==j) for i in range(n)] for j in range(n)]+[[F((-1)**(j+1)) for j in range(n)]]
        vs=[[F(0 if i==j else (1 if i>j else -1)) for j in range(n)] for i in range(n+1)]
        ts=[[[ident[i][j]-2*u[i]*v[j] for j in range(n)] for i in range(n)] for u,v in zip(us,vs)]
        ps=[zero]+[[[F(x) for x in row] for row in p] for p in data['constant_choice']]+[zero]
        signs=[1]+[-1 if i in data['inverse_interior_sites'] else 1 for i in range(1,n)]+[1]
        for p,t in zip(ps,ts):
            if mul(p,p)!=p or mul(p,t)!=mul(t,p):raise ValueError('A projector or compatibility equation fails')
        if any(sum(p[i][i] for i in range(n))!=2 for p in ps[1:-1]):raise ValueError('An interior rank differs from two')
        plus=word(ps,[ident]*(n+1),signs);minus=word(ps,ts,signs)
        if minus!={k:scale(a,-1) for k,a in plus.items()}:raise ValueError('The full Laurent query identity fails')
        changed=list(ts);changed[n//2]=ident
        if word(ps,changed,signs)==minus:raise ValueError('An omitted query was accepted')
        h=zero
        for degree,a in plus.items():h=add(h,scale(a,degree))
        if not any(mul(h,a)!=scale(a,degree) for degree,a in plus.items()):raise ValueError('The saved mixed-stage control unexpectedly has a monomial endpoint')
        results.append({'n':n,'query_identity':'PASS','omitted_query_rejected':True,'stage_endpoint':'NOT_A_CONSTANT_MONOMIAL_REPRESENTATION','stage_powers':sorted(plus)})
    print(json.dumps(results,indent=2));print('PASS_SIGNED_QUERY_WORDS_ONLY')


if __name__=='__main__':verify()
