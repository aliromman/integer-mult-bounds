"""Replay stored trace witnesses without a certificate search.

Only rational arithmetic is used. A changed witness must fail.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from trace_curvature import equations


def polynomial(system, certificate):
    out={}
    for term in certificate:
        prefix=tuple(term['left_word']);scalar=F(term['coefficient'])
        for word,value in system[term['equation']].items():
            key=prefix+word;out[key]=out.get(key,F(0))+scalar*value
    return {w:v for w,v in out.items() if v}


def commutator_remainder(poly):
    remaining=dict(poly);pairs=0
    for word,coefficient in poly.items():
        if not word:continue
        rotation=min(range(len(word)),key=lambda i:word[i:]+word[:i])
        if not rotation:continue
        a,b=word[:rotation],word[rotation:]
        remaining[a+b]=remaining.get(a+b,F(0))-coefficient
        remaining[b+a]=remaining.get(b+a,F(0))+coefficient
        pairs+=1
    return {w:v for w,v in remaining.items() if v},pairs


def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def trace_value(poly,width):
    variables={v for word in poly for v in word}
    matrices={v:[[F((v+2*i-j)%7-3) for j in range(width)] for i in range(width)] for v in variables}
    total=F(0)
    for word,coefficient in poly.items():
        matrix=[[F(i==j) for j in range(width)] for i in range(width)]
        for v in word:matrix=multiply(matrix,matrices[v])
        total+=coefficient*sum(matrix[i][i] for i in range(width))
    return total


def verify():
    receipts=[]
    for n,signs in ((4,'++'),(6,'++'),(8,'++'),(10,'++'),(8,'--'),(10,'--')):
        _,system,_=equations(n,signs)
        suffix='-minusminus' if signs=='--' else ''
        file=Path(__file__).with_name('discovery-results')/f'trace-n{n}{suffix}-d1.json'
        certificate=json.loads(file.read_text())['certificate']
        p=polynomial(system,certificate)
        remainder,count=commutator_remainder(p)
        if remainder!={():F(1)}:raise ValueError('The stored trace witness does not equal one modulo commutators')
        changed=[dict(t) for t in certificate]
        changed[0]['coefficient']=str(F(changed[0]['coefficient'])+1)
        if commutator_remainder(polynomial(system,changed))[0]=={():F(1)}:raise ValueError('The changed coefficient was not rejected')
        for width in (2,3):
            if trace_value(p,width)!=width:raise ValueError('The literal matrix trace check fails')
        receipts.append({'n':n,'boundary_signs':signs,'terms':len(certificate),'commutators':count,'width_controls':[2,3],'changed_coefficient_rejected':True})
    return {'certificates':receipts,'new_kappa':False}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
    print('PASS_ALL_LENGTH_BOUNDARY_TRACE_CERTIFICATES')
