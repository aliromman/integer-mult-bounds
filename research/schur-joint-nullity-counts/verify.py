"""Compare incidence inversion with complete matrix and factor domains.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
from collections import Counter
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import types

P=Path(__file__).resolve().parent
ROOT=P.parents[1]


def need(ok,message):
    if not ok:raise ValueError(message)


def load(name,path):
    m=types.ModuleType(name);m.__file__=str(path);sys.modules[name]=m
    exec(compile(path.read_bytes(),str(path),'exec'),m.__dict__)
    return m


class Field:
    def __init__(self,q):
        self.q=q
        self.poly={4:0b111,8:0b1011}.get(q)
        self.add=lambda a,b:(a^b) if self.poly else (a+b)%q
        self.neg=lambda a:a if self.poly else (-a)%q
        self.inv={a:next(b for b in range(1,q) if self.mul(a,b)==1) for a in range(1,q)}
    def mul(self,a,b):
        if not self.poly:return a*b%self.q
        out=0
        while b:
            if b&1:out^=a
            b>>=1;a<<=1
            if a&self.q:a^=self.poly
        return out
    def rank(self,a):
        a=[row[:] for row in a]
        if not a:return 0
        rank=0
        for j in range(len(a[0])):
            k=next((i for i in range(rank,len(a)) if a[i][j]),None)
            if k is None:continue
            a[k],a[rank]=a[rank],a[k]
            inv=self.inv[a[rank][j]];a[rank]=[self.mul(inv,v) for v in a[rank]]
            for i in range(rank+1,len(a)):
                c=a[i][j]
                a[i]=[self.add(x,self.neg(self.mul(c,y))) for x,y in zip(a[i],a[rank])]
            rank+=1
            if rank==len(a):break
        return rank
    def mm(self,x,y,d):
        out=[]
        for row in x:
            result=[]
            for j in range(d):
                value=0
                for k,a in enumerate(row):value=self.add(value,self.mul(a,y[k][j]))
                result.append(value)
            out.append(result)
        return out


def matrices(r,c,q):
    for v in product(range(q),repeat=r*c):yield [list(v[i*c:(i+1)*c]) for i in range(r)]


def fingerprint(value):
    raw=value.to_bytes(max(1,(value.bit_length()+7)//8),'big')
    return dict(bits=value.bit_length(),sha256=hashlib.sha256(raw).hexdigest())


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write-result',action='store_true');args=ap.parse_args()
    for row in json.loads((P/'SOURCES.json').read_text())['files']:
        if row['usage']=='input':need(hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],'source hash')
    c=load('counts',P/'counts.py')
    old=load('binary_source',ROOT/'research/schur-header-count/verify.py')
    domains=[];total=0
    for q,d in ((2,1),(2,2),(2,3),(3,1),(3,2),(3,3),(5,1),(5,2),(4,1),(4,2),(8,1),(8,2)):
        field=Field(q);observed=Counter()
        for h in matrices(d,d,q):
            w=[[field.add(int(i==j),field.neg(h[i][j])) for j in range(d)] for i in range(d)]
            observed[d-field.rank(h),d-field.rank(w)]+=1
        predicted=c.joint_nullities(q,d)
        need(dict(observed)=={k:v for k,v in predicted.items() if v},'complete joint-nullity census')
        domains.append(dict(q=q,d=d,matrices=sum(observed.values()),nonempty_strata=len(observed)))
        total+=sum(observed.values())
    factor_pairs=0;factor_targets=0
    for q,d,r in ((2,1,0),(2,1,1),(2,1,2),(2,2,2),(2,2,3),(3,1,2),(3,2,2),(4,1,2),(4,2,1)):
        field=Field(q);all_counts=Counter();full_counts=Counter()
        for x in matrices(d,r,q):
            for y in matrices(r,d,q):
                a=field.mm(x,y,d);key=tuple(map(tuple,a));all_counts[key]+=1;factor_pairs+=1
                if field.rank(x)==d and field.rank(y)==d:full_counts[key]+=1
        for a in matrices(d,d,q):
            key=tuple(map(tuple,a));rho=field.rank(a)
            need(all_counts[key]==c.unrestricted(q,r,d,rho),'unrestricted factor count')
            need(full_counts[key]==c.full_rank(q,r,d,rho),'full factor count')
            factor_targets+=1
    comparisons=[]
    for d,profile in ((1,(1,1,1,2)),(1,(0,2,2,1)),(2,(2,2,2,2)),(2,(1,3,3,2)),(3,(3,3,3,4))):
        got=c.headers(2,d,profile);want=old.predict(profile,d)
        need((got['sharp'],got['common'],got['total'])==(want['sharp_headers'],want['identity_headers'],want['total_headers']),'retained binary header law')
        comparisons.append(dict(d=d,profile=list(profile),sharp_probability=str(got['sharp_probability']),conditional_probability=str(got['conditional_probability'])))
    # Literal full headers use the predecessor's independent bit-matrix oracle.
    profile=(1,1,1,1);common=sharp=0
    for h in old.gl_headers(4):
        rank=old.inspect_header(h,profile,1)
        common+=rank is not None;sharp+=rank==2
    got=c.headers(2,1,profile)
    need((common,sharp)==(got['common'],got['sharp']),'all GL4 literal headers')
    large=[]
    for q,d in ((2,32),(4,16)):
        law=c.joint_nullities(q,d)
        need(law[d,0]==law[0,d]==1 and all(law[a,b]==law[b,a] for a,b in law),'large symmetry and unique endpoints')
        result=c.headers(q,d,(d,d,d,d))
        need(0<=result['sharp']<=result['common']<=result['total'],'large containment')
        large.append(dict(q=q,d=d,strata=len(law),matrix_space_bits=(q**(d*d)).bit_length(),sharp=fingerprint(result['sharp']),common=fingerprint(result['common']),total=fingerprint(result['total'])))
    need(c.joint_nullities(2,1).get((1,1),0)==0,'zero and one eigenspaces disjoint')
    need(c.headers(2,1,(1,0,0,0))['conditional_probability'] is None,'empty event has no conditional probability')
    need(c.headers(3,1,(1,1,1,1))['sharp']==16*c.gl(3,3),'odd-field scalar H2 count')
    need(c.headers(2,1,(0,2,2,1))['sharp']==181440,'unrestricted00 boundary')
    result=dict(status='PASS_SCHUR_JOINT_NULLITY_COUNTS',domains=domains,matrices=total,factor_pairs=factor_pairs,factor_targets=factor_targets,binary_comparisons=comparisons,literal_headers=c.gl(2,4),large=large,controls=4,counts_sha256=hashlib.sha256((P/'counts.py').read_bytes()).hexdigest())
    if args.write_result:(P/'RESULT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:need(json.dumps(result,sort_keys=True)==json.dumps(json.loads((P/'RESULT.json').read_text()),sort_keys=True),'typed complete receipt')
    print(json.dumps({k:result[k] for k in ('status','matrices','factor_pairs','factor_targets','literal_headers','controls')}))


if __name__=='__main__':main()
