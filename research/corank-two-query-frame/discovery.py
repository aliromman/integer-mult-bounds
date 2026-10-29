"""Test rational corank-two query frames and full-stage projectors.

The query frame is checked before any stage equation is used. Projector
elimination is scalar and commutative; it does not exclude matrix coupling.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import signal
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
import exact_replay as e


def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def expanded_frame(n,site,kernel):
    us,vs,_=e.canonical(n)
    pairs=[(list(map(F,u)),[-2*F(x) for x in v]) for u,v in zip(us,vs)]
    ua,va=pairs[site];ub,vb=pairs[site+1]
    a,b=next((a,b) for a in range(n) for b in range(a+1,n) if va[a]*vb[b]-va[b]*vb[a])
    free=[j for j in range(n) if j not in (a,b)]
    if not 0<=kernel<len(free):raise ValueError('The kernel index is outside the basis')
    j=free[kernel];correction=mm(e.inverse([[va[a],va[b]],[vb[a],vb[b]]]),[[va[j]],[vb[j]]])
    w=[F(i==j) for i in range(n)];w[a]-=correction[0][0];w[b]-=correction[1][0]
    if sum(x*y for x,y in zip(va,w)) or sum(x*y for x,y in zip(vb,w)):raise ValueError('The expansion vector is outside the common kernel')
    ident=[[F(i==j) for j in range(n)] for i in range(n)]
    inverse_a=[[ident[i][j]-ua[i]*va[j] for j in range(n)] for i in range(n)]
    inverse_b=[[ident[i][j]-ub[i]*vb[j] for j in range(n)] for i in range(n)]
    coupling=sum(x*y for x,y in zip(vb,ua))
    residue=[(1+coupling)*x+y for x,y in zip(va,vb)]
    cancel=[-x for x in mm([residue],mm(inverse_a,inverse_b))[0]]
    pairs[site:site+2]=[([x+y for x,y in zip(ua,w)],va),([x+y for x,y in zip(ub,w)],vb),(w,cancel)]
    queries=[];word=ident
    for u,v in pairs:
        if not any(u) or not any(v) or sum(x*y for x,y in zip(u,v)):raise ValueError('The expanded query is zero or not square zero')
        query=[[ident[i][j]+u[i]*v[j] for j in range(n)] for i in range(n)]
        word=mm(query,word);queries.append(query)
    if word!=[[-x for x in row] for row in ident]:raise ValueError('The expanded query word differs from minus identity')
    return [],pairs,queries


def frame(n,mode):
    if mode=='triangles':
        if n!=4:raise ValueError('The two-triangle control has four fields')
        us,vs,_=e.canonical(2);pairs=[];queries=[]
        for offset in (0,2):
            for u,v in zip(us,vs):
                col=[F(0)]*4;row=[F(0)]*4
                col[offset:offset+2]=list(map(F,u));row[offset:offset+2]=[-2*F(x) for x in v]
                pairs.append((col,row));queries.append([[F(i==j)+col[i]*row[j] for j in range(4)] for i in range(4)])
        return [],pairs,queries
    if mode.startswith('expand-'):
        _,site,kernel=mode.split('-');site,kernel=int(site),int(kernel)
        if not 0<=site<n:raise ValueError('The expansion site is outside the word')
        return expanded_frame(n,site,kernel)
    length=n+2
    t=[F(i+1) if mode=='integer' else F(2**i) for i in range(length-1)]
    t.append(-math.prod(t[::2])/math.prod(t[1::2]))
    g=[[F(0) if i==j else (1 if j>i else -1)*(t[i]+t[j])/t[j] for j in range(length)] for i in range(length)]
    k=[row[:n] for row in g[:n]];kinv=e.inverse(k)
    v=[row[:n] for row in g];u=mm(kinv,g[:n])
    if mm(v,u)!=g:raise ValueError('The Gram rank is not n')
    ident=[[F(i==j) for j in range(n)] for i in range(n)];word=ident;queries=[];pairs=[]
    for i in range(length):
        col=[row[i] for row in u];row=v[i]
        if not any(col) or not any(row) or sum(x*y for x,y in zip(col,row)):raise ValueError('A query is zero or not square zero')
        query=[[ident[a][b]+col[a]*row[b] for b in range(n)] for a in range(n)]
        word=mm(query,word);queries.append(query);pairs.append((col,row))
    if word!=[[-x for x in row] for row in ident]:raise ValueError('The complete query word is not minus identity')
    return t,pairs,queries


def construct(n,mode):
    t,pairs,queries=frame(n,mode);names=[];ps=[]
    def variable(name):
        i=len(names);names.append(name);return {(i,):F(1)}
    for i,(u,v) in enumerate(pairs):
        if i in (0,len(pairs)-1):ps.append(e.zero(n));continue
        pivot=next(j for j,x in enumerate(v) if x);drop=next(j for j,x in enumerate(u) if x and j!=pivot)
        h=[F(int(j==pivot),v[pivot]) for j in range(n)]
        columns=[u,h]+[[int(k==j)-v[j]*h[k] for k in range(n)] for j in range(n) if j not in (pivot,drop)]
        basis=list(map(list,zip(*columns)));a=[variable(f'a{i}_{j}') for j in range(n-2)];b=[variable(f'b{i}_{j}') for j in range(n-2)];c=variable(f'c{i}')
        core=e.zero(n);core[0][0]=core[1][1]=e.const(1);core[0][1]=c
        for j in range(n-2):core[j+2][1]=a[j];core[0][j+2]=b[j]
        ps.append(e.mm(e.mm(e.numeric(basis),core),e.numeric(e.inverse(basis))))
    return t,names,ps,[e.numeric(q) for q in queries]


def coefficient(ps,ts,degree):
    n=len(ps[0]);plus=[e.eye(n)];minus=[e.eye(n)]
    for p,t in zip(ps,ts):
        for word,query in ((plus,e.eye(n)),(minus,t)):
            old=[e.mm(query,a) for a in word]
            word[:]=[old[0]]+[e.ma(old[i],e.mm(p,old[i-1])) for i in range(1,len(old))]
            if len(old)<=degree:word.append(e.mm(p,old[-1]))
    return [p for row in e.ma(plus[degree],minus[degree]) for p in row]


def search(n,mode,degree,seconds):
    start=time.monotonic();t,names,ps,ts=construct(n,mode);steps=[];highest=0;status='UNRESOLVED'
    def expire(signum,frame):raise TimeoutError('The exact scalar elimination time limit expired')
    signal.signal(signal.SIGALRM,expire);signal.alarm(seconds)
    try:
        values=e.linear_solution(coefficient(ps,ts,1))
        if values is None:status='EXACT_SCALAR_FIRST_COEFFICIENT_CONTRADICTION'
        else:
            ps=[e.ms(p,values) for p in ps]
            system=[p for mat in ps for row in e.ma(e.mm(mat,mat),mat,-1) for p in row]
            for d in range(2,degree+1):
                highest=d;system+=coefficient(ps,ts,d)
                while True:
                    reduced=e.row_reduce(system);affine=[p for p in reduced if all(len(w)<=1 for w in p)]
                    step={'degree':d,'variables':len({v for p in system for w in p for v in w}),'rank':len(reduced),'affine':len(affine)};steps.append(step);print(json.dumps(step),flush=True)
                    if any(set(p)=={()} for p in affine):status='EXACT_SCALAR_CONTRADICTION';break
                    if not affine:status='UNRESOLVED' if any(system) else 'COEFFICIENT_FAMILY';break
                    values=e.linear_solution(affine);ps=[e.ms(p,values) for p in ps];system=[e.substitute(p,values) for p in system]
                if 'CONTRADICTION' in status:break
    except TimeoutError:status='UNKNOWN_TIMEOUT'
    finally:signal.alarm(0)
    result={'n':n,'mode':mode,'t':list(map(str,t)),'queries':n+2,'frame_rank':n,'query_word':'EXACT_MINUS_IDENTITY','initial_variables':len(names),'requested_degree':degree,'last_degree':highest,'status':status,'steps':steps,'seconds':time.monotonic()-start,'limit_seconds':seconds,'scope':'scalar commuting coefficients','new_kappa':False}
    if 'CONTRADICTION' not in status and status!='UNKNOWN_TIMEOUT':
        raw=[[[p.get((),F(0)) for p in row] for row in mat] for mat in ps];candidate=[e.numeric(p) for p in raw]
        result['constant_choice_projectors']=not any(p for mat in candidate for row in e.ma(e.mm(mat,mat),mat,-1) for p in row)
        result['constant_choice_full_query']=not any(p for d in range(n+1) for p in coefficient(candidate,ts,d))
        result['constant_choice']=[[[str(x) for x in row] for row in mat] for mat in raw]
        result['chart_variables']=sorted({v for mat in ps for row in mat for p in row for w in p for v in w})
        result['chart']=[[[[[list(w),str(c)] for w,c in sorted(p.items())] for p in row] for row in mat] for mat in ps[1:-1]]
        _,pairs,queries=frame(n,mode)
        result['query_u']=[[str(x) for x in u] for u,v in pairs]
        result['query_v']=[[str(x) for x in v] for u,v in pairs]
        result['raw_queries']=[[[str(x) for x in row] for row in mat] for mat in queries]
        result['inverse_interior_sites']=[]
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--n',type=int,default=4);p.add_argument('--mode',default='integer');p.add_argument('--degree',type=int);p.add_argument('--seconds',type=int,default=180);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.n<2 or a.n%2:raise ValueError('The input count must be positive and even')
    if a.mode not in ('integer','dyadic','triangles') and not a.mode.startswith('expand-'):raise ValueError('Unknown query frame mode')
    r=search(a.n,a.mode,a.degree or a.n,a.seconds);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ('constant_choice','chart','raw_queries','query_u','query_v')}))
