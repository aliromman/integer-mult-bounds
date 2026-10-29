"""Check the complete first stage coefficient of three exact query blocks.

This is a fixed-frame screen. Passing it does not give full projectors.
"""
import argparse
import json
from pathlib import Path
import signal
import sys
import time
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
from exact_replay import canonical


def build(block,order):
    n=3*block;frames=[]
    us,vs,_=canonical(block)
    sites=[(c,i) for i in range(block+1) for c in range(3)] if order=='round' else [(c,i) for c in range(3) for i in range(block+1)]
    variable_count=len(sites)*(2*n-3)
    ps=[];rs=[];labels=[];offset=0
    for c,i in sites:
        u=np.zeros(n,dtype=np.int64);v=u.copy();u[c*block:(c+1)*block]=us[i];v[c*block:(c+1)*block]=vs[i]
        pivot=next(j for j,x in enumerate(v) if x);drop=next(j for j,x in enumerate(u) if x and j!=pivot)
        h=np.eye(n,dtype=np.int64)[:,pivot]*v[pivot]
        cols=[u,h]+[np.eye(n,dtype=np.int64)[:,j]-v[j]*h for j in range(n) if j not in (pivot,drop)]
        basis=np.array(cols).T;inv=np.array(sp.Matrix(basis).inv(),dtype=np.int64)
        if not np.array_equal(basis@inv,np.eye(n,dtype=np.int64)):
            raise ValueError('The integral chart basis is not invertible')
        p=np.zeros((variable_count+1,n,n),dtype=np.int64)
        p[0]=np.outer(basis[:,0],inv[0])+np.outer(basis[:,1],inv[1])
        for j in range(n-2):
            p[1+offset+j]=np.outer(basis[:,j+2],inv[1])
            p[1+offset+n-2+j]=np.outer(basis[:,0],inv[j+2])
        p[1+offset+2*(n-2)]=np.outer(basis[:,0],inv[1]);offset+=2*n-3
        ps.append(p);rs.append(np.outer(u,v));labels.append(c)
    return ps,rs,labels


def run(block,order,seconds,save_chart=None):
    started=time.monotonic()
    def expire(*_):
        raise TimeoutError('The exact first-coefficient time limit expired')
    signal.signal(signal.SIGALRM,expire)
    signal.alarm(seconds)
    n=3*block
    try:
        ps,rs,labels=build(block,order);vars=len(ps[0])-1;ident=np.eye(n,dtype=np.int64)
        chars=[(1,1,1),(-1,1,-1),(1,-1,-1),(-1,-1,1)]
        coeff=[];out=[]
        for ch in chars:
            ts=[ident+(ch[c]-1)*r for c,r in zip(labels,rs)];prefix=[ident]
            for t in ts:prefix.append(t@prefix[-1])
            expected=np.diag(np.repeat(ch,block))
            if not np.array_equal(prefix[-1],expected):raise ValueError('The query control fails')
            first=np.zeros_like(ps[0]);suffix=ident
            for site in range(len(ps)-1,-1,-1):
                first+=suffix@ps[site]@prefix[site+1];suffix=suffix@ts[site]
            coeff.append(first);out.append(expected)
        rows=[]
        for ch in range(1,4):
            diff=coeff[ch]-out[ch]@coeff[0]
            rows.extend(diff[:,i,j].tolist() for i in range(n) for j in range(n) if np.any(diff[:,i,j]))
        rows.extend(coeff[0][:,i,j].tolist() for i in range(n) for j in range(n) if i//block!=j//block and np.any(coeff[0][:,i,j]))
        matrix=sp.Matrix(rows); augmented=matrix[:,1:].row_join(-matrix[:,0]);rref,pivots=augmented.rref()
        if vars in pivots: status='EXACT_FIRST_COEFFICIENT_CONTRADICTION';cross=None
        else:
            free=[i for i in range(vars) if i not in pivots];sol=sp.zeros(vars,len(free)+1)
            for j,f in enumerate(free):sol[f,j+1]=1
            for row,p in enumerate(pivots):
                sol[p,0]=rref[row,vars]
                for j,f in enumerate(free):sol[p,j+1]=-rref[row,f]
            if augmented[:,:vars]*sol[:,0]!=augmented[:,vars]:raise ValueError('The affine solution fails')
            if free and augmented[:,:vars]*sol[:,1:]!=sp.zeros(len(rows),len(free)):raise ValueError('The homogeneous solution fails')
            cross=0
            for p,c in zip(ps,labels):
                for i in range(n):
                    for j in range(n):
                        if i//block==j//block:continue
                        value=sp.Matrix([p[k,i,j] for k in range(1,vars+1)]).T*sol
                        value[0]+=int(p[0,i,j]);cross+=int(any(value))
            status='FIRST_COEFFICIENT_FAMILY'
            if save_chart:
                solution=np.vstack([np.r_[1.,np.zeros(len(free))],np.array(sol,float)])
                reduced=np.array([np.einsum('vij,vk->kij',p,solution) for p in ps])
                np.savez_compressed(save_chart,projectors=reduced,queries=np.array(rs),labels=np.array(labels),block=block)
        return {'block':block,'fields':n,'full_calls':len(ps),'order':order,'variables':vars,'rank':len(pivots),
                'free_variables':vars-len(pivots),'remaining_off_block_entries':cross,'status':status,
                'seconds':time.monotonic()-started,'limit_seconds':seconds,'complete_query_control':'PASS','new_kappa':False}
    except TimeoutError:return {'block':block,'order':order,'status':'UNKNOWN_TIMEOUT','limit_seconds':seconds,'new_kappa':False}
    finally:signal.alarm(0)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--block',type=int,choices=[2,4,6],default=2)
    p.add_argument('--order',choices=['round','block'],default='round');p.add_argument('--seconds',type=int,default=120)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--save-chart',type=Path);a=p.parse_args()
    result=run(a.block,a.order,a.seconds,a.save_chart);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
