"""Check favorable caller price upper bounds with rational intervals.

The source profile is retained unchanged. No physical caller is constructed.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def power_bounds(rank,kappa):
    if rank==1:return F(1),F(1)
    x=F(rank-1,rank+1);term=x;lower=F(0);terms=128
    for j in range(terms):lower+=2*term/(2*j+1);term*=x*x
    upper=lower+2*term/((2*terms+1)*(1-x*x))
    a,b=kappa*lower,kappa*upper
    if not 0<=a<=b<=1:raise ValueError('The exponential interval is outside its controlled range')
    def series(value,length):
        out=term=F(1)
        for j in range(1,length+1):term*=-value/j;out+=term
        return out
    return rank*series(b,9),rank*series(a,8)


def outward(interval):
    scale=2**160;lo,hi=interval
    return [str(F((lo*scale).__floor__(),scale)),str(F((hi*scale).__ceil__(),scale))]


def verify():
    pin=json.loads((ROOT/'SELECTED-LOCAL-PROFILE.json').read_text())
    data=(ROOT/'sources/selected6030-result.json').read_bytes()
    if hashlib.sha256(data).hexdigest()!=pin['source_sha256']:raise ValueError('The selected source hash differs')
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if blob!=pin['git_blob']:raise ValueError('The selected source Git blob differs')
    source=json.loads(data)
    for key in ('paid','helpers'):
        if source['local'][key]!=pin[key]:raise ValueError('The local profile differs from the retained source')
    kappa=F(pin['kappa']);prices={r:power_bounds(r,kappa) for r in range(1,26)}
    paid={int(k):v for k,v in pin['paid'].items()};helpers={int(k):v for k,v in pin['helpers'].items()}
    c=(sum(v*prices[r][0] for r,v in paid.items())-sum(v*prices[r][1] for r,v in helpers.items()),
       sum(v*prices[r][1] for r,v in paid.items())-sum(v*prices[r][0] for r,v in helpers.items()))
    if sum(r*v for r,v in paid.items())-sum(r*v for r,v in helpers.items())!=1720:raise ValueError('The source linear mass differs')
    results=[]
    for n in (2,4,6):
      for boundaries in (0,1,2):
        full=n//2+boundaries;thin=n+1-full
        low,high=divmod(2*full,n)
        targets={8*low+1:n-high}
        if high:targets[8*(low+1)+1]=high
        margin=(104*sum(v*prices[r][0] for r,v in targets.items())-full*c[1]-104*thin,
                104*sum(v*prices[r][1] for r,v in targets.items())-full*c[0]-104*thin)
        if margin[1]>=0:raise ValueError('The short-word price upper bound is not negative')
        results.append({'n':n,'full_boundaries':boundaries,'scope':'all boundary orientations; favorable balanced absolute-power bound','target_ranks_per_port':targets,'margin_interval':outward(margin)})
    for n in (8,10):
        full=n//2+2;thin=n//2-1
        target_upper=104*(n+8*(n+2))
        upper=target_upper-full*c[0]-104*thin
        if upper>=0:raise ValueError('The mixed-boundary rank-mass price bound is not negative')
        results.append({'n':n,'scope':'one positive and one inverse boundary','target_rank_mass_upper':target_upper,'margin_upper':outward((upper,upper))[1]})
    # Independent small controls for the alternating series convention.
    for rank in (2,9,17):
        if power_bounds(rank,F(0))!=(F(rank),F(rank)):raise ValueError('The zero-exponent-discount control fails')
    return {'kappa':str(kappa),'source_commit':pin['source_commit'],'full_net_price_interval':outward(c),'cases':results,'extra_fees':'omitted in the favorable bound','new_kappa':False}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
    print('PASS_SELECTED_ALTERNATING_FAMILY_PRICE_BOUNDS')
