"""Price the proposed eighteen-field boundary with exact rational intervals.

This omits unknown additional fees. It is not an attained caller invoice.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import runpy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def check():
    manifest=json.loads((HERE/'SOURCES.json').read_text())
    for item in manifest['inputs']:
        raw=(ROOT/item['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=item['sha256']:
            raise ValueError('A scientific source differs')
    bounds=runpy.run_path(str(ROOT/'research/alternating-caller-closure/price_check.py'),run_name='price_bounds')
    source=json.loads((ROOT/'research/alternating-caller-closure/sources/selected6030-result.json').read_text())['local']
    c={int(k):-21*v for k,v in source['paid'].items()}
    for k,v in source['helpers'].items():c[int(k)]=c.get(int(k),0)+21*v
    for k,v in [(17,104*12),(25,104*6)]:c[k]=c.get(k,0)+v
    records=[]
    for k in ['0.003125804','0.0034','0.00343211','0.00343212']:
        p={r:bounds['power_bounds'](r,Fraction(k)) for r in c}
        interval=(sum(v*p[r][0 if v>0 else 1] for r,v in c.items()),sum(v*p[r][1 if v>0 else 0] for r,v in c.items()))
        records.append({'kappa':k,'raw_margin':bounds['outward'](interval),'positive':interval[0]>0,'negative':interval[1]<0})
    return {'source':source,'full_calls':21,'main_targets':{'17':1248,'25':624},'records':records,
            'rank_scale':'A common factor (2m)^(1-kappa) is removed from both sides',
            'other_fees':'unknown and omitted','attained':False,'new_kappa':False}


if __name__=='__main__':print(json.dumps(check(),indent=2))
