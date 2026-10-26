"""Check the exact composition, all physical maps and complete kappa.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
import argparse
from collections import Counter
from copy import deepcopy
from hashlib import sha256
import gzip,json,sys,tempfile,shutil
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from sources import sources,load,need,authenticate
import construct

def rejected(action):
    try:action()
    except ValueError as error:return str(error)
    raise ValueError('An invalid construction was accepted')

def difference(a,b):
    out=Counter(a);out.subtract(b);return {r:n for r,n in out.items()if n}

def run():
    pins=authenticate()
    with sources() as(parent,contact,base,root):
        L=load(base/'local.py','local')
        old_constructor=load(parent/'construct.py','composition_interior_constructor')
        source=L.read(parent/'WORD.json.gz')
        rebuilt,_=old_constructor.construct(L.read(contact/'WORD.json.gz'),L.read(parent/'WITNESS.json'),L)
        need(rebuilt==source,'Every paid-interior source change is rebuilt')
        reorder_source=L.read(HERE/'reorder/WORD.json.gz')
        need(construct.construct(L.read(contact/'WORD.json.gz'),reorder_source['schema'])==reorder_source,'Released R22 source word is reproduced')
        word=construct.construct(source);old=L.replay(source);local=L.replay(word)
        trace=construct.trace_equivalent(source,word)
        need(difference(local['paid'],old['paid'])=={1:-2,3:1},'Exact changed paid requests')
        need(difference(local['helpers'],old['helpers'])=={8:-1,9:1},'Exact changed independent targets')
        need(local['public_inputs']==743 and local['owned_fields']==9 and local['owned_lifetimes']==11 and local['peak_owned_fields']==2,'Full input and owned inventory')
        need(local['real_scalar_allowance']==old['real_scalar_allowance'],'No uncharged scalar operation')
        columns=L.reverse_columns(word)
        print('PASS exact composition, trace and all 565504 endpoint entries',flush=True)
        census=load(HERE/'cuts.py','composition_cut_census');cuts=census.run(word)
        positive=L.read(base/'inputs/PR5580-WORD.json.gz')
        main=set(positive['live'])|{r+112 for r in positive['live']}
        group={r for r in range(positive['public'])if r not in main and (positive['initial'][r],positive['final'][r])==(3,127)}
        why=rejected(lambda:census.scan(positive,group,3,127))
        need(why=='No independent helper pairing remains','The retained positive split fixture must be detected')
        print('PASS full cut census',cuts['contexts'],flush=True)
        native=L.native(word)
        print('PASS native phases',native['distinct_maps'],flush=True)
        C=load(root/'research/p104-productive-controller/caller.py','composition_caller')
        caller=C.construct();need(caller==L.read(root/'research/p104-productive-controller/CONTROLLER.json'),'Retained physical controller')
        physical=C.physical(caller,word)
        price_source=(parent/'pricing.py').read_text()
        for before,after in L.read(HERE/'PRICING-EDITS.json'):
            need(price_source.count(before)==1,'One declared price edit');price_source=price_source.replace(before,after)
        need(price_source==(HERE/'pricing.py').read_text(),'Only declared price edits')
        P=load(HERE/'pricing.py','composition_price');price=P.run(local,native,caller,physical)
        controls=C.controls(caller)
        bad=deepcopy(source);bad['initial'][473]=3
        controls['wrong_source_start']=rejected(lambda:construct.construct(bad))
        bad=deepcopy(word);bad['initial'][473]=594
        controls['lost_endpoint_change']=rejected(lambda:L.replay(bad))
        bad=deepcopy(word);bad['program'].pop(next(i for i,e in enumerate(bad['program'])if e[:3]==['ADD',213,473]))
        controls['lost_reordered_contact']=rejected(lambda:L.replay(bad))
        bad=deepcopy(local);bad['paid'][1]-=1
        controls['unpaid_request']=rejected(lambda:P.assemble(bad,native,caller,physical))
        bad=deepcopy(local);bad['helpers'][9]+=1
        controls['extra_target']=rejected(lambda:P.assemble(bad,native,caller,physical))
        result=dict(status='PASS_P104_REORDER_INTERIOR_KAPPA_00310819637',authenticated_inputs=pins,source_local=old,local=local,trace=trace,separate_columns=columns,cut_census=cuts,native=native,physical_caller=physical,complete_cost=price,negative_controls=controls,limits=['The five inherited framework assumptions and helper endpoint contract remain.','The huge uniform word and rounded multiplier are not materialized.','Author checks are not independent scientific acceptance.','The cut result is scoped to this emitted word and common-frame helper boundary groups.'])
        return P.encode(result),word

def source_controls():
    entries=json.loads((HERE/'INPUTS.json').read_bytes())
    with tempfile.TemporaryDirectory(prefix='p104-composition-corruption-')as directory:
        root=Path(directory)/'package';shutil.copytree(HERE,root,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for item in entries:
            path=root/item['local'];raw=path.read_bytes();need(raw,'Nonempty source fixture')
            path.write_bytes(bytes([raw[0]^1])+raw[1:])
            need(rejected(lambda:authenticate(root))=='Frozen input '+item['local'],'Actual retained source corruption rejects')
            path.write_bytes(raw)
    print('PASS source corruption controls',len(entries),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true');args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='p104-composition-bytecode-')as cache:
        sys.pycache_prefix=cache;result,word=run();source_controls()
    for name,value in [('RESULT.json',result),('WORD.json.gz',word)]:
        path=HERE/name;raw=(json.dumps(value,sort_keys=True,indent=2)+'\n').encode()
        if args.write:path.write_bytes(gzip.compress(raw,mtime=0)if name.endswith('.gz')else raw)
        else:
            saved=path.read_bytes();need(json.loads(gzip.decompress(saved)if name.endswith('.gz')else saved)==value,'Exact rebuilt '+name)
    authenticate()
    print(result['status'],flush=True)
if __name__=='__main__':main()
