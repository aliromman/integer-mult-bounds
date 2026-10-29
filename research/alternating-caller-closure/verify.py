"""Replay the all-length proof inputs and exact finite certificates."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
PINNED={'inputs/exact_replay.py': '83352900ee81cc79b1ba438cb87762ccc73ee323dd027be1389ef2eec245535d', 'inputs/free_replay.py': '85ca11c5a752680c6a3a97e1aec55f06dab9a497371126f70670c333d99dadba', 'sources/selected6030-result.json': '84fcbe2a0b64307b5703c363408042bdd91555b0419509aa50daa43b3b1924ad'}


def verify():
    for path,expected in PINNED.items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError('A frozen proof input differs: '+path)
    # A fresh cache path prevents old bytecode from replacing checked source.
    with tempfile.TemporaryDirectory(prefix='imb-closure-bytecode-') as cache:
        previous=sys.pycache_prefix;sys.pycache_prefix=cache
        try:
            sys.path.insert(0,str(ROOT))
            import trace_verify
            import rank_curvature_check
            import price_check
            result={'trace':trace_verify.verify(),'curvature':rank_curvature_check.verify(),
                    'prices':price_check.verify(),'new_kappa':False}
        finally:sys.pycache_prefix=previous
    print(json.dumps(result,indent=2))
    print('PASS_ALL_LENGTH_ALTERNATING_CALLER_CLOSURE')


if __name__=='__main__':verify()
