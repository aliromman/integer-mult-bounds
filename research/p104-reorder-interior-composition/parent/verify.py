"""Check the basis certificates and the complete conditional price.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
import argparse
from collections import Counter
from copy import deepcopy
import gzip
import json
from pathlib import Path
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sources import sources, load, need
import construct


def rejected(action):
    try:
        action()
    except ValueError as error:
        return str(error)
    raise ValueError('An invalid certificate was accepted')


def difference(after, before):
    result = Counter(after)
    result.subtract(before)
    return {r: count for r, count in result.items() if count}


def run():
    witnesses = json.loads((HERE / 'WITNESS.json').read_bytes())
    need(len(witnesses) == 4, 'Four fixed basis changes')
    with sources() as (parent, outer_root):
        retained = load(parent / 'sources.py', 'interior_retained_sources')
        with retained.sources() as (base, root):
            L = load(base / 'local.py', 'local')
            source = L.read(parent / 'WORD.json.gz')
            old = L.replay(source)
            word, steps = construct.construct(source, witnesses, L)
            local = L.replay(word)
            expected = [{1: -1, 2: -1, 3: 1}] + [{1: -2, 2: 1}] * 3
            previous = old
            for receipt, delta in zip(steps, expected):
                need(difference(receipt['paid'], previous['paid']) == delta,
                     'Exact paid delta at each basis change')
                need(receipt['helpers'] == old['helpers'], 'All productive targets retained')
                previous = receipt
            need(difference(local['paid'], old['paid']) == {1: -7, 2: 2, 3: 1},
                 'Complete local paid delta')
            need((local['public_inputs'], local['owned_fields'], local['owned_lifetimes'],
                  local['peak_owned_fields']) == (743, 9, 11, 2), 'Public and owned inventories')
            columns = L.reverse_columns(word)
            print('PASS four basis changes, eight scalar blocks and all 565504 endpoint entries', flush=True)
            native = L.native(word)
            C = load(root / 'research/p104-productive-controller/caller.py', 'interior_caller')
            caller = C.construct()
            need(caller == L.read(root / 'research/p104-productive-controller/CONTROLLER.json'),
                 'Same complete caller')
            physical = C.physical(caller, word)
            print('PASS all signed native maps and the complete physical caller', flush=True)
            price_source = (parent / 'pricing.py').read_text()
            for before, after in L.read(HERE / 'PRICING-EDITS.json'):
                need(price_source.count(before) == 1, 'One declared price edit')
                price_source = price_source.replace(before, after)
            need(price_source == (HERE / 'pricing.py').read_text(), 'Only declared price edits')
            P = load(HERE / 'pricing.py', 'interior_price')
            price = P.run(local, native, caller, physical)
            controls = C.controls(caller)
            bad = deepcopy(witnesses[0])
            bad['coefficient'] = str(-L.Q(bad['coefficient']))
            controls['wrong_basis_sign'] = rejected(lambda: construct.step(source, bad, L))
            bad = deepcopy(witnesses[0])
            bad['blocks'][0]['after'][0][-2] = '17'
            controls['wrong_scalar_block'] = rejected(lambda: construct.step(source, bad, L))
            bad = deepcopy(witnesses[0])
            bad['blocks'].pop()
            controls['missing_conjugated_block'] = rejected(lambda: construct.step(source, bad, L))
            bad = deepcopy(word)
            bad['program'].remove(next(e for e in bad['program'] if e[0] == 'MOVE'))
            controls['missing_native_transport'] = rejected(lambda: L.replay(bad))
            bad = deepcopy(word)
            bad['program'].pop(max(i for i, e in enumerate(bad['program']) if e[0] == 'CLEAR'))
            controls['missing_owned_erasure'] = rejected(lambda: L.replay(bad))
            bad = deepcopy(local)
            bad['paid'][1] -= 1
            controls['unpaid_query'] = rejected(lambda: P.assemble(bad, native, caller, physical))
            bad = deepcopy(local)
            bad['helpers'][1] += 1
            controls['extra_target'] = rejected(lambda: P.assemble(bad, native, caller, physical))
            result = dict(
                status='PASS_P104_PAID_INTERIOR_BASES_KAPPA_00310815849',
                source_commit='3b0b4f37d5005f381ed9cfcccebfbe1ab35c71aa',
                authenticated_outer_files=23, authenticated_nested_files=36,
                source_local=old, steps=steps, local=local, separate_columns=columns,
                local_paid_delta={1: -7, 2: 2, 3: 1}, local_helper_delta={},
                native=native, physical_caller=physical, complete_cost=price,
                negative_controls=controls,
                limits=[
                    'The five inherited assumptions and helper endpoint contract remain.',
                    'The huge uniform word and complete rounded multiplier are not materialized.',
                    'Author checks do not establish independent scientific acceptance.',
                    'This gain is small. The requested large gain and a global optimum remain unproved.'])
            return P.encode(result), word


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='p104-interior-bytecode-') as cache:
        sys.pycache_prefix = cache
        result, word = run()
    for name, value in (('RESULT.json', result), ('WORD.json.gz', word)):
        path = HERE / name
        payload = (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()
        if args.write:
            path.write_bytes(gzip.compress(payload, mtime=0) if name.endswith('.gz') else payload)
        else:
            raw = path.read_bytes()
            need(json.loads(gzip.decompress(raw) if name.endswith('.gz') else raw) == value,
                 'Exact retained output: ' + name)
    print(result['status'])


if __name__ == '__main__':
    main()
