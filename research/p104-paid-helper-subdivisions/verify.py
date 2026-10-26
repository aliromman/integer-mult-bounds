"""Verify seventeen paid helper subdivisions and the complete conditional price.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
The source archive retains all predecessor notices.
"""
import argparse
from collections import Counter
from contextlib import ExitStack
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
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


def reject(action):
    try:
        action()
    except ValueError as error:
        return str(error)
    raise ValueError('Invalid construction accepted')


def declared_edits(source, name, edits):
    text = source.read_text(encoding='utf-8')
    for before, after in json.loads((HERE/edits).read_bytes()):
        need(text.count(before) == 1, 'Unique declared edit: '+name)
        text = text.replace(before, after)
    need(text == (HERE/name).read_text(encoding='utf-8'), 'Only declared edits: '+name)


def controls(source, plan, word, L):
    request = deepcopy(plan[0]['request'])
    request['cut_correction']['146'] = '-1/2'
    reasons = {'wrong_retirement_sign': reject(lambda: construct.construct(source, request))}
    request = deepcopy(plan[0]['request'])
    request['cut'] = 0
    reasons['crossed_scalar_dependency'] = reject(lambda: construct.construct(source, request))
    events = {
        'missing_paid_subdivision': ['MOVE', 744, 531, 127, False],
        'missing_output_restoration': ['ADD', 115, 744, '4', 127],
    }
    for label, event in events.items():
        need(word['program'].count(event) == 1, 'Unique actual control event')
        bad = deepcopy(word)
        bad['program'].remove(event)
        reasons[label] = reject(lambda: L.replay(bad))
    bad = deepcopy(word)
    bad['initial'][744] = 529
    reasons['wrong_fresh_input_chart'] = reject(lambda: L.replay(bad))
    bad = deepcopy(word)
    bad['final'][731] = 127
    reasons['uncharged_original_H9_target'] = reject(lambda: L.replay(bad))
    bad = deepcopy(word)
    event = next(e for e in reversed(bad['program']) if e[0] == 'ADD')
    event[3] = str(-Q(event[3]))
    reasons['wrong_terminal_correction_sign'] = reject(lambda: L.replay(bad))
    bad = deepcopy(word)
    event = next(e for e in reversed(bad['program']) if e[0] == 'CLEAR')
    bad['program'].remove(event)
    reasons['missing_owned_clear'] = reject(lambda: L.replay(bad))
    return reasons


def run():
    with ExitStack() as contexts:
        parent, _ = contexts.enter_context(sources())
        predecessor_sources = load(parent/'sources.py', 'paid_split_predecessor_sources')
        predecessor, _ = contexts.enter_context(predecessor_sources.sources())
        cursor = predecessor
        for depth in range(3):
            module = load(cursor/'sources.py', 'paid_split_sources_'+str(depth))
            cursor, root = contexts.enter_context(module.sources())
        pins = json.loads((HERE/'SOURCES.json').read_bytes())
        need(sha256((cursor/'local.py').read_bytes()).hexdigest() == pins['derived_local_sha256'],
             'Exact authenticated owned-aware local engine')
        declared_edits(cursor/'local.py', 'local.py', 'LOCAL-EDITS.json')
        declared_edits(parent/'pricing.py', 'pricing.py', 'PRICING-EDITS.json')
        L = load(HERE/'local.py', 'local')
        L.ROOT = root
        prior = load(parent/'construct.py', 'paid_split_predecessor_constructor')
        source, prior_witness = prior.construct(L.read(predecessor/'WORD.json.gz'))
        need(source == L.read(parent/'WORD.json.gz')
             and construct.encode(prior_witness) == L.read(parent/'WITNESS.json'),
             'Reconstructed exact PR5718 source and witness')
        plan = L.read(HERE/'PLAN.json')
        word, witnesses = construct.build(source, plan)
        local = L.replay(word)
        columns = L.reverse_columns(word)
        need(local['public_inputs'] == 761 and local['owned_fields'] == 9
             and local['owned_lifetimes'] == 11 and local['peak_owned_fields'] == 2
             and len(word['program']) == 27502, 'Exact final finite inventory')
        need(local['paid'] == {1:2202, 2:458, 3:495, 4:92, 5:175, 6:6, 7:15}
             and local['helpers'] == {1:19, 2:5, 4:91, 5:6, 7:6, 8:32, 9:394}
             and local['real_scalar_allowance'] == 244338, 'Exact local invoices')
        previous = L.replay(source)
        paid_delta, target_delta = Counter(local['paid']), Counter(local['helpers'])
        paid_delta.subtract(previous['paid'])
        target_delta.subtract(previous['helpers'])
        need({r:n for r,n in paid_delta.items() if n} == {1:18, 2:-1, 5:16, 6:-16}
             and {r:n for r,n in target_delta.items() if n} == {1:17, 8:17, 9:-17},
             'Seventeen signed paid and productive deltas')
        faults = controls(source, plan, word, L)
        print('PASS 17 constructions, 770 columns, actual charts, lifetimes and fault controls', flush=True)
        native = L.native(word)
        need(native['distinct_maps'] == 2912, 'Complete native map inventory')
        C = load(root/'research/p104-productive-controller/caller.py', 'paid_split_caller')
        caller = C.construct()
        need(caller == L.read(root/'research/p104-productive-controller/CONTROLLER.json'),
             'Retained complete physical caller')
        physical = C.physical(caller, word)
        print('PASS 2912 native maps with phases and complete physical caller', flush=True)
        P = load(HERE/'pricing.py', 'paid_split_price')
        price = P.run(local, native, caller, physical)
        faults.update(C.controls(caller))
        bad = deepcopy(local)
        bad['paid'][1] -= 1
        faults['unpaid_rank1_query'] = reject(lambda: P.assemble(bad, native, caller, physical))
        bad = deepcopy(local)
        bad['helpers'][9] += 1
        faults['invented_extra_H9_target'] = reject(lambda: P.assemble(bad, native, caller, physical))
        result = {
            'status': 'PASS_P104_PAID_HELPER_SUBDIVISIONS_KAPPA_00311224341',
            'source_commit': pins['source_export']['commit'],
            'canonical_word_sha256': construct.canonical_hash(word),
            'steps': len(witnesses), 'instructions': len(word['program']),
            'local': local, 'separate_columns': columns, 'native': native,
            'physical_caller': physical, 'complete_cost': price,
            'delta_paid_local': {r:n for r,n in paid_delta.items() if n},
            'delta_targets_local': {r:n for r,n in target_delta.items() if n},
            'negative_controls': faults,
            'limits': [
                'The five inherited framework assumptions remain.',
                'The improvement is modest and conditional; no selected-frontier admission is claimed.',
                'The huge uniform word and complete rounded multiplier are not materialized.',
                'The final package has not received independent human or formal acceptance.',
                'The fixed discovery search is not a global optimality proof.']}
        return P.encode(result), witnesses, word


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='p104-paid-split-bytecode-') as cache:
        sys.pycache_prefix = cache
        result, witness, word = run()
    for name, value in (('RESULT.json', result), ('WITNESS.json', witness), ('WORD.json.gz', word)):
        path = HERE/name
        if args.write:
            raw = (json.dumps(value, sort_keys=True, indent=2)+'\n').encode()
            path.write_bytes(gzip.compress(raw, mtime=0) if name.endswith('.gz') else raw)
        else:
            raw = path.read_bytes()
            need(json.loads(gzip.decompress(raw) if name.endswith('.gz') else raw) == value,
                 'Complete retained output: '+name)
    print(result['status'])


if __name__ == '__main__':
    main()
