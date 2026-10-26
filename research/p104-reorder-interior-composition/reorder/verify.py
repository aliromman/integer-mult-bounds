"""Check the R22 reorder+reframe word on PR5660's P104 word on PR5660's own complete stack.

Joel Pulikkan (GamingPuzzled) with Claude Opus 5.5 (Anthropic) via Claude Code. Apache-2.0.
Every retained checker is loaded unmodified from ../p104-contact-rewrite (Ali Romman, zxv and
the contributors listed in its NOTICE); only the four declared pricing literals change.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse, gzip, importlib.util, json, sys, tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE = HERE.parent/'p104-contact-rewrite'
sys.path.insert(0, str(HERE))
from construct import construct

PINS = {   # PR5660 head 3b0b4f37d5005f381ed9cfcccebfbe1ab35c71aa
    'WORD.json.gz': '9e78864cd46fc443c6183674588c99c52cfb8e8e408c7ba09bb50cf4694463f6',
}
PUBLIC = 743
PARENT_KAPPA = Q('.00310676333')
LITERALS = [
    ("KAPPA, TESTED = Q('.00310676333'), Q('.003106763335')", "KAPPA, TESTED = Q('.00310680117'), Q('.003106801175')"),
    ("PREVIOUS = Q('.00310644469')", "PREVIOUS = Q('.00310676333')"),
    ("sum(paid.values()) == 66474", "sum(paid.values()) == 66456"),
    ("'kappa_decimal':'.00310676333'", "'kappa_decimal':'.00310680117'"),
]
STATUS = 'PASS_P104_R22_REFRAME_KAPPA_00310680117'


def need(ok, why):
    if not ok:
        raise ValueError(why)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module
    spec.loader.exec_module(module); return module


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def kinds(e):
    """(register, type) per touched register; T/T and S/S commute exactly, everything else conflicts.
    Owned registers keep every operation exclusive; INIT/CLEAR also share the owned-lifetime token."""
    op = e[0]
    if op == 'ADD':
        out = [(e[1], 'T'), (e[2], 'S')]
    elif op == 'COPY':
        out = [(e[1], 'X'), (e[2], 'S')]
    elif op in ('INIT', 'CLEAR'):
        return [(e[1], 'X'), ('owned-lifetime', 'X')]
    else:
        out = [(e[1], 'X')]
    return [(r, 'X' if r >= PUBLIC else t) for r, t in out]


def strip(e):
    return tuple(e[:2]) if e[0] in ('INIT', 'CLEAR') else tuple(e[:-1])


def trace_equivalent(old, new):
    """Every scalar/owned operation appears once more; every dependent pair keeps its order."""
    need(Counter(map(strip, old)) == Counter(map(strip, new)), 'Same multiset of scalar and owned operations')
    slots = defaultdict(list)
    for j, e in enumerate(new):
        slots[strip(e)].append(j)
    for v in slots.values():
        v.reverse()
    match = [slots[strip(e)].pop() for e in old]
    groups = defaultdict(list)
    for i, e in enumerate(old):
        for r, t in kinds(e):
            g = groups[r]
            if g and g[-1][0] == t and t != 'X':
                g[-1][1].append(i)
            else:
                g.append([t, [i]])
    for g in groups.values():
        for (_, a), (_, b) in zip(g, g[1:]):
            need(max(match[i] for i in a) < min(match[i] for i in b), 'Dependent operations keep their order')
    return match


def changes(a, b):
    d = Counter(a); d.subtract(b); return {str(r): c for r, c in sorted(d.items()) if c}


def reject(f):
    try:
        f()
    except (ValueError, AssertionError, KeyError, IndexError) as e:
        return str(e) or type(e).__name__
    raise ValueError('Malformed control accepted')


def run():
    for name, digest in PINS.items():
        need(sha256((BASE/name).read_bytes()).hexdigest() == digest, 'Exact PR5650 source ' + name)
    base, word = read(BASE/'WORD.json.gz'), read(HERE/'WORD.json.gz')
    need(construct(base) == word, 'Exact regeneration from PR5660')
    for key in base:
        if key not in ('program', 'initial', 'schema'):
            need(word[key] == base[key], 'Unchanged word field ' + key)
    starts = {str(r): [a, b] for r, (a, b) in enumerate(zip(base['initial'], word['initial'])) if a != b}
    need(starts == {'473': [594, 3]}, 'Only helper 473 start changes')
    need(473 not in set(base['live']), 'Changed start is a free helper endpoint')
    scalar = lambda w: [e for e in w['program'] if e[0] != 'MOVE']
    match = trace_equivalent(scalar(base), scalar(word))
    old = scalar(base)
    charts = [[i, old[i][:-1], old[i][-1], scalar(word)[match[i]][-1]] for i in range(len(old))
              if scalar(word)[match[i]] != old[i]]
    print('PASS source pin, regeneration, trace equivalence; changed charts', charts, flush=True)

    sys.path.insert(0, str(BASE))
    S = load(BASE/'sources.py', 'r22_sources')
    original = (BASE/'pricing.py').read_text()
    for before, after in LITERALS:
        need(original.count(before) == 1, 'Unique pricing literal ' + before); original = original.replace(before, after)
    need((HERE/'pricing.py').read_text() == original, 'Only declared pricing literals change')
    with S.sources() as (parent, framework):
        L = load(parent/'local.py', 'local')
        L.ROOT = framework
        parent_local, local = L.replay(base), L.replay(word)
        need(parent_local['real_scalar_allowance'] == local['real_scalar_allowance'], 'Same scalar allowance')
        paid_delta = changes(local['paid'], parent_local['paid'])
        helper_delta = changes(local['helpers'], parent_local['helpers'])
        need(paid_delta == {'1': -2, '3': 1} and helper_delta == {'8': -1, '9': 1}, 'Exact local invoice change')
        print('PASS exact local replay', local['paid'], local['helpers'], flush=True)
        transpose = L.reverse_columns(word)
        print('PASS every dense transpose column', transpose['columns'], flush=True)
        native = L.native(word)
        print('PASS every native map and exact phase', native['distinct_maps'], flush=True)
        C = load(framework/'research/p104-productive-controller/caller.py', 'r22_caller')
        caller = C.construct()
        need(caller == read(framework/'research/p104-productive-controller/CONTROLLER.json'), 'Exact retained caller')
        physical = C.physical(caller, word)
        print('PASS full physical caller', flush=True)
        P = load(HERE/'pricing.py', 'r22_pricing')
        price = P.run(local, native, caller, physical)
        print('PASS complete pricing', price['raw_root_interval'], flush=True)
        lo = Q(price['raw_root_interval'][0])
        need(lo > PARENT_KAPPA, 'Strict raw gain over PR5660 claim')
        controls = {}
        bad = deepcopy(word); bad['initial'][473] = 594
        controls['old_helper_start'] = reject(lambda: L.replay(bad))
        bad = deepcopy(word); i = bad['program'].index(['ADD', 213, 473, '1/2', 3]); bad['program'][i][-1] = 594
        controls['read_on_old_chart'] = reject(lambda: L.replay(bad))
        bad = deepcopy(word); bad['program'].remove(['MOVE', 473, 3, 716, False])
        controls['missing_helper_move'] = reject(lambda: L.replay(bad))
        bad = deepcopy(local); bad['paid'][1] -= 1
        controls['unpaid_query'] = reject(lambda: P.assemble(bad, native, caller, physical))
        bad = deepcopy(local); bad['helpers'][9] += 1
        controls['extra_target_credit'] = reject(lambda: P.assemble(bad, native, caller, physical))
        bad = scalar(word)
        i, j = bad.index(['ADD', 213, 473, '1/2', 3]), bad.index(['ADD', 473, 100, '1', 716])
        bad[i], bad[j] = bad[j], bad[i]
        controls['helper_read_after_write'] = reject(lambda: trace_equivalent(scalar(base), bad))
    return P.encode({'status': STATUS, 'source': {'pr5660': '3b0b4f37d5005f381ed9cfcccebfbe1ab35c71aa'},
                     'helper_starts': starts, 'changed_charts': charts, 'local': local,
                     'local_paid_delta': paid_delta, 'helper_credit_delta': helper_delta,
                     'global_paid_delta': {r: 18*n for r, n in paid_delta.items()},
                     'transpose': transpose, 'native': native, 'physical': physical,
                     'complete_cost': price, 'controls': controls})


def main():
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='p104-r22-bytecode-') as cache:
        sys.pycache_prefix = cache; result = run()
    target = HERE/'RESULT.json'
    if args.write:
        target.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    else:
        need(read(target) == result, 'Exact saved complete receipt')
    print(result['status'], flush=True)


if __name__ == '__main__':
    main()
