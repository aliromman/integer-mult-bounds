"""Reconstruct the selected word with eight arbitrary owned inputs.

This changes the source contract. Every COPY becomes ADD; INIT and CLEAR
only mark boundaries. The complete joint matrix is retained. It is not a
productive completion or a complete caller.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def add(dst, src, factor):
    for j, a in tuple(src.items()):
        b = dst.get(j, 0) + factor*a
        if b:
            dst[j] = b
        else:
            dst.pop(j, None)


def rank(rows):
    pivots = {}
    for source in rows:
        row = dict(source)
        while row:
            j = min(row)
            if j not in pivots:
                a = row[j]
                pivots[j] = {k: b/a for k, b in row.items()}
                break
            add(row, pivots[j], -row[j])
    return len(pivots)


def encoded(rows):
    return [[[j, str(a)] for j, a in sorted(row.items())] for row in rows]


def check():
    pin = next(p for p in json.loads((HERE/'SOURCE-PINS.json').read_text())
               if p['local_path'] == 'sources/selected-word.json.gz')
    raw = (HERE/pin['local_path']).read_bytes()
    need(sha256(raw).hexdigest() == pin['sha256'], 'Source bytes differ')
    w = json.loads(gzip.decompress(raw))
    public, count = w['public'], w['public']+8
    rows = [{i: Q(1)} for i in range(count)]
    starts, ends, copies = {}, {}, []
    paid = Counter()
    for event, e in enumerate(w['program']):
        op, r = e[:2]
        if op == 'INIT':
            need(r not in starts, 'Use one entrance per owned role')
            starts[r] = e[2]
        elif op == 'CLEAR':
            need(r not in ends, 'Use one exit per owned role')
            ends[r] = e[2]
        elif op == 'MOVE':
            paid[abs(len(w['frames'][e[2]])-len(w['frames'][e[3]]))] += 1
        elif op == 'SCALE':
            rows[r] = {j: Q(e[-2])*a for j, a in rows[r].items()}
        else:
            need(op in ('COPY', 'ADD'), 'Known scalar instruction')
            if op == 'COPY':
                copies.append({'event': event, 'owned': r, 'source': e[2],
                               'coefficient': str(Q(e[-2])), 'frame': e[-1]})
            add(rows[r], rows[e[2]], Q(e[-2]))
    need(sorted(starts) == sorted(ends) == list(range(public, count)),
         'Eight distinct owned boundaries')
    live = set(w['live'])
    # The old public block must be the exact original shear.
    for i, row in enumerate(rows[:public]):
        expected = {i: Q(1)}
        if i-112 in live:
            expected[i-112] = Q(1)
        need({j:a for j,a in row.items() if j < public} == expected,
             'Original public block')
    b = [{j:a for j,a in row.items() if j >= public} for row in rows[:public]]
    c = [{j:a for j,a in row.items() if j < public} for row in rows[public:]]
    d = [{j:a for j,a in row.items() if j >= public} for row in rows[public:]]
    tinvb = [dict(row) for row in b]
    for i in live:
        add(tinvb[i+112], b[i], -1)
    ctinv = [dict(row) for row in c]
    for row in ctinv:
        for i in live:
            add(row, {i:Q(1)}, -row.get(i+112, 0))
    schur = [dict(row) for row in d]
    for row, ci in zip(schur, c):
        for j, a in ci.items():
            add(row, tinvb[j], -a)
    need(rank(schur) == 8, 'Invertible joint operator')
    initial = w['initial'] + [starts[r] for r in range(public, count)]
    final = w['final'] + [ends[r] for r in range(public, count)]
    need(schur == d, 'The cross feedback is zero')
    for i, row in enumerate(d):
        need(set(row) == {public+i}, 'Diagonal owned block')
    for i, row in enumerate(tinvb):
        if row:
            need(not w['frames'][initial[i]], 'Preparation uses the zero frame')
    for row in ctinv:
        for i in row:
            need(len(w['frames'][final[i]]) == 9, 'Readout uses the full frame')
    # A separate reverse interpreter gives every column of the transpose.
    dual = [{i: Q(1)} for i in range(count)]
    for e in reversed(w['program']):
        op, r = e[:2]
        if op in ('MOVE', 'INIT', 'CLEAR'):
            continue
        if op == 'SCALE':
            dual[r] = {j: Q(e[-2])*a for j, a in dual[r].items()}
        else:
            add(dual[e[2]], dual[r], Q(e[-2]))
    transpose = [{} for _ in range(count)]
    for i, row in enumerate(rows):
        for j, a in row.items():
            transpose[j][i] = a
    need(dual == transpose, 'Independent complete transpose endpoint')
    main = live | {i+112 for i in live}
    support = []
    for j in range(public, count):
        affected = [i for i, row in enumerate(b) if j in row]
        support.append({'input':j, 'public_outputs':affected,
                        'main_outputs':[i for i in affected if i in main],
                        'helper_outputs':[i for i in affected if i not in main]})
    return {'status':'PASS_EXACT_JOINT_OWNED_MATRIX', 'new_kappa':False,
            'source_commit':pin['commit'], 'source_sha256':pin['sha256'],
            'contract':'COPY to ADD; omit zero initialization and erasure',
            'public':public, 'independent_inputs':count,
            'initial':initial, 'final':final, 'frames':w['frames'],
            'changed_copy_events':copies, 'paid':dict(sorted(paid.items())),
            'matrix':encoded(rows), 'B':encoded(b), 'C':encoded(c),
            'D':encoded(d), 'T_inverse_B':encoded(tinvb),
            'C_T_inverse':encoded(ctinv), 'schur':encoded(schur),
            'B_rank':rank(b), 'C_rank':rank(c), 'schur_rank':rank(schur),
            'transpose_coefficients_checked':count*count,
            'direct_schur_completion':{'added_paid':{'8':16},
                                       'added_targets':{'7':8},
                                       'linear_margin_change':-72},
            'dirty_input_support':support,
            'scope':'Exact scalar endpoint and literal boundary charts; no native caller or price gain'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = check()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','independent_inputs','B_rank','C_rank','schur_rank','new_kappa')}))
