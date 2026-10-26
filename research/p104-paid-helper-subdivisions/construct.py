"""Emit one exact paid helper split from a fixed chronological witness."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from hashlib import sha256

def need(ok, why):
    if not ok:
        raise ValueError(why)

def add(dst, src, q):
    for r, v in src.items():
        x = dst.get(r, Q(0))+q*v
        if x:
            dst[r] = x
        else:
            dst.pop(r, None)

def prefix(word, cut):
    n = word['public']+word['owned_fields']
    rows = [{r: Q(1)} for r in range(n)]
    charts = word['initial']+[None]*word['owned_fields']
    for e in word['program'][:cut]:
        op, r = e[:2]
        if op == 'MOVE':
            charts[r] = e[3]
        elif op in ('INIT', 'CLEAR'):
            rows[r] = {}
            charts[r] = e[2] if op == 'INIT' else None
        elif op == 'SCALE':
            rows[r] = {s: Q(e[-2])*v for s, v in rows[r].items()}
        else:
            if op == 'COPY':
                rows[r] = {}
            add(rows[r], rows[e[2]], Q(e[-2]))
    return rows, charts

def suffix(word, cut):
    n = word['public']
    cols = [{r: Q(1)} if r < n else {} for r in range(n+word['owned_fields'])]
    for e in reversed(word['program'][cut:]):
        op, r = e[:2]
        if op == 'MOVE':
            continue
        if op in ('INIT', 'CLEAR'):
            cols[r] = {}
        elif op == 'SCALE':
            cols[r] = {s: Q(e[-2])*v for s, v in cols[r].items()}
        else:
            add(cols[e[2]], cols[r], Q(e[-2]))
            if op == 'COPY':
                cols[r] = {}
    return cols

def construct(source, request):
    w = deepcopy(source)
    n = w['public']
    r, cut, move, f = (request[key] for key in ('role', 'cut', 'move', 'frame'))
    a, b = request['old_frames']
    need(w['program'][move] == ['MOVE', r, a, b, False], 'Exact paid source transport')
    for e in w['program'][min(cut, move):max(cut, move)]:
        need(e[0] == 'MOVE' or (e[1] != r and (e[0] not in ('ADD', 'COPY') or e[2] != r)),
             'No crossed scalar read or write')
    need(r < n and r not in set(w['live']) | {x+112 for x in w['live']}, 'Public helper pivot')
    for i in range(min(cut, move), max(cut, move)):
        need(i == move or w['program'][i][0] != 'MOVE' or w['program'][i][1] != r, 'No crossed second pivot transport')
    rows, charts = prefix(w, cut)
    need(charts[r] == (a if cut <= move else b), 'Actual source pivot chart')
    cols = suffix(w, cut)
    v = {int(s): Q(q) for s, q in request['cut_correction'].items()}
    u = {int(s): Q(q) for s, q in request['fresh_decoder'].items()}
    need(r not in v and r not in u, 'Cut corrections do not read the changed pivot')
    need(all(charts[s] == f for s in v.keys() | u.keys()), 'Actual cut source coframes')
    h, outrow = dict(rows[r]), dict(rows[r])
    for s, q in v.items():
        add(h, rows[s], q)
    for s, q in u.items():
        add(outrow, rows[s], q)
    d = cols[r]
    D = dict(d)
    for x in w['live']:
        if x in d:
            add(D, {x+112: d[x]}, Q(-1))
    g = dict(outrow)
    for x in w['live']:
        if x+112 in outrow:
            add(g, {x: outrow[x+112]}, Q(-1))
    add(g, {r: Q(1)}, Q(-1))
    first, last = w['initial'][r], w['final'][r]
    def span(basis):
        values = {0}
        for x in basis:
            values |= {x ^ y for y in tuple(values)}
        return values
    need(span(w['frames'][a]) < span(w['frames'][f]) < span(w['frames'][b]), 'Strict new native subdivisions')
    need(h.get(r) == outrow.get(r) == d.get(r) == D.get(r) == 1,
         'Unit helper pivots at both boundaries')
    need(sum(q*D.get(s, 0) for s, q in h.items()) == 1
         and sum(q*D.get(s, 0) for s, q in outrow.items()) == 1,
         'Both boundary functionals have unit pairing')
    need(all(s < n and w['initial'][s] == first for s in h.keys() | D.keys()),
         'Legal input basis coframe')
    need(r not in g and all(s < n and w['final'][s] == last for s in g.keys() | d.keys()),
         'Legal output normalization coframe')
    fresh = n
    w['public'] += 1
    w['initial'].append(f)
    w['final'].append(last)
    w['final'][r] = f
    w['schema'] = 'p104-paid-helper-subdivision/1'
    program = []
    for s, q in sorted(h.items()):
        if s != r:
            program.append(['ADD', r, s, str(-q), first])
    for s, q in sorted(D.items()):
        if s != r:
            program.append(['ADD', s, r, str(q), first])
    def rename(role, after):
        if role == r and after:
            return fresh
        return role+1 if role >= n else role
    for i, original in enumerate(source['program']):
        if i == cut:
            program.append(['MOVE', r, a, f, False])
            for s, q in sorted(v.items()):
                program.append(['ADD', r, rename(s, False), str(q), f])
            for s, q in sorted(u.items()):
                program.append(['ADD', fresh, rename(s, False), str(-q), f])
            program.append(['MOVE', fresh, f, b, False])
        if i == move:
            continue
        e = list(original)
        e[1] = rename(e[1], i >= cut)
        if e[0] in ('ADD', 'COPY'):
            e[2] = rename(e[2], i >= cut)
        program.append(e)
    for s, q in sorted(d.items()):
        if s != r:
            program.append(['ADD', s, fresh, str(-q), last])
    for s, q in sorted(g.items()):
        program.append(['ADD', fresh, s, str(q), last])
    w['program'] = program
    result = {'role': r, 'fresh': fresh, 'cut': cut, 'move': move,
              'cut_frame': f, 'old_frames': [a, b],
              'input_functional': h, 'output_functional': outrow,
              'suffix_column': d, 'input_direction': D,
              'terminal_correction': g,
              'cut_correction': v, 'fresh_decoder': u}
    return w, result

def audit(word):
    rows, charts = prefix(word, len(word['program']))
    n = word['public']
    for r, row in enumerate(rows):
        want = {r: Q(1)} if r < n else {}
        if r-112 in word['live']:
            want[r-112] = Q(1)
        need(row == want, 'Exact endpoint row '+str(r))
    need(charts[:n] == word['final'], 'Declared final frames')
    return n+word['owned_fields']

def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def canonical_hash(word):
    return sha256(json.dumps(word, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def build(source, plan):
    need(len(plan) == 17, 'Exact seventeen paid cuts')
    word = source
    witnesses = []
    for index, item in enumerate(plan):
        need(canonical_hash(word) == item['source_sha256'], 'Exact step source '+str(index))
        word, witness = construct(word, item['request'])
        need(canonical_hash(word) == item['word_sha256'], 'Exact step output '+str(index))
        witness['endpoint_columns'] = audit(word)
        witnesses.append(encode(witness))
    return word, witnesses
