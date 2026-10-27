"""Build the complete selected P104 contact and helper graph.

The Lagrangian convention and exact distance follow Cinna Davis's contact
model and Dugongue's native-clifford-gauge-compiler package. The source
word and all prior attribution remain inputs. Apache-2.0.
"""
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            bit = row.bit_length()-1
            if bit not in basis:
                basis[bit] = row
                break
            row ^= basis[bit]
    return len(basis)


def nullspace(rows, width):
    rows = list(rows)
    pivots = []
    k = 0
    for j in range(width):
        found = next((i for i in range(k, len(rows)) if rows[i] >> j & 1), None)
        if found is None:
            continue
        rows[k], rows[found] = rows[found], rows[k]
        for i in range(len(rows)):
            if i != k and rows[i] >> j & 1:
                rows[i] ^= rows[k]
        pivots.append(j)
        k += 1
    return [(1 << j) | sum(((rows[i] >> j) & 1) << v for i, v in enumerate(pivots))
            for j in range(width) if j not in pivots]


def lagrangian(basis, width=9):
    return [x | (x << width) for x in basis] + [x << width for x in nullspace(basis, width)]


def transvection(basis, axis, width=9):
    mask = (1 << width)-1
    def paired(x):
        return (((x & mask) & (axis >> width)).bit_count() +
                ((x >> width) & (axis & mask)).bit_count()) & 1
    return [x ^ (axis if paired(x) else 0) for x in basis]


def distance(a, b, width=9):
    return rank(a+b)-width


def build(word):
    frames, program, n = word['frames'], word['program'], word['public']
    parent, frame = [], []
    def new(f):
        i = len(parent)
        parent.append(i)
        frame.append(f)
        return i
    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    def union(a, b):
        a, b = root(a), root(b)
        if frame[a] != frame[b]:
            raise ValueError('Source scalar contact has different frames')
        parent[b] = a
    current = [new(f) for f in word['initial']] + [None]*8
    entry = current[:n]
    edges, events, owned_nodes = [], [], []
    counts = Counter()
    for ei, op in enumerate(program):
        name, r = op[:2]
        if name == 'INIT':
            if current[r] is not None:
                raise ValueError('Owned field is already live')
            current[r] = new(op[2])
            owned_nodes.append(current[r])
        elif name == 'CLEAR':
            owned_nodes.append(current[r])
            events.append(current[r])
            current[r] = None
            continue
        elif name == 'MOVE':
            old, target = current[r], new(op[3])
            if frame[old] != op[2]:
                raise ValueError('Wrong MOVE entrance')
            grade = abs(len(frames[op[3]])-len(frames[op[2]]))
            edges.append((old, target, 1, grade, ei))
            counts[grade] += 1
            current[r] = target
        elif name in ('ADD', 'COPY'):
            if frame[current[r]] != op[-1] or frame[current[op[2]]] != op[-1]:
                raise ValueError('Wrong source contact chart')
            union(current[r], current[op[2]])
        elif name != 'SCALE':
            raise ValueError('Unknown source operation')
        events.append(current[r])
    helpers = []
    main = set(word['live']) | {r+112 for r in word['live']}
    for r in range(n):
        if frame[current[r]] != word['final'][r]:
            raise ValueError('Wrong source exit chart')
        if r not in main:
            grade = len(frames[word['final'][r]])-len(frames[word['initial'][r]])
            helpers.append((entry[r], current[r], -1, grade, r))
    roots = sorted({root(i) for i in range(len(parent))})
    ids = {a: i for i, a in enumerate(roots)}
    convert = lambda a: ids[root(a)]
    table = {kind: [[convert(a), convert(b), sign, grade, source]
                   for a, b, sign, grade, source in rows]
             for kind, rows in [('paid', edges), ('helpers', helpers)]}
    fixed = sorted({convert(x) for r in main for x in (entry[r], current[r])})
    nodes = [frame[a] for a in roots]
    lags = [lagrangian(b) for b in frames]
    for a, b, sign, grade, _ in table['paid'] + table['helpers']:
        if distance(lags[nodes[a]], lags[nodes[b]]) != grade:
            raise ValueError('Native chart distance differs from source grade')
    return {'nodes': nodes, 'fixed': fixed,
            'owned_contacts': sorted({convert(x) for x in owned_nodes}),
            'events': [convert(x) for x in events], 'lags': lags,
            'entry': [convert(x) for x in entry],
            'exit': [convert(x) for x in current[:n]],
            'paid_profile': dict(sorted(counts.items())),
            'helper_profile': dict(sorted(Counter(x[3] for x in helpers).items())), **table}


def labels(graph, axes):
    out = []
    for lag in graph['lags']:
        options = []
        for mask in range(1 << len(axes)):
            value = lag
            for i, axis in enumerate(axes):
                if mask >> i & 1:
                    value = transvection(value, axis)
            options.append(value)
        out.append(options)
    return out


def tables(graph, axes):
    options = labels(graph, axes)
    cache, result = {}, []
    for edge in graph['paid'] + graph['helpers']:
        a, b, sign, grade, source = edge[:5]
        fa, fb = graph['nodes'][a], graph['nodes'][b]
        key = (fa, fb)
        if key not in cache:
            cache[key] = [[distance(x, y) for y in options[fb]] for x in options[fa]]
        table = cache[key]
        if table[0][0] != grade:
            raise ValueError('Source tariff changed')
        result.append((a, b, sign, table, source))
    return result


def invoice(graph, axes, chosen):
    p, h = Counter(), Counter()
    for a, b, sign, costs, _ in tables(graph, axes):
        r = costs[chosen[a]][chosen[b]]
        if sign == -1 and r == 0:
            raise ValueError('An independent target has zero Gaussian rank')
        if r:
            (p if sign == 1 else h)[r] += 1
    return {'paid': dict(sorted(p.items())), 'helpers': dict(sorted(h.items()))}


if __name__ == '__main__':
    import sys
    source, output = map(Path, sys.argv[1:])
    raw = source.read_bytes()
    word = json.loads(gzip.decompress(raw))
    g = build(word)
    g['source_sha256'] = hashlib.sha256(raw).hexdigest()
    output.write_text(json.dumps(g, separators=(',', ':'))+'\n')
    print(json.dumps({k: g[k] for k in ['paid_profile', 'helper_profile']} |
                     {k: len(g[k]) for k in ['nodes', 'fixed', 'owned_contacts', 'paid', 'helpers']}))
