"""Expose a native transport position at every scalar contact.

The scalar order and coefficients are fixed. Extra transports have a cost.
"""
from collections import Counter, deque
import gzip
import hashlib
import json
from pathlib import Path

import gauge_graph as gg
from binary_subspaces import span, join


def build(word):
    n = word['public']
    frames = word['frames']
    nodes = list(word['initial'])
    entry = list(range(n))
    current = entry+[None]*8
    paid, owned, events = [], [], []
    def new(f):
        nodes.append(f)
        return len(nodes)-1
    def link(r, target, ei):
        old = current[r]
        grade = abs(len(frames[nodes[old]])-len(frames[nodes[target]]))
        paid.append([old, target, 1, grade, ei, r])
        current[r] = target
    for ei, op in enumerate(word['program']):
        name, r = op[:2]
        if name == 'INIT':
            if current[r] is not None:
                raise ValueError('Owned field is already live')
            current[r] = new(op[2])
            owned.append(current[r])
        elif name == 'CLEAR':
            owned.append(current[r])
            events.append(current[r])
            current[r] = None
            continue
        elif name == 'MOVE':
            if nodes[current[r]] != op[2]:
                raise ValueError('Wrong source transport')
            link(r, new(op[3]), ei)
        elif name in ('ADD', 'COPY'):
            s = op[2]
            if r == s or nodes[current[r]] != op[-1] or nodes[current[s]] != op[-1]:
                raise ValueError('Wrong scalar contact')
            target = new(op[-1])
            link(r, target, ei)
            link(s, target, ei)
        elif name != 'SCALE':
            raise ValueError('Unknown event')
        events.append(current[r])
    main = set(word['live']) | {r+112 for r in word['live']}
    helpers = []
    for r in range(n):
        if nodes[current[r]] != word['final'][r]:
            raise ValueError('Wrong public exit')
        if r not in main:
            a, b = entry[r], current[r]
            grade = len(frames[nodes[b]])-len(frames[nodes[a]])
            helpers.append([a, b, -1, grade, r])
    profile = Counter(e[3] for e in paid if e[3])
    return {'nodes': nodes, 'entry': entry, 'exit': current[:n],
            'fixed': sorted({v for r in main for v in (entry[r], current[r])}),
            'owned_contacts': sorted(set(owned)), 'events': events,
            'paid': paid, 'helpers': helpers,
            'paid_profile': dict(sorted(profile.items())),
            'helper_profile': dict(sorted(Counter(e[3] for e in helpers).items())),
            'lags': [gg.lagrangian(b) for b in frames]}


def bounds(word, graph):
    spaces = [span(b) for b in word['frames']]
    original = [spaces[f] for f in graph['nodes']]
    fixed = set(graph['fixed'])
    directed = []
    for a, b, _, _, _, role in graph['paid']:
        if role >= word['public']:
            a, b = b, a
        if not original[a] <= original[b]:
            raise ValueError('Wrong source edge direction')
        directed.append((a, b))
    incoming = [[] for _ in original]
    outgoing = [[] for _ in original]
    for a, b in directed:
        outgoing[a].append(b)
        incoming[b].append(a)
    lower = [s if i in fixed else frozenset({0}) for i, s in enumerate(original)]
    upper = [s if i in fixed else frozenset(range(512)) for i, s in enumerate(original)]
    for values, neighbors, combine in ((lower, outgoing, join),
                                        (upper, incoming, lambda a, b: a & b)):
        queue = deque(fixed)
        pending = set(fixed)
        while queue:
            a = queue.popleft()
            pending.remove(a)
            for b in neighbors[a]:
                value = combine(values[a], values[b])
                if value != values[b]:
                    values[b] = value
                    if b not in pending:
                        queue.append(b)
                        pending.add(b)
    if any(not lo <= old <= hi for lo, old, hi in zip(lower, original, upper)):
        raise ValueError('Bounds exclude the source')
    return original, lower, upper, directed, fixed


if __name__ == '__main__':
    import sys
    source, output = map(Path, sys.argv[1:])
    raw = source.read_bytes()
    word = json.loads(gzip.decompress(raw))
    graph = build(word)
    graph['source_sha256'] = hashlib.sha256(raw).hexdigest()
    original, lower, upper, directed, fixed = bounds(word, graph)
    graph['intervals'] = [{'lower': sorted(lo), 'upper': sorted(hi)} for lo, hi in zip(lower, upper)]
    data = (json.dumps(graph, separators=(',', ':'))+'\n').encode()
    output.write_bytes(gzip.compress(data, mtime=0) if output.suffix == '.gz' else data)
    print(json.dumps({'nodes': len(original), 'edges': len(directed),
                      'paid': graph['paid_profile'], 'helpers': graph['helper_profile'],
                      'gaps': dict(sorted(Counter((len(hi)//len(lo)).bit_length()-1
                                                for lo, hi in zip(lower, upper)).items()))}, indent=2))
