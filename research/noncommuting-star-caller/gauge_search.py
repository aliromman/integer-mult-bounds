"""Bounded multi-label contact-chart discovery with the full signed bill.

Floating MILP output proposes an assignment. It does not certify kappa.
The native word and its exact compiler obligations remain separate.
"""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array

import gauge_graph as gg


def canonical(rows):
    rows = list(rows)
    k = 0
    for j in range(18):
        found = next((i for i in range(k, len(rows)) if rows[i] >> j & 1), None)
        if found is None:
            continue
        rows[k], rows[found] = rows[found], rows[k]
        for i in range(len(rows)):
            if i != k and rows[i] >> j & 1:
                rows[i] ^= rows[k]
        k += 1
    return tuple(rows[:k])


def score(profile, saving):
    t = 1-saving
    return (sum(int(r)**t*n for r, n in profile['paid'].items()) -
            sum(int(r)**t*n for r, n in profile['helpers'].items()))


def search(graph, axes, saving, seconds, freeze_owned=False):
    start = time.monotonic()
    options = gg.labels(graph, axes)
    fixed = set(graph['fixed'])
    if freeze_owned:
        fixed.update(graph['owned_contacts'])
    allowed = []
    for v, f in enumerate(graph['nodes']):
        seen, keep = set(), []
        for label, lag in enumerate(options[f]):
            if v in fixed and label != 0:
                continue
            key = canonical(lag)
            if key not in seen:
                seen.add(key)
                keep.append(label)
        allowed.append(keep)
    t = 1-saving
    edges = {}
    for a, b, sign, ranks, _ in gg.tables(graph, axes):
        if a > b:
            a, b = b, a
            ranks = list(map(list, zip(*ranks)))
        if a == b:
            if any(ranks[l][l] for l in allowed[a]):
                raise ValueError('A self edge has nonzero grade')
            continue
        data = edges.setdefault((a, b), {'cost': defaultdict(float), 'forbidden': set()})
        for i in allowed[a]:
            for j in allowed[b]:
                r = ranks[i][j]
                data['cost'][i, j] += sign*(r**t if r else 0)
                if sign < 0 and r == 0:
                    data['forbidden'].add((i, j))
    objective, integrality, low, high = [], [], [], []
    def variable(cost, binary):
        k = len(objective)
        objective.append(cost)
        integrality.append(int(binary))
        low.append(0)
        high.append(1)
        return k
    x = [{l: variable(0, True) for l in ls} for ls in allowed]
    ri, ci, val, lower, upper = [], [], [], [], []
    def constraint(terms, lo, hi):
        r = len(lower)
        for c, v in terms:
            ri.append(r)
            ci.append(c)
            val.append(v)
        lower.append(lo)
        upper.append(hi)
    for node in x:
        constraint([(v, 1) for v in node.values()], 1, 1)
    constant = 0
    for (a, b), data in edges.items():
        base = data['cost'][0, 0]
        constant += base
        y = {(i, j): variable(data['cost'][i, j]-base, False)
             for i in allowed[a] for j in allowed[b] if (i, j) not in data['forbidden']}
        for i in allowed[a]:
            constraint([(k, 1) for (l, _), k in y.items() if l == i] + [(x[a][i], -1)], 0, 0)
        for j in allowed[b]:
            constraint([(k, 1) for (_, l), k in y.items() if l == j] + [(x[b][j], -1)], 0, 0)
    matrix = coo_array((np.array(val, dtype=float),
                        (np.array(ri, dtype=np.int32), np.array(ci, dtype=np.int32))),
                       shape=(len(lower), len(objective))).tocsc()
    result = milp(np.array(objective), integrality=np.array(integrality, dtype=np.int32),
                  bounds=Bounds(low, high), constraints=LinearConstraint(matrix, lower, upper),
                  options={'time_limit': seconds, 'mip_rel_gap': 1e-8})
    out = {'axes': axes, 'saving': saving, 'freeze_owned': freeze_owned,
           'variables': len(objective), 'binary_variables': sum(integrality),
           'constraints': len(lower), 'status': int(result.status), 'message': result.message,
           'elapsed': time.monotonic()-start,
           'mip_dual_bound': getattr(result, 'mip_dual_bound', None)}
    if result.x is not None:
        chosen = []
        for node in x:
            label = max(node, key=lambda l: result.x[node[l]])
            if abs(result.x[node[label]]-1) > 1e-6:
                raise ValueError('No integral assignment')
            chosen.append(label)
        if any(chosen[v] for v in fixed):
            raise ValueError('A fixed boundary changed')
        profile = gg.invoice(graph, axes, chosen)
        baseline = {'paid': graph['paid_profile'], 'helpers': graph['helper_profile']}
        delta = score(profile, saving)-score(baseline, saving)
        if abs(delta-result.fun) > 1e-6:
            raise ValueError('Recomputed whole invoice differs from MILP objective')
        out.update(chosen=chosen, changed_nodes=sum(bool(x) for x in chosen),
                   profile=profile, signed_cost_delta=delta)
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--graph', type=Path, required=True)
    p.add_argument('--axes', type=int, nargs='+', required=True)
    p.add_argument('--saving', type=float, default=.0034)
    p.add_argument('--seconds', type=float, default=30)
    p.add_argument('--freeze-owned', action='store_true')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    graph_bytes = a.graph.read_bytes()
    if a.graph.suffix == '.gz':
        graph_bytes = gzip.decompress(graph_bytes)
    result = search(json.loads(graph_bytes), a.axes, a.saving, a.seconds, a.freeze_owned)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['chosen', 'profile']}, indent=2))
