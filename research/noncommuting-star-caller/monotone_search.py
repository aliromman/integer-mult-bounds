"""Search every monotone frame on the fixed public-boundary contact graph.

The solver output is a discovery result. It is not an exact optimality proof.
"""
from collections import Counter
import argparse
import gzip
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array

import gauge_graph as gg


from binary_subspaces import span, basis, join, all_between


def intervals(word, graph, public_fixed=True):
    spaces = [span(x) for x in word['frames']]
    original = [spaces[f] for f in graph['nodes']]
    fixed = set(graph['entry']+graph['exit']) if public_fixed else set(graph['fixed'])
    directed = []
    for a, b, _, _, ei in graph['paid']:
        if word['program'][ei][1] >= word['public']:
            a, b = b, a
        if not original[a] < original[b]:
            raise ValueError('Source edge is not strict in the stated direction')
        directed.append((a, b))
    lower = [s if i in fixed else frozenset({0}) for i, s in enumerate(original)]
    upper = [s if i in fixed else frozenset(range(512)) for i, s in enumerate(original)]
    order = sorted(directed, key=lambda e: len(original[e[0]]))
    for a, b in order:
        lower[b] = join(lower[a], lower[b])
    for a, b in reversed(order):
        upper[a] &= upper[b]
    if any(not lo <= old <= hi for lo, old, hi in zip(lower, original, upper)):
        raise ValueError('Propagation rejected the source assignment')
    return original, lower, upper, directed, fixed


def search(word, graph, saving, seconds, cap=None, free_helpers=False, full_palette=False):
    start = time.monotonic()
    if full_palette:
        from full_contact_graph import bounds
        original, lower, upper, directed, fixed = bounds(word, graph)
        free_helpers = True
    else:
        original, lower, upper, directed, fixed = intervals(word, graph, not free_helpers)
    if free_helpers:
        # A component with one paid edge and the same target has zero cost.
        incident = Counter(v for e in graph['paid']+graph['helpers'] for v in e[:2])
        paid_pairs = Counter(tuple(e[:2]) for e in graph['paid'])
        for a, b, _, _, _ in graph['helpers']:
            if paid_pairs[a, b] == 1 and incident[a] == incident[b] == 2:
                for v in (a, b):
                    lower[v] = upper[v] = original[v]
                    fixed.add(v)
    if full_palette:
        candidates = [{lo, old, hi} for lo, old, hi in zip(lower, original, upper)]
        for a, b in directed:
            if lower[a] <= original[b] <= upper[a]:
                candidates[a].add(original[b])
            if lower[b] <= original[a] <= upper[b]:
                candidates[b].add(original[a])
        choices = [all_between(lo, hi) if len(hi)//len(lo) <= 4
                   else sorted(ls, key=lambda s: (len(s), sorted(s)))
                   for lo, hi, ls in zip(lower, upper, candidates)]
    else:
        choices = [all_between(lo, hi) for lo, hi in zip(lower, upper)]
    dimensions = [[len(s).bit_length()-1 for s in ls] for ls in choices]
    price = (lambda r: min(r, cap)) if cap is not None else (lambda r: r**(1-saving) if r else 0)
    objective, integers = [], []
    def variable(cost, binary):
        i = len(objective)
        objective.append(cost)
        integers.append(int(binary))
        return i
    x = [[variable(0, True) for _ in ls] for ls in choices]
    ri, ci, values, low, high = [], [], [], [], []
    def row(terms, lo=0, hi=0):
        r = len(low)
        for c, v in terms:
            ri.append(r)
            ci.append(c)
            values.append(v)
        low.append(lo)
        high.append(hi)
    for ls in x:
        row([(k, 1) for k in ls], 1, 1)
    edges = [(a, b, count, False) for (a, b), count in Counter(directed).items()]
    if free_helpers:
        edges += [(a, b, -count, True)
                  for (a, b), count in Counter(tuple(e[:2]) for e in graph['helpers']).items()]
    for a, b, count, strict in edges:
        y = {}
        for i, sa in enumerate(choices[a]):
            for j, sb in enumerate(choices[b]):
                if sa <= sb and (not strict or sa != sb):
                    grade = dimensions[b][j]-dimensions[a][i]
                    y[i, j] = variable(count*price(grade), False)
        for i in range(len(x[a])):
            row([(k, 1) for (ii, _), k in y.items() if ii == i]+[(x[a][i], -1)])
        for j in range(len(x[b])):
            row([(k, 1) for (_, jj), k in y.items() if jj == j]+[(x[b][j], -1)])
    matrix = coo_array((np.array(values, dtype=float),
                        (np.array(ri, dtype=np.int32), np.array(ci, dtype=np.int32))),
                       shape=(len(low), len(objective))).tocsc()
    result = milp(np.array(objective), integrality=np.array(integers, dtype=np.int32),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, low, high),
                  options={'time_limit': seconds, 'mip_rel_gap': 1e-10})
    scope = 'all monotone frames; main boundaries fixed' if free_helpers else 'all monotone frames; all public boundaries fixed'
    if full_palette:
        scope = 'all scalar contacts; bounded monotone frame palette; main boundaries fixed'
    out = {'saving': saving, 'cap': cap, 'scope': scope,
           'status': int(result.status), 'message': result.message,
           'elapsed_seconds': time.monotonic()-start,
           'interval_dimensions': dict(sorted(Counter((len(hi)//len(lo)).bit_length()-1
                                                       for lo, hi in zip(lower, upper)).items())),
           'choices': sum(map(len, choices)), 'variables': len(objective),
           'constraints': len(low)}
    if result.x is not None:
        chosen = [max(range(len(ls)), key=lambda j: result.x[ls[j]]) for ls in x]
        if any(abs(result.x[ls[j]]-1) > 1e-6 for ls, j in zip(x, chosen)):
            raise ValueError('Solver did not give an integral assignment')
        selected = [ls[j] for ls, j in zip(choices, chosen)]
        if any(selected[i] != original[i] for i in fixed):
            raise ValueError('Fixed boundary changed')
        if any(not selected[a] <= selected[b] for a, b in directed):
            raise ValueError('Selected edge is not nested')
        paid = Counter((len(selected[b])//len(selected[a])).bit_length()-1 for a, b in directed)
        paid.pop(0, None)
        helpers = Counter((len(selected[b])//len(selected[a])).bit_length()-1
                          for a, b, _, _, _ in graph['helpers'])
        if helpers.get(0):
            raise ValueError('An independent target is not productive')
        cost = sum(n*price(r) for r, n in paid.items())
        baseline = sum(n*price(int(r)) for r, n in graph['paid_profile'].items())
        if free_helpers:
            cost -= sum(n*price(r) for r, n in helpers.items())
            baseline -= sum(n*price(int(r)) for r, n in graph['helper_profile'].items())
        if abs(cost-result.fun) > 1e-6:
            raise ValueError('Complete invoice differs from the solver objective')
        out.update(paid=dict(sorted(paid.items())), helpers=dict(sorted(helpers.items())), paid_cost_delta=cost-baseline,
                   changed=[{'node': i, 'old': basis(original[i]), 'new': basis(s)}
                            for i, s in enumerate(selected) if s != original[i]],
                   lower_bound_delta=float(result.mip_dual_bound)-baseline)
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--word', type=Path, required=True)
    p.add_argument('--graph', type=Path, required=True)
    p.add_argument('--saving', type=float, default=.0034)
    p.add_argument('--seconds', type=float, default=60)
    p.add_argument('--cap', type=int)
    p.add_argument('--free-helpers', action='store_true')
    p.add_argument('--full-palette', action='store_true')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    word = json.loads(gzip.decompress(a.word.read_bytes()))
    graph_bytes = a.graph.read_bytes()
    if a.graph.suffix == '.gz':
        graph_bytes = gzip.decompress(graph_bytes)
    graph = json.loads(graph_bytes)
    result = search(word, graph, a.saving, a.seconds, a.cap, a.free_helpers, a.full_palette)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'changed'}, indent=2))
