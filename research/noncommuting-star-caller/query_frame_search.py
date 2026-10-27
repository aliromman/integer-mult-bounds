"""Vary full-query nilpotents in a paired-stage binary caller.

Each full image has dimension two. Its query can be any of the three
rank-one nilpotents over F2. Every stage coefficient is checked.
"""
import argparse
from functools import reduce
import json
from pathlib import Path
import random
import time

import z3
import screen as s


def xor(terms):
    return reduce(lambda a, b: a ^ b, terms, z3.BitVecVal(0, 1))


def choices(a):
    zero = (0,)*len(a)
    return (((zero, a), (zero, zero)),
            ((zero, zero), (a, zero)),
            ((a, a), (a, a)))


def run(g, thin_bound, trials, timeout, seed, output, rank_per_gap=1):
    d = 2
    rng = random.Random(seed)
    order = [(j, 0) for j in range(g)]+[(j, 1) for j in reversed(range(g))]
    gaps = list(range(2*g+1))
    n, m = (g+1)*d, len(gaps)
    receipt = {'scope': 'nested B=0 projectors; all rank-one full query frames; rank-factored thin maps at all gaps',
               'g': g, 'thin_rank_bound': thin_bound, 'seed': seed,
               'rank_per_gap': rank_per_gap, 'timeout_ms': timeout, 'trials': []}
    for trial in range(trials):
        start = time.monotonic()
        projections = [s.random_projection(d, rng) for _ in range(g)]
        pref, inv, _ = s.prepare(projections, order)
        query_options = [choices(projections[j] if half == 0 else s.mxor(s.eye(d), projections[j]))
                         for j, half in order]
        solver = z3.Solver()
        solver.set(timeout=timeout, random_seed=seed+trial)
        u = [[[z3.BitVec(f'u_{t}_{k}_{i}', 1) for i in range(n)]
              for k in range(rank_per_gap)] for t in gaps]
        v = [[[z3.BitVec(f'v_{t}_{k}_{i}', 1) for i in range(n)]
              for k in range(rank_per_gap)] for t in gaps]
        z = [[[xor([u[t][k][i] & v[t][k][j] for k in range(rank_per_gap)])
               for j in range(n)] for i in range(n)] for t in gaps]
        active = [z3.If(z3.And(z3.Or(*[x != 0 for x in u[t][k]]),
                               z3.Or(*[x != 0 for x in v[t][k]])), 1, 0)
                  for t in gaps for k in range(rank_per_gap)]
        solver.add(z3.Sum(active) <= thin_bound)
        full = [[z3.BitVec(f'f_{t}_{k}', 1) for k in range(3)] for t in range(2*g)]
        for fs in full:
            solver.add(z3.Sum([z3.If(x != 0, 1, 0) for x in fs]) == 1)
        nv_thin = m*d*d
        nv = nv_thin+6*g
        zero = (0,)*g
        neq = 0
        linear_failure = False
        for x in range(g+1):
            for y in range(g+1):
                target = {zero: s.eye(d)} if x == y else {}
                eq, old_nv = s.block_equations([inv[t][x] for t in gaps],
                                               [pref[t][y] for t in gaps], target, d)
                eq = {k: (v & ((1 << old_nv)-1)) | ((v >> old_nv) << nv) for k, v in eq.items()}
                for t, (leaf, _) in enumerate(order):
                    roles = (0, leaf+1)
                    if x not in roles or y not in roles:
                        continue
                    a, b = roles.index(x), roles.index(y)
                    for kind in range(3):
                        block = query_options[t][kind][a][b]
                        contribution = s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y])
                        for mon, rows in contribution.items():
                            for i, row in enumerate(rows):
                                for j in range(d):
                                    if row >> j & 1:
                                        key = mon, i, j
                                        eq[key] = eq.get(key, 0) ^ (1 << (nv_thin+3*t+kind))
                pivots = s.reduce_equations(eq, nv)
                if pivots is None:
                    linear_failure = True
                    break
                for row in pivots.values():
                    value = row >> nv
                    bits = row & ((1 << nv)-1)
                    terms = []
                    while bits:
                        bit = bits & -bits
                        index = bit.bit_length()-1
                        if index < nv_thin:
                            t, rem = divmod(index, d*d)
                            i, j = divmod(rem, d)
                            terms.append(z[t][x*d+i][y*d+j])
                        else:
                            t, kind = divmod(index-nv_thin, 3)
                            terms.append(full[t][kind])
                        bits ^= bit
                    solver.add(xor(terms) == value)
                    neq += 1
            if linear_failure:
                break
        answer = 'linear_inconsistent' if linear_failure else str(solver.check())
        result = {'trial': trial, 'projections': projections, 'order': order,
                  'constraints': neq, 'status': answer, 'seconds': time.monotonic()-start}
        if answer == 'unknown':
            result['reason'] = solver.reason_unknown()
        if answer == 'sat':
            model = solver.model()
            selected = [next(k for k in range(3) if model.eval(fs[k]).as_long()) for fs in full]
            thin = [tuple(sum(model.eval(z[t][i][j], model_completion=True).as_long() << j
                              for j in range(n)) for i in range(n)) for t in gaps]
            for x in range(g+1):
                for y in range(g+1):
                    actual = {}
                    for t in gaps:
                        block = tuple((thin[t][x*d+i] >> (y*d)) & ((1 << d)-1) for i in range(d))
                        actual = s.padd(actual, s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y]))
                    for t, (leaf, _) in enumerate(order):
                        roles = (0, leaf+1)
                        if x in roles and y in roles:
                            block = query_options[t][selected[t]][roles.index(x)][roles.index(y)]
                            actual = s.padd(actual, s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y]))
                    expected = {zero: s.eye(d)} if x == y else {}
                    if actual != expected:
                        raise ValueError('Model fails a complete Laurent block')
            ranks = list(map(s.rank, thin))
            if max(ranks) > rank_per_gap or sum(ranks) > thin_bound:
                raise ValueError('Model exceeds the thin rank bound')
            result.update(full_types=selected, thin=thin, thin_ranks=ranks, thin_rank=sum(ranks))
        receipt['trials'].append(result)
        output.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps({k: v for k, v in result.items() if k not in ['projections', 'order', 'thin']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--g', type=int, default=5)
    p.add_argument('--thin-rank', type=int, default=3)
    p.add_argument('--trials', type=int, default=3)
    p.add_argument('--timeout', type=int, default=45000)
    p.add_argument('--seed', type=int, default=5747)
    p.add_argument('--rank-per-gap', type=int, default=1)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.g, a.thin_rank, a.trials, a.timeout, a.seed, a.output, a.rank_per_gap)
