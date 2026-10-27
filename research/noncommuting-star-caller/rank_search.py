"""Bounded low-rank discovery. Exact endpoint replay follows every model."""
import argparse
from functools import reduce
import json
from pathlib import Path
import random
import time

import z3
import screen as s


def run(g, d, rank_bound, trials, timeout, seed, output, total_rank=None):
    rng = random.Random(seed)
    order = [(j, 0) for j in range(g)] + [(j, 1) for j in reversed(range(g))]
    gaps = list(range(2*g+1))
    n, m = (g+1)*d, len(gaps)
    receipt = {'domain': {'g': g, 'd': d, 'rank_per_gap': rank_bound,
                'gaps': gaps, 'order': order, 'seed': seed,
                'trials': trials, 'timeout_ms_per_trial': timeout,
                'total_rank_bound': total_rank}, 'trials': []}
    for trial in range(trials):
        start = time.monotonic()
        a = [s.random_projection(d, rng) for _ in range(g)]
        data = s.prepare(a, order)
        if data is None:
            raise ValueError('Pure stage failure')
        pref, inv, rhs = data
        solver = z3.Solver()
        solver.set(timeout=timeout, random_seed=seed+trial)
        u = [[[z3.BitVec(f'u_{t}_{i}_{k}', 1) for k in range(rank_bound)]
              for i in range(n)] for t in range(m)]
        v = [[[z3.BitVec(f'v_{t}_{k}_{j}', 1) for j in range(n)]
              for k in range(rank_bound)] for t in range(m)]
        z = [[[z3.BitVec(f'z_{t}_{i}_{j}', 1) for j in range(n)]
              for i in range(n)] for t in range(m)]
        if total_rank is not None:
            used = [z3.If(z3.And(z3.Or(*[u[t][i][k] != 0 for i in range(n)]),
                                z3.Or(*[v[t][k][j] != 0 for j in range(n)])), 1, 0)
                    for t in range(m) for k in range(rank_bound)]
            solver.add(z3.Sum(used) <= total_rank)
        for t in range(m):
            for i in range(n):
                for j in range(n):
                    solver.add(z[t][i][j] == reduce(lambda x, y: x ^ y,
                        [u[t][i][k] & v[t][k][j] for k in range(rank_bound)]))
        neq = 0
        inconsistent = False
        for x in range(g+1):
            for y in range(g+1):
                equations, nv = s.block_equations([inv[t][x] for t in gaps],
                    [pref[t][y] for t in gaps], rhs[x][y], d)
                pivots = s.reduce_equations(equations, nv)
                if pivots is None:
                    inconsistent = True
                    break
                for row in pivots.values():
                    value = row >> nv
                    bits = row & ((1 << nv)-1)
                    terms = []
                    while bits:
                        low = bits & -bits
                        t, rem = divmod(low.bit_length()-1, d*d)
                        i, j = divmod(rem, d)
                        terms.append(z[t][x*d+i][y*d+j])
                        bits ^= low
                    solver.add(reduce(lambda x, y: x ^ y, terms) == value)
                    neq += 1
            if inconsistent:
                break
        answer = 'linear_inconsistent' if inconsistent else str(solver.check())
        record = {'trial': trial, 'projections': a, 'constraints': neq,
                  'status': answer, 'seconds': time.monotonic()-start}
        if answer == 'unknown':
            record['reason'] = solver.reason_unknown()
        if answer == 'sat':
            model = solver.model()
            thin = [[sum(model.eval(z[t][i][j], model_completion=True).as_long() << j
                         for j in range(n)) for i in range(n)] for t in range(m)]
            zero = (0,)*g
            for x in range(g+1):
                for y in range(g+1):
                    actual = {}
                    for t in gaps:
                        block = tuple((thin[t][x*d+i] >> (y*d)) & ((1 << d)-1)
                                      for i in range(d))
                        actual = s.padd(actual, s.pmul(s.pmul(inv[t][x],
                            {zero: block}), pref[t][y]))
                    if actual != rhs[x][y]:
                        raise ValueError('Model fails the complete endpoint')
            ranks = list(map(s.rank, thin))
            if max(ranks) > rank_bound:
                raise ValueError('Model exceeds its rank bound')
            if total_rank is not None and sum(ranks) > total_rank:
                raise ValueError('Model exceeds its total rank bound')
            record.update(thin=thin, ranks=ranks, total_rank=sum(ranks))
        receipt['trials'].append(record)
        output.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps({k: v for k, v in record.items() if k not in ['thin', 'projections']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--g', type=int, default=3)
    p.add_argument('--d', type=int, default=2)
    p.add_argument('--rank', type=int, default=1)
    p.add_argument('--trials', type=int, default=2)
    p.add_argument('--timeout', type=int, default=60000)
    p.add_argument('--seed', type=int, default=5747)
    p.add_argument('--total-rank', type=int)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.g, a.d, a.rank, a.trials, a.timeout, a.seed, a.output, a.total_rank)
