"""Search all binary full-query frames on a four-dimensional image.

The stage projectors and their nested order are fixed in each trial.
Thin maps can occur at every gap. Each accepted model gets a full Laurent
coefficient replay. An unsuccessful solver run is only a bounded screen.
"""
import argparse
from functools import reduce
import itertools
import json
from pathlib import Path
import random
import time

import z3
import screen as s


def xor(values):
    return reduce(lambda a, b: a ^ b, values, z3.BitVecVal(0, 1))


def image_coordinates(a):
    """Return column and row factors B,C with BC=A and CB=I."""
    d = len(a)
    columns = []
    for column in s.mt(a):
        if s.rank(columns+[column]) > len(columns):
            columns.append(column)
    q = len(columns)
    lookup = {reduce(int.__xor__, [columns[j] for j in range(q) if m >> j & 1], 0): m
              for m in range(1 << q)}
    coefficients = [lookup[column] for column in s.mt(a)]
    b = tuple(sum(((columns[j] >> i) & 1) << j for j in range(q)) for i in range(d))
    c = tuple(sum(((coefficients[j] >> i) & 1) << j for j in range(d)) for i in range(q))
    if s.mm(b, c) != a or s.mm(c, b) != s.eye(q):
        raise ValueError('Invalid image coordinates')
    return b, c


def units(a):
    b, c = image_coordinates(a)
    q, d = len(c), len(a)
    out = []
    for i in range(2*q):
        for j in range(2*q):
            rows = tuple(c[j % q] if (b[k] >> (i % q)) & 1 else 0 for k in range(d))
            blocks = [[(0,)*d for _ in range(2)] for _ in range(2)]
            blocks[i//q][j//q] = rows
            out.append(blocks)
    return out


def nilpotents(d):
    result = []
    for rows in itertools.product(range(1 << d), repeat=d):
        if s.rank(rows) == d//2 and not any(s.mm(rows, rows)):
            result.append(sum(row << (i*d) for i, row in enumerate(rows)))
    if d == 4 and len(result) != 210:
        raise ValueError('Wrong four-dimensional nilpotent count')
    return result


def run(args):
    g, d, cap = args.g, 4, args.rank_per_gap
    n = (g+1)*d
    order = [(j, 0) for j in range(g)]+[(j, 1) for j in reversed(range(g))]
    gaps = list(range(2*g+1)) if args.gaps is None else [int(x) for x in args.gaps.split(',')]
    if len(set(gaps)) != len(gaps) or any(t < 0 or t > 2*g for t in gaps):
        raise ValueError('Invalid thin gap list')
    options = nilpotents(d)
    rng = random.Random(args.seed)
    receipt = {'scope': 'nested B=0 stage projections; all 210 rank-two full-query nilpotents',
               'parameters': vars(args) | {'output': str(args.output)}, 'trials': []}
    for trial in range(args.trials):
        started = time.monotonic()
        projections = [s.random_projection(d, rng) for _ in range(g)]
        pref, inv, _ = s.prepare(projections, order)
        full_units = [units(projections[j] if half == 0 else s.mxor(s.eye(d), projections[j]))
                      for j, half in order]
        solver = z3.Solver()
        solver.set(timeout=args.timeout, random_seed=args.seed+trial)
        full_words = [z3.BitVec(f'f_{t}', d*d) for t in range(2*g)]
        full = [[z3.Extract(k, k, w) for k in range(d*d)] for w in full_words]
        for w in full_words:
            solver.add(z3.Or(*[w == value for value in options]))
        u = [[[z3.BitVec(f'u_{t}_{k}_{i}', 1) for i in range(n)] for k in range(cap)] for t in gaps]
        v = [[[z3.BitVec(f'v_{t}_{k}_{i}', 1) for i in range(n)] for k in range(cap)] for t in gaps]
        z = [[[z3.BitVec(f'z_{t}_{i}_{j}', 1) for j in range(n)] for i in range(n)] for t in gaps]
        for t in range(len(gaps)):
            for i in range(n):
                for j in range(n):
                    solver.add(z[t][i][j] == xor([u[t][k][i] & v[t][k][j] for k in range(cap)]))
        used = [z3.If(z3.And(z3.Or(*[x != 0 for x in u[t][k]]),
                            z3.Or(*[x != 0 for x in v[t][k]])), 1, 0)
                for t in range(len(gaps)) for k in range(cap)]
        solver.add(z3.Sum(used) <= args.thin_rank)
        nv_thin, nv = len(gaps)*d*d, (len(gaps)+2*g)*d*d
        zero = (0,)*g
        neq = 0
        for x in range(g+1):
            for y in range(g+1):
                target = {zero: s.eye(d)} if x == y else {}
                eq, old_nv = s.block_equations([inv[t][x] for t in gaps],
                                               [pref[t][y] for t in gaps], target, d)
                eq = {key: (row & ((1 << old_nv)-1)) | ((row >> old_nv) << nv) for key, row in eq.items()}
                for t, (leaf, _) in enumerate(order):
                    roles = (0, leaf+1)
                    if x not in roles or y not in roles:
                        continue
                    for kind, blocks in enumerate(full_units[t]):
                        block = blocks[roles.index(x)][roles.index(y)]
                        contribution = s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y])
                        for mon, rows in contribution.items():
                            for i, row in enumerate(rows):
                                for j in range(d):
                                    if row >> j & 1:
                                        key = mon, i, j
                                        eq[key] = eq.get(key, 0) ^ (1 << (nv_thin+t*d*d+kind))
                pivots = s.reduce_equations(eq, nv)
                if pivots is None:
                    solver.add(False)
                    continue
                for row in pivots.values():
                    terms, bits = [], row & ((1 << nv)-1)
                    while bits:
                        bit = bits & -bits
                        index = bit.bit_length()-1
                        if index < nv_thin:
                            t, rem = divmod(index, d*d)
                            i, j = divmod(rem, d)
                            terms.append(z[t][x*d+i][y*d+j])
                        else:
                            t, kind = divmod(index-nv_thin, d*d)
                            terms.append(full[t][kind])
                        bits ^= bit
                    solver.add(xor(terms) == row >> nv)
                    neq += 1
        answer = str(solver.check())
        record = {'trial': trial, 'projections': projections, 'status': answer,
                  'constraints': neq, 'seconds': time.monotonic()-started}
        if answer == 'unknown':
            record['reason'] = solver.reason_unknown()
        if answer == 'sat':
            model = solver.model()
            chosen = [model.eval(w).as_long() for w in full_words]
            thin = [tuple(sum(model.eval(z[t][i][j], model_completion=True).as_long() << j
                              for j in range(n)) for i in range(n)) for t in range(len(gaps))]
            for x in range(g+1):
                for y in range(g+1):
                    actual = {}
                    for ti, t in enumerate(gaps):
                        block = tuple((thin[ti][x*d+i] >> (y*d)) & ((1 << d)-1) for i in range(d))
                        actual = s.padd(actual, s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y]))
                    for t, (leaf, _) in enumerate(order):
                        roles = (0, leaf+1)
                        if x in roles and y in roles:
                            block = (0,)*d
                            for kind in range(d*d):
                                if chosen[t] >> kind & 1:
                                    block = s.mxor(block, full_units[t][kind][roles.index(x)][roles.index(y)])
                            actual = s.padd(actual, s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y]))
                    if actual != ({zero: s.eye(d)} if x == y else {}):
                        raise ValueError('Model fails the complete Laurent endpoint')
            ranks = list(map(s.rank, thin))
            if max(ranks) > cap or sum(ranks) > args.thin_rank:
                raise ValueError('Model exceeds the thin rank bound')
            record.update(full_queries=chosen, thin=thin, thin_ranks=ranks, thin_rank=sum(ranks))
        receipt['trials'].append(record)
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps({key: value for key, value in record.items() if key not in ['projections', 'thin']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--g', type=int, default=3)
    p.add_argument('--thin-rank', type=int, default=8)
    p.add_argument('--rank-per-gap', type=int, default=2)
    p.add_argument('--gaps', help='Comma-separated gap indices; default is every gap')
    p.add_argument('--trials', type=int, default=1)
    p.add_argument('--timeout', type=int, default=45000)
    p.add_argument('--seed', type=int, default=5747)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args())
