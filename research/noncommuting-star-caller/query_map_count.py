"""Eliminate thin coefficients before a finite full-query frame search.

The test allows unrestricted thin ranks at specified gaps. It retains every
Laurent coefficient. The solver selects full queries from the exact finite
nilpotent list. A separate matrix replay checks each accepted word.
"""
import argparse
import itertools
import json
from pathlib import Path
import random
import time

import z3
import screen as s
from coupled_query_search import units, nilpotents, xor


def thin_elimination(equations, count):
    pivots, remaining = {}, []
    mask = (1 << count)-1
    for row in equations.values():
        while row & mask:
            low = row & -row
            pivot = low.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
        else:
            if row:
                remaining.append(row >> count)
    return remaining


def search(projections, order, gaps, timeout):
    g, d = len(projections), len(projections[0])
    pref, inv, _ = s.prepare(projections, order)
    full_units = [units(projections[j] if half == 0 else s.mxor(s.eye(d), projections[j]))
                  for j, half in order]
    nv_thin, nv_full = len(gaps)*d*d, len(order)*d*d
    nv = nv_thin+nv_full
    zero = (0,)*g
    solver = z3.Solver()
    solver.set(timeout=timeout)
    words = [z3.BitVec(f'full_{t}', d*d) for t in range(len(order))]
    bits = [z3.Extract(k, k, word) for word in words for k in range(d*d)]
    options = nilpotents(d)
    for word in words:
        solver.add(z3.Or(*[word == value for value in options]))
    rows = []
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
                    for mon, matrix in contribution.items():
                        for i, row in enumerate(matrix):
                            for j in range(d):
                                if row >> j & 1:
                                    key = mon, i, j
                                    eq[key] = eq.get(key, 0) ^ (1 << (nv_thin+t*d*d+kind))
            rows.extend(thin_elimination(eq, nv_thin))
    pivots = s.reduce_equations(dict(enumerate(rows)), nv_full)
    if pivots is None:
        return {'status': 'linear_inconsistent'}
    for row in pivots.values():
        solver.add(xor([bits[k] for k in range(nv_full) if row >> k & 1]) == row >> nv_full)
    answer = str(solver.check())
    result = {'status': answer, 'full_linear_constraints': len(pivots)}
    if answer == 'unknown':
        result['reason'] = solver.reason_unknown()
    if answer == 'sat':
        model = solver.model()
        chosen = [model.eval(word).as_long() for word in words]
        thin = [[0]*((g+1)*d) for _ in gaps]
        for x in range(g+1):
            for y in range(g+1):
                target = {zero: s.eye(d)} if x == y else {}
                for t, (leaf, _) in enumerate(order):
                    roles = (0, leaf+1)
                    if x in roles and y in roles:
                        block = (0,)*d
                        for kind in range(d*d):
                            if chosen[t] >> kind & 1:
                                block = s.mxor(block, full_units[t][kind][roles.index(x)][roles.index(y)])
                        target = s.padd(target, s.pmul(s.pmul(inv[t][x], {zero: block}), pref[t][y]))
                solution = s.solve_block([inv[t][x] for t in gaps], [pref[t][y] for t in gaps], target, d)
                if solution is None:
                    raise ValueError('Full-query model fails thin solvability')
                for k, block in enumerate(solution['blocks']):
                    for i, row in enumerate(block):
                        thin[k][x*d+i] |= row << (y*d)
        result.update(full_queries=chosen, thin=thin, thin_ranks=list(map(s.rank, thin)))
    return result


def run(args):
    rng = random.Random(args.seed)
    order = [(j, 0) for j in range(args.g)]+[(j, 1) for j in reversed(range(args.g))]
    receipt = {'scope': 'nested B=0 stages; all binary full-query frames; unrestricted thin ranks',
               'parameters': vars(args) | {'output': str(args.output)}, 'trials': []}
    for trial in range(args.trials):
        projections = [s.random_projection(args.d, rng) for _ in range(args.g)]
        record = {'trial': trial, 'projections': projections, 'systems': []}
        receipt['trials'].append(record)
        for gaps in itertools.combinations(range(2*args.g+1), args.thin):
            started = time.monotonic()
            result = search(projections, order, gaps, args.timeout)
            result.update(gaps=gaps, seconds=time.monotonic()-started)
            record['systems'].append(result)
            args.output.write_text(json.dumps(receipt, indent=2)+'\n')
            print(json.dumps({key: value for key, value in result.items() if key != 'thin'}), flush=True)
            if result['status'] == 'sat':
                break


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--g', type=int, default=5)
    p.add_argument('--d', type=int, choices=(2, 4), default=4)
    p.add_argument('--thin', type=int, default=2)
    p.add_argument('--trials', type=int, default=1)
    p.add_argument('--timeout', type=int, default=3000)
    p.add_argument('--seed', type=int, default=5747)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args())
