"""Exact Laurent endpoint screens for block-diagonal paired stages.

Discovery only. The field is F2. Thin maps have no rank restriction here.
All scalar stage variables are independent. Apache-2.0.
"""
import argparse
import itertools
import json
import random
import time
from pathlib import Path


def eye(d):
    return tuple(1 << i for i in range(d))


def mxor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def mm(a, b):
    out = []
    for row in a:
        z = 0
        while row:
            low = row & -row
            z ^= b[low.bit_length() - 1]
            row ^= low
        out.append(z)
    return tuple(out)


def mt(a):
    return tuple(sum(((row >> j) & 1) << i for i, row in enumerate(a))
                 for j in range(len(a)))


def rank(a):
    pivots = {}
    for row in a:
        while row:
            j = row.bit_length() - 1
            if j not in pivots:
                pivots[j] = row
                break
            row ^= pivots[j]
    return len(pivots)


def inverse(a):
    d = len(a)
    rows = [r | (1 << (d+i)) for i, r in enumerate(a)]
    for j in range(d):
        k = next(k for k in range(j, d) if (rows[k] >> j) & 1)
        rows[j], rows[k] = rows[k], rows[j]
        for i in range(d):
            if i != j and (rows[i] >> j) & 1:
                rows[i] ^= rows[j]
    return tuple(r >> d for r in rows)


def padd(*polys):
    out = {}
    for p in polys:
        for mon, a in p.items():
            out[mon] = mxor(out[mon], a) if mon in out else a
    return {m: a for m, a in out.items() if any(a)}


def pmul(p, q):
    out = {}
    for m, a in p.items():
        for n, b in q.items():
            mon = tuple(x+y for x, y in zip(m, n))
            ab = mm(a, b)
            out[mon] = mxor(out[mon], ab) if mon in out else ab
    return {m: a for m, a in out.items() if any(a)}


def factor(a, j, g, sign=1):
    mon = tuple(sign if k == j else 0 for k in range(g))
    return padd({(0,)*g: mxor(eye(len(a)), a)}, {mon: a})


def prepare(projections, order):
    g, d = len(projections), len(projections[0])
    unit = {(0,)*g: eye(d)}
    pref = [[unit for _ in range(g+1)]]
    inv = [[unit for _ in range(g+1)]]
    rhs = [[{} for _ in range(g+1)] for _ in range(g+1)]
    for i in range(g+1):
        rhs[i][i] = unit
    for j, half in order:
        a = projections[j] if half == 0 else mxor(eye(d), projections[j])
        ap = {(0,)*g: a}
        for x in (0, j+1):
            for y in (0, j+1):
                term = pmul(pmul(inv[-1][x], ap), pref[-1][y])
                rhs[x][y] = padd(rhs[x][y], term)
        p, q = list(pref[-1]), list(inv[-1])
        for x in (0, j+1):
            p[x] = pmul(factor(a, j, g), p[x])
            q[x] = pmul(q[x], factor(a, j, g, -1))
        pref.append(p)
        inv.append(q)
    expected = {(1,)*g: eye(d)}
    if pref[-1][0] != expected:
        return None
    for j in range(g):
        if pref[-1][j+1] != {tuple(int(k == j) for k in range(g)): eye(d)}:
            return None
    return pref, inv, rhs


def block_equations(left, right, rhs, d):
    """Build every Laurent coefficient equation for constant thin blocks."""
    equations = {}
    nvar = len(left)*d*d
    for mon, a in rhs.items():
        for i, row in enumerate(a):
            while row:
                low = row & -row
                equations[(mon, i, low.bit_length()-1)] = 1 << nvar
                row ^= low
    for t, (lp, rp) in enumerate(zip(left, right)):
        for m, a in lp.items():
            for n, b in rp.items():
                mon = tuple(x+y for x, y in zip(m, n))
                bt = mt(b)
                for i, row in enumerate(a):
                    for j, col in enumerate(bt):
                        mask = 0
                        bits = row
                        while bits:
                            low = bits & -bits
                            mask ^= col << ((low.bit_length()-1)*d)
                            bits ^= low
                        if mask:
                            key = (mon, i, j)
                            equations[key] = equations.get(key, 0) ^ (mask << (t*d*d))
    return equations, nvar


def reduce_equations(equations, nvar):
    pivots = {}
    mask = (1 << nvar)-1
    for row in equations.values():
        while row & mask:
            j = (row & mask).bit_length()-1
            if j not in pivots:
                pivots[j] = row
                break
            row ^= pivots[j]
        else:
            if row >> nvar:
                return None
    return pivots


def solve_block(left, right, rhs, d):
    """Solve every Laurent coefficient and replay the chosen solution."""
    equations, nvar = block_equations(left, right, rhs, d)
    pivots = reduce_equations(equations, nvar)
    if pivots is None:
        return None
    answer = 0
    for j in sorted(pivots):
        row = pivots[j]
        if ((row & answer).bit_count() ^ (row >> nvar)) & 1:
            answer |= 1 << j
    blocks = [tuple((answer >> (t*d*d+i*d)) & ((1 << d)-1)
                    for i in range(d)) for t in range(len(left))]
    # Replay the chosen solution without the equation builder.
    actual = {}
    zero = (0,)*len(next(iter(left[0])))
    for lp, z, rp in zip(left, blocks, right):
        actual = padd(actual, pmul(pmul(lp, {zero: z}), rp))
    if actual != rhs:
        raise ValueError('Chosen solution fails the complete polynomial endpoint')
    return {'blocks': blocks, 'free_variables': nvar-len(pivots)}


def solve(projections, order, gaps, prepared=None):
    data = prepared or prepare(projections, order)
    if data is None:
        return {'status': 'pure_stage_failure'}
    p, q, rhs = data
    g, d = len(projections), len(projections[0])
    chosen = {}
    pairs = ([(0, j) for j in range(1, g+1)] +
             [(j, 0) for j in range(1, g+1)] +
             [(j, j) for j in range(g+1)])
    for x, y in pairs:
        result = solve_block([q[k][x] for k in gaps],
                             [p[k][y] for k in gaps], rhs[x][y], d)
        if result is None:
            return {'status': 'inconsistent', 'block': [x, y]}
        chosen[f'{x},{y}'] = result
    thin = [[0]*(d*(g+1)) for _ in gaps]
    for key, result in chosen.items():
        x, y = map(int, key.split(','))
        for t, block in enumerate(result['blocks']):
            for i, row in enumerate(block):
                thin[t][x*d+i] ^= row << (y*d)
    return {'status': 'endpoint', 'thin_ranks': list(map(rank, thin)),
            'thin': thin, 'free_variables': sum(v['free_variables'] for v in chosen.values())}


def random_projection(d, rng, graph=False):
    q = d//2
    if graph:
        return tuple((1 << i) | (rng.randrange(1 << q) << q) if i < q else 0
                     for i in range(d))
    a = eye(d)
    for _ in range(8*d):
        i, j = rng.sample(range(d), 2)
        a = tuple(x ^ a[j] if k == i else x for k, x in enumerate(a))
    p = tuple(1 << i if i < q else 0 for i in range(d))
    return mm(mm(a, p), inverse(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--g', type=int, default=5)
    ap.add_argument('--d', type=int, default=4)
    ap.add_argument('--thin', type=int, default=2)
    ap.add_argument('--trials', type=int, default=10)
    ap.add_argument('--seed', type=int, default=5747)
    ap.add_argument('--graph', action='store_true')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    start = time.monotonic()
    order = [(j, 0) for j in range(args.g)] + [(j, 1) for j in reversed(range(args.g))]
    receipt = {'parameters': vars(args) | {'output': str(args.output)},
               'status': 'bounded_screen', 'systems': 0, 'successes': []}
    for trial in range(args.trials):
        projections = [random_projection(args.d, rng, args.graph) for _ in range(args.g)]
        if any(mm(a, a) != a or rank(a) != args.d//2 for a in projections):
            raise ValueError('Invalid projection')
        data = prepare(projections, order)
        if data is None:
            raise ValueError('Nested pure-stage endpoint failed')
        for gaps in itertools.combinations(range(2*args.g+1), args.thin):
            ans = solve(projections, order, gaps, data)
            receipt['systems'] += 1
            if ans['status'] == 'endpoint':
                receipt['successes'].append({'trial': trial, 'gaps': gaps,
                    'projections': projections, 'order': order, 'result': ans})
                print(json.dumps({'trial': trial, 'gaps': gaps, 'ranks': ans['thin_ranks']}), flush=True)
                args.output.write_text(json.dumps(receipt, indent=2)+'\n')
                break
        print(json.dumps({'trial': trial, 'systems': receipt['systems'],
                          'successes': len(receipt['successes'])}), flush=True)
    receipt['seconds'] = time.monotonic()-start
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')


if __name__ == '__main__':
    main()
