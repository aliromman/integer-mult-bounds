"""Test initialized work with no condition on discarded outputs.

The query tuple is fixed. Projector charts are complete for that tuple.
Only the retained block must have the specified degree-two/three endpoint.
An unresolved system is not a caller or a cost certificate.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
import exact_replay as e


def construct(r, work, thin=()):
    n = r + work
    oldus, oldvs, _ = e.canonical(r)
    us = [u + [0]*work for u in oldus]
    vs = [v + [0]*work for v in oldvs]
    names, ps, qs, ts = [], [], [], []

    def variable(name):
        index = len(names)
        names.append(name)
        return {(index,): F(1)}

    for i, (u, v) in enumerate(zip(us, vs)):
        ts.append(e.numeric([[int(a == b)-2*u[a]*v[b] for b in range(n)]
                             for a in range(n)]))
        if i in thin:
            ps.append(e.zero(n))
            continue
        pivot = next(j for j in range(n) if v[j])
        drop = next(j for j in range(n) if u[j] and j != pivot)
        h = [F(int(j == pivot), v[pivot]) for j in range(n)]
        columns = [u, h] + [[int(k == j)-v[j]*h[k] for k in range(n)]
                            for j in range(n) if j not in (pivot, drop)]
        basis = list(map(list, zip(*columns)))
        av = [variable(f'a{i}_{j}') for j in range(n-2)]
        bv = [variable(f'b{i}_{j}') for j in range(n-2)]
        c = variable(f'c{i}')
        core = e.zero(n)
        core[0][0] = core[1][1] = e.const(1)
        core[0][1] = c
        for j in range(n-2):
            core[j+2][1], core[0][j+2] = av[j], bv[j]
        ps.append(e.mm(e.mm(e.numeric(basis), core), e.numeric(e.inverse(basis))))
        qs.append(e.add(c, *(e.mul(a, b) for a, b in zip(av, bv))))
    return names, ps, qs, ts


def coefficient(ps, ts, degree):
    n = len(ps[0])
    words = [[e.eye(n)], [e.eye(n)]]
    for p, t in zip(ps, ts):
        for channel, query in enumerate([e.eye(n), t]):
            old = [e.mm(query, x) for x in words[channel]]
            word = [old[0]] + [e.ma(old[j], e.mm(p, old[j-1]))
                                for j in range(1, len(old))]
            if len(old) <= degree:
                word.append(e.mm(p, old[-1]))
            words[channel] = word
    return words[0][degree], words[1][degree]


def retained(matrix, r):
    return [row[:r] for row in matrix[:r]]


def expression(poly, names):
    return '+'.join('('+str(value)+')' + ''.join('*'+names[i] for i in monomial)
                    for monomial, value in sorted(poly.items())) or '0'


def check(args):
    started = time.monotonic()
    r, work = args.r, args.work
    thin = tuple(range(2, r, 2)) if args.alternating else ()
    base = 1 if args.alternating else 2
    higher = 2*(r+1-len(thin))-base*r
    names, ps, qs, ts = construct(r, work, thin)
    plus, minus = coefficient(ps, ts, 1)
    linear = [x for row in retained(e.ma(plus, minus), r) for x in row]
    values = e.linear_solution(linear)
    if values is None:
        raise ValueError('The first-coefficient chart is inconsistent')
    ps = [e.ms(p, values) for p in ps]
    qs = [e.substitute(p, values) for p in qs]
    h = retained(e.ms(plus, values), r)
    k = e.ma(h, [[e.scale(x, base) for x in row] for row in e.eye(r)], -1)
    if e.add(*(k[i][i] for i in range(r))) != e.const(higher):
        raise ValueError('The target rank control fails')
    plus2, minus2 = coefficient(ps, ts, 2)
    expected = e.ma([[e.scale(x, base*(base-1)//2) for x in row] for row in e.eye(r)],
                    [[e.scale(x, base) for x in row] for row in k])
    constraints = [e.ma(e.mm(k, k), k, -1),
                   e.ma(retained(plus2, r), expected, -1),
                   e.ma(retained(minus2, r), expected)]
    system = qs + [x for matrix in constraints for row in matrix for x in row]
    steps = []
    while True:
        reduced = e.row_reduce(system)
        affine = [p for p in reduced if all(len(m) <= 1 for m in p)]
        row = {'variables': len({v for p in system for m in p for v in m}),
               'equations': len([p for p in system if p]),
               'affine_consequences': len(affine)}
        steps.append(row)
        print(json.dumps(row), flush=True)
        if any(set(p) == {()} for p in affine):
            status = 'EXACT_INITIALIZED_SECOND_COEFFICIENT_CONTRADICTION'
            break
        if not affine:
            status = 'UNRESOLVED_INITIALIZED_SYSTEM'
            break
        values = e.linear_solution(affine)
        system = [e.substitute(p, values) for p in system]
        ps = [e.ms(p, values) for p in ps]
        qs = [e.substitute(p, values) for p in qs]
        k = e.ms(k, values)
    if args.chart and status.startswith('UNRESOLVED'):
        used = sorted({v for p in ps for row in p for x in row for m in x for v in m})
        serial = lambda matrix: [[expression(x, names) for x in row] for row in matrix]
        chart = {'r': r, 'work': work, 'target_base': base,
                 'higher_target_count': higher, 'free': [names[i] for i in used],
                 'projectors': [serial(p) for p in ps],
                 'quadratics': [expression(p, names) for p in qs],
                 'target_projector': serial(k),
                 'queries': [[[int(x.get((), 0)) for x in row] for row in t] for t in ts]}
        args.chart.write_text(json.dumps(chart, indent=2)+'\n')
    result = {'r': r, 'initialized_work': work, 'full_calls': r+1-len(thin),
              'thin_sites': list(thin),
              'target_stage_powers': [base]*(r-higher)+[base+1]*higher, 'status': status,
              'steps': steps, 'seconds': time.monotonic()-started, 'new_kappa': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, choices=(4, 6, 8, 10, 12, 16), default=6)
    parser.add_argument('--work', type=int, choices=(1, 2, 3, 4, 6, 8), default=1)
    parser.add_argument('--alternating', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--chart', type=Path)
    check(parser.parse_args())
