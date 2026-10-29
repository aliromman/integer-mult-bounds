"""Test the shortest-query equations over a free associative algebra.

Coefficients are rational. Variable products keep their order. A constant
contradiction applies to every matrix block width. A surviving system is
not a representation or a complete caller.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
import exact_replay as e


def mul(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            out[a+b] = out.get(a+b, 0)+x*y
    return {m: v for m, v in out.items() if v}


def mm(left, right):
    return [[e.add(*(mul(a, b) for a, b in zip(row, col)))
             for col in zip(*right)] for row in left]


def substitute(poly, values):
    terms = []
    for word, coefficient in poly.items():
        term = e.const(coefficient)
        for variable in word:
            term = mul(term, values.get(variable, {(variable,): F(1)}))
        terms.append(term)
    return e.add(*terms)


def ms(matrix, values):
    return [[substitute(x, values) for x in row] for row in matrix]


def coefficient(ps, ts, degree):
    n = len(ts[0])
    plus, minus = [e.eye(n)], [ts[0]]
    for site, p in enumerate(ps, 1):
        for word, query in ((plus, e.eye(n)), (minus, ts[site])):
            old = [mm(query, x) for x in word]
            word[:] = [old[0]] + [e.ma(old[j], mm(p, old[j-1]))
                                  for j in range(1, len(old))]
            if len(old) <= degree:
                word.append(mm(p, old[-1]))
    return [x for row in e.ma(mm(ts[-1], minus[degree]), plus[degree]) for x in row]


def check(args):
    started = time.monotonic()
    names, ps, _, ts = e.chart(args.n, tuple(args.thin))
    linear = coefficient(ps, ts, 1)
    values = e.linear_solution(linear)
    if values is None:
        result = {'n': args.n, 'thin_sites': args.thin,
                  'status': 'EXACT_FREE_ALGEBRA_LINEAR_CONTRADICTION',
                  'steps': [], 'seconds': time.monotonic()-started, 'new_kappa': False}
        args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result), flush=True)
        return
    ps = [ms(p, values) for p in ps]
    system = [x for p in ps for row in e.ma(mm(p, p), p, -1) for x in row]
    system += coefficient(ps, ts, 2)
    steps = []
    for extension in range(args.extensions+1):
        while True:
            reduced = e.row_reduce(system)
            affine = [p for p in reduced if all(len(m) <= 1 for m in p)]
            row = {'variables': len({v for p in system for m in p for v in m}),
                   'equations': len([p for p in system if p]),
                   'rank': len(reduced), 'affine_consequences': len(affine),
                   'extension': extension}
            steps.append(row)
            print(json.dumps(row), flush=True)
            if any(set(p) == {()} for p in affine):
                status = 'EXACT_FREE_ALGEBRA_CONTRADICTION'
                break
            if not affine:
                status = ('UNRESOLVED_FREE_ALGEBRA_SYSTEM' if any(system)
                          else 'EXACT_FREE_ALGEBRA_SECOND_COEFFICIENT_CONSISTENT')
                system = reduced
                break
            values = e.linear_solution(affine)
            system = [substitute(p, values) for p in system]
            ps = [ms(p, values) for p in ps]
        if status.startswith('EXACT'):
            break
        if extension < args.extensions:
            variables = sorted({v for p in system for m in p for v in m})
            system += [new for p in list(system) for v in variables
                       for new in (mul({(v,): F(1)}, p), mul(p, {(v,): F(1)}))]
    if args.chart and status.startswith('UNRESOLVED'):
        variables = sorted({v for p in ps for row in p for x in row for m in x for v in m})
        def serial(poly):
            return [{'word': [names[i] for i in word], 'coefficient': str(value)}
                    for word, value in sorted(poly.items())]
        chart = {'n': args.n, 'free': [names[i] for i in variables],
                 'projectors': [[[serial(x) for x in row] for row in p] for p in ps],
                 'queries': [[[int(x.get((), 0)) for x in row] for row in t] for t in ts],
                 'system': [serial(p) for p in system]}
        args.chart.write_text(json.dumps(chart, indent=2)+'\n')
    result = {'n': args.n, 'thin_sites': args.thin, 'status': status, 'steps': steps,
              'seconds': time.monotonic()-started, 'new_kappa': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, choices=(4, 6, 8), default=4)
    parser.add_argument('--extensions', type=int, default=0)
    parser.add_argument('--thin', type=int, nargs='*', default=[])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--chart', type=Path)
    check(parser.parse_args())
