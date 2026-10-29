"""Test full stages with coefficients in an ordinary involution J.

The query S does not commute with J. The S3 mode uses all three rational
irreducible representations. A search result is not a native caller.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
import exact_replay as e


def construct(n):
    us, vs, _ = e.canonical(n)
    names, sectors = [], []

    def variable(name):
        index = len(names)
        names.append(name)
        return {(index,): F(1)}

    for sector in ('plus', 'minus'):
        ps = []
        for site, (u, v) in enumerate(zip(us, vs)):
            pivot = next(j for j, x in enumerate(v) if x)
            drop = next(j for j, x in enumerate(u) if x and j != pivot)
            h = [F(int(j == pivot), v[pivot]) for j in range(n)]
            columns = [u, h] + [[int(k == j) - v[j]*h[k] for k in range(n)]
                                for j in range(n) if j not in (pivot, drop)]
            basis = list(map(list, zip(*columns)))
            a = [variable(f'{sector}_a{site}_{j}') for j in range(n-2)]
            b = [variable(f'{sector}_b{site}_{j}') for j in range(n-2)]
            # The rank-two chart has no c variable when n=2.
            c = variable(f'{sector}_c{site}') if n > 2 else {}
            core = e.zero(n)
            core[0][0] = core[1][1] = e.const(1)
            core[0][1] = c
            for j in range(n-2):
                core[j+2][1], core[0][j+2] = a[j], b[j]
            ps.append(e.mm(e.mm(e.numeric(basis), core), e.numeric(e.inverse(basis))))
        sectors.append(ps)
    rs = [e.numeric([[u[i]*v[j] for j in range(n)] for i in range(n)])
          for u, v in zip(us, vs)]
    return names, sectors, rs


def diagonal(a, b):
    n = len(a)
    return [row + [{} for _ in range(n)] for row in a] + \
           [[{} for _ in range(n)] + row for row in b]


def tensor(a, b):
    return [[e.scale(b[i][j], a[r][s]) for s in range(len(a)) for j in range(len(b))]
            for r in range(len(a)) for i in range(len(b))]


def coefficient(ps, ts, degree):
    word = [e.eye(len(ps[0]))]
    for p, t in zip(ps, ts):
        old = [e.mm(t, a) for a in word]
        word[:] = [old[0]] + [e.ma(old[j], e.mm(p, old[j-1])) for j in range(1, len(old))]
        if len(old) <= degree:
            word.append(e.mm(p, old[-1]))
    return word[degree] if degree < len(word) else e.zero(len(ps[0]))


def equations(sectors, rs, degree, mode):
    plus, minus = sectors
    n = len(plus[0])
    standard_s = [[F(-1, 2), F(3, 2)], [F(1, 2), F(1, 2)]]
    if mode == 'swap':
        standard_s = [[F(0), F(1)], [F(1), F(0)]]
    s_minus_i = [[x - int(i == j) for j, x in enumerate(row)]
                 for i, row in enumerate(standard_s)]
    trivial = coefficient(plus, [e.eye(n) for _ in plus], degree)
    standard = coefficient([diagonal(a, b) for a, b in zip(plus, minus)],
                           [e.ma(e.eye(2*n), tensor(s_minus_i, r)) for r in rs], degree)
    matrices = [e.ma(standard, tensor(standard_s, trivial), -1)]
    if mode == 's3':
        signed = coefficient(minus, [e.ma(e.eye(n), r, -2) for r in rs], degree)
        matrices.append(e.ma(signed, trivial))
    return [entry for matrix in matrices for row in matrix for entry in row]


def search(n, mode, degree, seconds):
    start = time.monotonic()
    names, sectors, rs = construct(n)
    if any(equations(sectors, rs, 0, mode)):
        raise ValueError('The complete query control fails at z=1')
    steps = []
    status = 'UNRESOLVED'

    def expire(signum, frame):
        raise TimeoutError('The scalar coefficient time limit expired')

    signal.signal(signal.SIGALRM, expire)
    signal.alarm(seconds)
    try:
        values = e.linear_solution(equations(sectors, rs, 1, mode))
        if values is None:
            status = 'EXACT_FIRST_COEFFICIENT_CONTRADICTION'
        else:
            sectors = [[e.ms(p, values) for p in ps] for ps in sectors]
            system = [entry for ps in sectors for p in ps
                      for row in e.ma(e.mm(p, p), p, -1) for entry in row]
            for d in range(2, degree+1):
                system += equations(sectors, rs, d, mode)
                while True:
                    reduced = e.row_reduce(system)
                    affine = [p for p in reduced if all(len(w) <= 1 for w in p)]
                    step = {'degree': d, 'variables': len({v for p in system for w in p for v in w}),
                            'rank': len(reduced), 'affine': len(affine)}
                    steps.append(step)
                    print(json.dumps(step), flush=True)
                    if any(set(p) == {()} for p in affine):
                        status = 'EXACT_SCALAR_CONTRADICTION'
                        break
                    if not affine:
                        status = 'UNRESOLVED' if any(system) else 'COEFFICIENT_FAMILY'
                        break
                    values = e.linear_solution(affine)
                    sectors = [[e.ms(p, values) for p in ps] for ps in sectors]
                    system = [e.substitute(p, values) for p in system]
                if 'CONTRADICTION' in status:
                    break
    except TimeoutError:
        status = 'UNKNOWN_TIMEOUT'
    finally:
        signal.alarm(0)
    result = {'n': n, 'queries': n+1, 'mode': mode, 'initial_variables': len(names),
              'requested_degree': degree, 'status': status, 'steps': steps,
              'seconds': time.monotonic()-start, 'limit_seconds': seconds,
              'complete_query_control': 'PASS', 'new_kappa': False,
              'scope': 'Fixed canonical queries; scalar coefficients on two J eigenspaces; common stage endpoint'}
    if status in ('UNRESOLVED', 'COEFFICIENT_FAMILY'):
        result['chart_variables'] = sorted({v for ps in sectors for p in ps for row in p
                                            for entry in row for w in entry for v in w})
        result['chart'] = [[[[[[list(w), str(c)] for w, c in sorted(entry.items())]
                             for entry in row] for row in p] for p in ps] for ps in sectors]
        candidate = [[[[e.const(entry.get((), F(0))) for entry in row] for row in p]
                      for p in ps] for ps in sectors]
        result['constant_choice_projectors'] = not any(entry for ps in candidate for p in ps
            for row in e.ma(e.mm(p, p), p, -1) for entry in row)
        result['constant_choice_full_query'] = not any(entry for d in range(n+2)
            for entry in equations(candidate, rs, d, mode))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--mode', choices=('s3', 'swap'), default='s3')
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--seconds', type=int, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.n < 2 or args.n % 2:
        parser.error('The field count must be positive and even')
    result = search(args.n, args.mode, args.degree, args.seconds)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'chart'}))
