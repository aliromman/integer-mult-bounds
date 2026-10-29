"""Test a four-call extension of two disjoint full calls.

The six-field query endpoint is constructed first. The stage equations
then use the complete compatible rank-two projector charts.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
import exact_replay as e


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def work_frame(b, t):
    ident = [[F(i == j) for j in range(2)] for i in range(2)]
    first = [[F(1), F(1)], [F(0), F(1)]]
    second = [[F(1), F(0)], [b, F(1)]]
    denominator = 1-b*t*t-b*t
    if not denominator:
        raise ValueError('The chosen rational work chart has a zero denominator')
    coefficient = -(4+b)/denominator
    third = [[1+coefficient*t, -coefficient*t*t], [coefficient, 1-coefficient*t]]
    fourth = [[-x for x in row] for row in e.inverse(mm(third, mm(second, first))) ]
    pairs = []
    product = ident
    for query in (first, second, third, fourth):
        update = [[query[i][j]-ident[i][j] for j in range(2)] for i in range(2)]
        pivot = next((i, j) for i in range(2) for j in range(2) if update[i][j])
        i, j = pivot
        u = [row[j] for row in update]
        v = [-x/(2*update[i][j]) for x in update[i]]
        if sum(x*y for x, y in zip(u, v)) or \
           any(query[a][c] != ident[a][c]-2*u[a]*v[c] for a in range(2) for c in range(2)):
            raise ValueError('The work factor is not a rank-one transvection')
        pairs.append((u, v))
        product = mm(query, product)
    if product != [[-x for x in row] for row in ident]:
        raise ValueError('The four-query work product fails')
    return pairs


def frame(b, t):
    pairs = work_frame(b, t)
    work = [list(row) for row in zip(*(u for u, v in pairs))]
    covectors = [v for u, v in pairs]
    gram = mm(covectors, work)
    d = [[F(i == j)+2*gram[i][j]*int(i > j) for j in range(4)] for i in range(4)]
    z = [[d[i][j]-gram[i][j] for j in range(4)] for i in range(4)]
    inverse_d = e.inverse(d)
    right = mm(inverse_d, covectors)
    for columns in combinations(range(4), 2):
        a = [[row[j] for j in columns] for row in z]
        left = mm(inverse_d, a)
        try:
            inverse = e.inverse([r+l for r, l in zip(right, left)])
        except StopIteration:
            continue
        bmap = inverse[2:]
        break
    else:
        raise ValueError('No rank-two complement basis was found')
    images = bmap+work+[[F(0)]*4 for _ in range(2)]
    rows = [[F(0), F(0)]+v+aa for v, aa in zip(covectors, a)]
    full_pairs = [([row[j] for row in images], rows[j]) for j in range(4)]
    ident = [[F(i == j) for j in range(6)] for i in range(6)]
    target = [row[:] for row in ident]
    target[0][4] = target[1][5] = F(-2)
    target[2][2] = target[3][3] = F(-1)
    product = ident
    queries = []
    for u, v in full_pairs:
        query = [[ident[i][j]-2*u[i]*v[j] for j in range(6)] for i in range(6)]
        queries.append(query)
        product = mm(query, product)
    if product != target:
        raise ValueError('The complete old-transfer and new-query endpoint fails')
    return full_pairs, queries, target


def construct(b, t):
    pairs, queries, target = frame(b, t)
    names, projectors = [], []

    def variable(name):
        index = len(names)
        names.append(name)
        return {(index,): F(1)}

    for site, (u, v) in enumerate(pairs):
        pivot = next(j for j, x in enumerate(v) if x)
        drop = next(j for j, x in enumerate(u) if x and j != pivot)
        h = [F(int(j == pivot), v[pivot]) for j in range(6)]
        columns = [u, h]+[[int(k == j)-v[j]*h[k] for k in range(6)]
                          for j in range(6) if j not in (pivot, drop)]
        basis = list(map(list, zip(*columns)))
        a = [variable(f'a{site}_{j}') for j in range(4)]
        bb = [variable(f'b{site}_{j}') for j in range(4)]
        core = e.zero(6)
        core[0][0] = core[1][1] = e.const(1)
        core[0][1] = variable(f'c{site}')
        for j in range(4):
            core[j+2][1], core[0][j+2] = a[j], bb[j]
        projectors.append(e.mm(e.mm(e.numeric(basis), core), e.numeric(e.inverse(basis))))
    return names, projectors, [e.numeric(q) for q in queries], e.numeric(target)


def coefficient(ps, ts, degree):
    word = [e.eye(6)]
    for p, t in zip(ps, ts):
        old = [e.mm(t, a) for a in word]
        word[:] = [old[0]]+[e.ma(old[j], e.mm(p, old[j-1])) for j in range(1, len(old))]
        if len(old) <= degree:
            word.append(e.mm(p, old[-1]))
    return word[degree] if degree < len(word) else e.zero(6)


def equations(ps, ts, target, degree):
    if degree == 0:
        stage = e.eye(6)
    elif degree == 1:
        stage = e.numeric([[int(i == j)*(2 if i in (2, 3) else 1) for j in range(6)] for i in range(6)])
    elif degree == 2:
        stage = e.numeric([[int(i == j and i in (2, 3)) for j in range(6)] for i in range(6)])
    else:
        stage = e.zero(6)
    plus = e.ma(coefficient(ps, [e.eye(6) for _ in ps], degree), stage, -1)
    minus = e.ma(coefficient(ps, ts, degree), e.mm(stage, target), -1)
    return [entry for matrix in (plus, minus) for row in matrix for entry in row]


def search(b, t, seconds):
    start = time.monotonic()
    names, ps, ts, target = construct(b, t)
    if any(equations(ps, ts, target, 0)):
        raise ValueError('The degree-zero endpoint control fails')
    steps = []
    status = 'UNRESOLVED'

    def expire(signum, frame):
        raise TimeoutError('The exact module time limit expired')

    signal.signal(signal.SIGALRM, expire)
    signal.alarm(seconds)
    try:
        values = e.linear_solution(equations(ps, ts, target, 1))
        if values is None:
            status = 'EXACT_FIRST_COEFFICIENT_CONTRADICTION'
        else:
            ps = [e.ms(p, values) for p in ps]
            system = [entry for p in ps for row in e.ma(e.mm(p, p), p, -1) for entry in row]
            for degree in range(2, 5):
                system += equations(ps, ts, target, degree)
                while True:
                    reduced = e.row_reduce(system)
                    affine = [p for p in reduced if all(len(w) <= 1 for w in p)]
                    step = {'degree': degree, 'variables': len({v for p in system for w in p for v in w}),
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
                    ps = [e.ms(p, values) for p in ps]
                    system = [e.substitute(p, values) for p in system]
                if 'CONTRADICTION' in status:
                    break
    except TimeoutError:
        status = 'UNKNOWN_TIMEOUT'
    finally:
        signal.alarm(0)
    result = {'b': str(b), 't': str(t), 'fields': 6, 'full_calls': 4,
              'initial_variables': len(names), 'status': status, 'steps': steps,
              'query_endpoint': 'PASS_EXACT_COMPLETE', 'seconds': time.monotonic()-start,
              'limit_seconds': seconds, 'new_kappa': False,
              'scope': 'Fixed rational query frame; all compatible scalar full-stage projectors'}
    if status in ('UNRESOLVED', 'COEFFICIENT_FAMILY'):
        result['chart'] = [[[[list(w), str(c)] for w, c in sorted(entry.items())]
                            for row in p for entry in row] for p in ps]
        candidate = [[ [e.const(entry.get((), F(0))) for entry in row] for row in p] for p in ps]
        result['constant_choice_complete'] = not any(entry for degree in range(5)
            for entry in equations(candidate, ts, target, degree)) and not any(entry for p in candidate
            for row in e.ma(e.mm(p, p), p, -1) for entry in row)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b', type=F, default=F(1))
    parser.add_argument('--t', type=F, default=F(1))
    parser.add_argument('--seconds', type=int, default=90)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = search(args.b, args.t, args.seconds)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'chart'}))
