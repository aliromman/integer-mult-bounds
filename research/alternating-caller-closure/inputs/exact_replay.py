"""Rebuild the projector equations with rational sparse polynomials.

This file uses the Python standard library. It does not read the symbolic
discovery output. All matrix products use the stated chronological order.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time


def add(*polys):
    out = {}
    for poly in polys:
        for monomial, value in poly.items():
            out[monomial] = out.get(monomial, 0) + value
    return {m: v for m, v in out.items() if v}


def scale(poly, value):
    return {m: v * value for m, v in poly.items() if v * value}


def mul(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            key = tuple(sorted(a + b))
            out[key] = out.get(key, 0) + x * y
    return {m: v for m, v in out.items() if v}


def const(value):
    return {(): F(value)} if value else {}


def substitute(poly, values):
    out = {}
    for monomial, coefficient in poly.items():
        term = const(coefficient)
        for variable in monomial:
            term = mul(term, values.get(variable, {(variable,): F(1)}))
        out = add(out, term)
    return out


def zero(n):
    return [[{} for _ in range(n)] for _ in range(n)]


def eye(n):
    return [[const(int(i == j)) for j in range(n)] for i in range(n)]


def numeric(matrix):
    return [[const(x) for x in row] for row in matrix]


def mm(left, right):
    return [[add(*(mul(a, b) for a, b in zip(row, column)))
             for column in zip(*right)] for row in left]


def ma(left, right, sign=1):
    return [[add(a, scale(b, sign)) for a, b in zip(x, y)]
            for x, y in zip(left, right)]


def ms(matrix, values):
    return [[substitute(x, values) for x in row] for row in matrix]


def inverse(matrix):
    n = len(matrix)
    rows = [[F(x) for x in row] + [F(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next(i for i in range(column, n) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [x / divisor for x in rows[column]]
        for i in range(n):
            if i != column and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[column])]
    return [row[n:] for row in rows]


def row_reduce(equations, divisors=None):
    """Use polynomial monomials as independent columns for row operations."""
    columns = sorted({m for p in equations for m in p},
                     key=lambda m: (len(m), m), reverse=True)
    rows = [dict(p) for p in equations if p]
    rank = 0
    for column in columns:
        pivot = next((i for i in range(rank, len(rows)) if rows[i].get(column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        if divisors is not None:
            divisors.append(rows[rank][column])
        rows[rank] = scale(rows[rank], 1 / rows[rank][column])
        for i in range(len(rows)):
            if i != rank and rows[i].get(column):
                rows[i] = add(rows[i], scale(rows[rank], -rows[i][column]))
        rank += 1
        if rank == len(rows):
            break
    return rows[:rank]


def linear_solution(equations):
    reduced = row_reduce(equations)
    values = {}
    for row in reduced:
        if set(row) == {()}:
            return None
        pivot = max(m[0] for m in row if m)
        values[pivot] = {m: -v for m, v in row.items() if m != (pivot,)}
    return values


def canonical(n):
    us = [[int(i == j) for i in range(n)] for j in range(n)]
    us.append([(-1)**(j + 1) for j in range(n)])
    vs = [[0 if i == j else (1 if i > j else -1) for j in range(n)]
          for i in range(n + 1)]
    ts = [[[int(i == j) - 2*u[i]*v[j] for j in range(n)] for i in range(n)]
          for u, v in zip(us, vs)]
    return us, vs, ts


def chart(n, thin=()):
    us, vs, ts = canonical(n)
    names, ps, constraints = [], [], []

    def variable(name):
        index = len(names)
        names.append(name)
        return {(index,): F(1)}

    for site in range(1, n):
        if site in thin:
            ps.append(zero(n))
            continue
        u, v = us[site], vs[site]
        pivot = next(j for j in range(n) if v[j])
        a = [F(int(j == pivot), v[pivot]) for j in range(n)]
        columns = [u, a] + [[int(i == j) - v[j]*a[i] for i in range(n)]
                            for j in range(n) if j not in (site, pivot)]
        basis = [list(row) for row in zip(*columns)]
        av = [variable(f'a{site}_{j}') for j in range(n - 2)]
        bv = [variable(f'b{site}_{j}') for j in range(n - 2)]
        c = variable(f'c{site}')
        core = zero(n)
        core[0][0] = core[1][1] = const(1)
        core[0][1] = c
        for j in range(n - 2):
            core[j + 2][1], core[0][j + 2] = av[j], bv[j]
        ps.append(mm(mm(numeric(basis), core), numeric(inverse(basis))))
        constraints.append(add(c, *(mul(a, b) for a, b in zip(av, bv))))
    return names, ps, constraints, [numeric(t) for t in ts]


def coefficient_equations(ps, ts, degree):
    """Multiply the original words, truncated after the requested degree."""
    n = len(ts[0])
    plus, minus = [eye(n)], [ts[0]]
    for site, p in enumerate(ps, 1):
        for word, query in ((plus, eye(n)), (minus, ts[site])):
            old = [mm(query, x) for x in word]
            word[:] = [old[0]] + [ma(old[j], mm(p, old[j - 1]))
                                      for j in range(1, len(old))]
            if len(old) <= degree:
                word.append(mm(p, old[-1]))
    if degree >= len(plus):
        return [{} for _ in range(n*n)]
    return [entry for row in ma(mm(ts[-1], minus[degree]), plus[degree]) for entry in row]


def replay(n, thin=(), verbose=True, inverse_sites=()):
    started = time.monotonic()
    if n < 2 or n % 2:
        raise ValueError('Use a positive even field count')
    names, ps, constraints, ts = chart(n, thin)
    ps = [ma(eye(n),p,-1) if site in inverse_sites else p
          for site,p in enumerate(ps,1)]
    word = eye(n)
    for t in ts:
        word = mm(t, word)
    if ma(word, eye(n)) != zero(n):
        raise ValueError('The canonical query word does not equal minus identity')
    linear = coefficient_equations(ps, ts, 1)
    values = linear_solution(linear)
    steps = []
    if values is None:
        status = 'FIRST_COEFFICIENT_CONTRADICTION'
    else:
        ps = [ms(p, values) for p in ps]
        constraints = [substitute(p, values) for p in constraints]
        if any(coefficient_equations(ps, ts, 1)):
            raise ValueError('The first-coefficient solution does not replay')
        equations = constraints + coefficient_equations(ps, ts, 2)
        while True:
            reduced = row_reduce(equations)
            affine = [p for p in reduced if all(len(m) <= 1 for m in p)]
            variables = {v for p in equations for m in p for v in m}
            row = {'variables': len(variables), 'equations': sum(bool(p) for p in equations),
                   'linear_consequences': len(affine)}
            steps.append(row)
            if verbose:
                print(json.dumps(row), flush=True)
            if any(set(p) == {()} for p in affine):
                status = 'EXACT_SECOND_COEFFICIENT_CONTRADICTION'
                break
            if not affine:
                status = 'UNRESOLVED' if any(equations) else 'SECOND_COEFFICIENT_CONSISTENT'
                break
            values = linear_solution(affine)
            equations = [substitute(p, values) for p in equations]
    return {'status': status, 'n': n, 'thin_sites': list(thin),
            'initial_variables': len(names), 'steps': steps,
            'seconds': time.monotonic() - started, 'new_kappa': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--thin', type=int, nargs='*', default=[])
    parser.add_argument('--inverse', type=int, nargs='*', default=[])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = replay(args.n, tuple(args.thin), inverse_sites=tuple(args.inverse))
    if args.inverse:
        result['inverse_sites'] = args.inverse
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
