"""Rebuild the query equations over a free associative algebra.

Words retain their order. Rational row operations and complete affine
substitutions prove contradictions for every matrix block width.
Prepared by Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
from fractions import Fraction as F
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


def first_chart(n, thin):
    names, ps, _, ts = e.chart(n, thin)
    values = e.linear_solution(coefficient(ps, ts, 1))
    return names, None if values is None else [ms(p, values) for p in ps], ts


def replay(n, thin=()):
    names, ps, ts = first_chart(n, thin)
    steps = []
    if ps is None:
        status = 'FIRST_COEFFICIENT_CONTRADICTION'
    else:
        system = [x for p in ps for row in e.ma(mm(p, p), p, -1) for x in row]
        system += coefficient(ps, ts, 2)
        while True:
            reduced = e.row_reduce(system)
            affine = [p for p in reduced if all(len(m) <= 1 for m in p)]
            steps.append({'variables': len({v for p in system for m in p for v in m}),
                          'equations': sum(bool(p) for p in system),
                          'rank': len(reduced), 'linear_consequences': len(affine)})
            if any(set(p) == {()} for p in affine):
                status = 'EXACT_FREE_ALGEBRA_CONTRADICTION'
                break
            if not affine:
                status = 'UNRESOLVED' if any(system) else 'SECOND_COEFFICIENT_CONSISTENT'
                break
            values = e.linear_solution(affine)
            system = [substitute(p, values) for p in system]
    return {'n': n, 'thin_sites': list(thin), 'initial_matrix_variables': len(names),
            'status': status, 'steps': steps}


def numeric_product(left, right):
    return [[sum(a*b for a, b in zip(row, col)) for col in zip(*right)] for row in left]


def numeric_sum(left, right, sign=1):
    return [[a+sign*b for a, b in zip(x, y)] for x, y in zip(left, right)]


def value(poly, variables, width):
    out = [[F(0) for _ in range(width)] for _ in range(width)]
    for word, scalar in poly.items():
        term = [[scalar*int(i == j) for j in range(width)] for i in range(width)]
        for variable in word:
            term = numeric_product(term, variables[variable])
        out = numeric_sum(out, term)
    return out


def lift(matrix, variables, width):
    blocks = [[value(p, variables, width) for p in row] for row in matrix]
    return [[x for block in row for x in block[i]] for row in blocks for i in range(width)]


def matrix_controls():
    """Compare free words with literal finite matrix multiplication."""
    checks = []
    for n, width in [(4, 2), (6, 3)]:
        names, ps, ts = first_chart(n, ())
        variables = {}
        for v in range(len(names)):
            variables[v] = [[F((v+1)*(i+1)-j) if i == j else F((v+i+2*j) % 5-2)
                             for j in range(width)] for i in range(width)]
        size = n*width
        ident = [[F(i == j) for j in range(size)] for i in range(size)]
        pnums = [lift(p, variables, width) for p in ps]
        tnums = [lift(t, variables, width) for t in ts]
        plus, minus = [ident], [tnums[0]]
        for site, p in enumerate(pnums, 1):
            for word, query in ((plus, ident), (minus, tnums[site])):
                old = [numeric_product(query, x) for x in word]
                word[:] = [old[0]] + [numeric_sum(old[j], numeric_product(p, old[j-1]))
                                      for j in range(1, len(old))]
                if len(old) <= 2:
                    word.append(numeric_product(p, old[-1]))
        direct = numeric_sum(numeric_product(tnums[-1], minus[2]), plus[2])
        symbolic = coefficient(ps, ts, 2)
        symbolic_matrix = [symbolic[i*n:(i+1)*n] for i in range(n)]
        if lift(symbolic_matrix, variables, width) != direct:
            raise ValueError('The free-word coefficient disagrees with literal matrix products')
        constraints = symbolic + [x for p in ps for row in e.ma(mm(p, p), p, -1) for x in row]
        changed = False
        for p in constraints:
            commuted = {}
            for word, scalar in p.items():
                key = tuple(sorted(word))
                commuted[key] = commuted.get(key, 0)+scalar
            if value(p, variables, width) != value(commuted, variables, width):
                changed = True
                break
        if not changed:
            raise ValueError('The noncommutation control does not reject commutative multiplication')
        checks.append({'n': n, 'width': width, 'literal_second_coefficient': 'PASS',
                       'commuted_word_mutation_rejected': True})
    return checks


def corner_controls():
    checks = []
    for n in (2, 4, 8, 16):
        scalar = [[F(0 if i == j else (1 if i > j else -1)) for j in range(n)]
                  for i in range(n)]
        inv = e.inverse(scalar)
        for width in (2, 3):
            ident = [[F(i == j) for j in range(width)] for i in range(width)]
            x = [[F(i+2*j+1) for j in range(width)] for i in range(width)]
            y = [[F(2*i-j+(i == j)) for j in range(width)] for i in range(width)]
            if numeric_product(x, y) == numeric_product(y, x):
                raise ValueError('The corner control needs noncommuting blocks')
            row = [x[i]+ident[i]*(n-1) for i in range(width)]
            col = [[-v for v in line] for block in [y]+[ident]*(n-1) for line in block]
            inverse_block = [[inv[i//width][j//width]*int(i % width == j % width)
                              for j in range(n*width)] for i in range(n*width)]
            product = numeric_product(numeric_product(row, inverse_block), col)
            schur = [[-v for v in line] for line in product]
            if schur != numeric_sum(x, y, -1):
                raise ValueError('The block corner Schur identity fails')
            checks.append({'n': n, 'width': width, 'schur_complement': 'X-Y'})
    return checks
