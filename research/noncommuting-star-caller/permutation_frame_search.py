"""Search full-call frames that depend on an ordinary block permutation.

The coefficient algebra is F2[(C2)^3 semidirect C3]. The three Hadamard
generators commute and have order two. R cycles them. U=R H_0 advances
one block. The central query is 1+e, with e^2=0.

The fixed query vectors come from Cinna Davis's simplex frame. The
complementary vectors have coefficients in F2[R]. A positive result
needs a characteristic-zero lift and a complete physical caller.
"""
import argparse
import json
from pathlib import Path
import time

import z3
import screen as s


ELEMENTS = [(h, r) for h in range(8) for r in range(3)]
INDEX = {g: i for i, g in enumerate(ELEMENTS)}


def rotate(h, r):
    return ((h << r) | (h >> (3-r))) & 7


def product(a, b):
    h, r = a
    k, t = b
    return h ^ rotate(k, r), (r+t) % 3


TABLE = [[INDEX[product(a, b)] for b in ELEMENTS] for a in ELEMENTS]
ONE = 1 << INDEX[0, 0]
R = INDEX[0, 1]
U = INDEX[2, 1]


def integer_left(g, x):
    value = 0
    for j in range(24):
        if x >> j & 1:
            value ^= 1 << TABLE[g][j]
    return value


def integer_mul(a, b):
    value = 0
    for i in range(24):
        if a >> i & 1:
            value ^= integer_left(i, b)
    return value


def ring_element(poly):
    return sum(((poly >> r) & 1) << INDEX[0, r] for r in range(3))


def power(g, e):
    value = 0
    for _ in range(e):
        value = TABLE[g][value]
    return value


def xor_all(xs, zero=0):
    for x in xs:
        zero ^= x
    return zero


def replay(frame, aa, bb, exponents):
    """Execute each full call in the literal finite group algebra."""
    n = len(aa[0])
    mat = [[ONE if i == j else 0 for j in range(n)] for i in range(n)]
    qry = [[0]*n for _ in range(n)]
    for (u, v), a, b in zip(frame, aa, bb):
        a, b = list(map(ring_element, a)), list(map(ring_element, b))
        if (xor_all(a[j] for j in range(n) if v >> j & 1) != ONE
                or xor_all(b[j] for j in range(n) if u >> j & 1) != ONE
                or xor_all(integer_mul(b[j], a[j]) for j in range(n))):
            raise ValueError('The two frame vectors are not dual')
        def apply(x):
            t0 = [xor_all(integer_mul(b[j], x[j][k]) for j in range(n)) for k in range(n)]
            t1 = [xor_all(x[j][k] for j in range(n) if v >> j & 1) for k in range(n)]
            d0 = [integer_left(U, y) ^ y for y in t0]
            d1 = [integer_left(U, y) ^ y for y in t1]
            return [[x[i][k] ^ (d0[k] if u >> i & 1 else 0)
                     ^ integer_mul(a[i], d1[k]) for k in range(n)] for i in range(n)]
        new_qry = apply(qry)
        source = [xor_all(mat[j][k] for j in range(n) if v >> j & 1) for k in range(n)]
        for i in range(n):
            if u >> i & 1:
                new_qry[i] = [x ^ integer_left(U, y) for x, y in zip(new_qry[i], source)]
        mat, qry = apply(mat), new_qry
    if mat != qry:
        raise ValueError('The full query endpoint fails')
    coefficients = [tuple(sum(((mat[i][j] >> g) & 1) << j for j in range(n))
                          for i in range(n)) for g in range(24)]
    target_groups = [power(U, e) for e in exponents]
    total = (0,)*n
    for g, coefficient in enumerate(coefficients):
        if s.rank(coefficient) != target_groups.count(g) or s.mm(coefficient, coefficient) != coefficient:
            raise ValueError('A target coefficient has the wrong rank or square')
        total = s.mxor(total, coefficient)
        for h, other in enumerate(coefficients):
            if g != h and any(s.mm(coefficient, other)):
                raise ValueError('The target images are not independent')
    if total != s.eye(n):
        raise ValueError('The target coefficients do not sum to the identity')
    return coefficients


def search(args):
    n = args.n
    frame = [(1 << i, ((1 << n)-1) ^ (1 << i)) for i in range(n)]
    frame.append(((1 << n)-1, (1 << n)-1))
    exponents = [int(x) for x in args.exponents.split(',')]
    if len(exponents) != n or any(x < 0 or x >= 6 for x in exponents):
        raise ValueError('Use one exponent in [0,6) per target field')
    solver = z3.Solver()
    solver.set(timeout=args.timeout, random_seed=args.seed)
    zero = z3.BitVecVal(0, 24)
    def bit(x, k):
        return z3.Extract(k, k, x) == 1
    def left(g, x):
        inverse = [TABLE[g].index(j) for j in range(24)]
        return z3.Concat(*(z3.Extract(inverse[j], inverse[j], x) for j in range(23, -1, -1)))
    def poly_left(p, x):
        return xor_all((z3.If(bit(p, c), left(INDEX[0, c], x), zero) for c in range(3)), zero)
    def poly_mul(p, q):
        out = z3.BitVecVal(0, 3)
        for c in range(3):
            out ^= z3.If(bit(p, c), z3.RotateLeft(q, c), z3.BitVecVal(0, 3))
        return out
    mat = [[z3.BitVecVal(ONE if i == j else 0, 24) for j in range(n)] for i in range(n)]
    qry = [[zero]*n for _ in range(n)]
    variables = []
    for stage, (u, v) in enumerate(frame):
        a = [z3.BitVec(f'a_{stage}_{i}', 3) for i in range(n)]
        b = [z3.BitVec(f'b_{stage}_{i}', 3) for i in range(n)]
        variables.append((a, b))
        solver.add(xor_all(a[j] for j in range(n) if v >> j & 1) == 1)
        solver.add(xor_all(b[j] for j in range(n) if u >> j & 1) == 1)
        solver.add(xor_all(poly_mul(b[j], a[j]) for j in range(n)) == 0)
        def apply(x, label):
            t0 = [xor_all((poly_left(b[j], x[j][k]) for j in range(n)), zero) for k in range(n)]
            t1 = [xor_all((x[j][k] for j in range(n) if v >> j & 1), zero) for k in range(n)]
            d0 = [left(U, y) ^ y for y in t0]
            d1 = [left(U, y) ^ y for y in t1]
            result = [[z3.BitVec(f'{label}_{stage}_{i}_{k}', 24) for k in range(n)] for i in range(n)]
            for i in range(n):
                for k in range(n):
                    value = x[i][k] ^ (d0[k] if u >> i & 1 else zero) ^ poly_left(a[i], d1[k])
                    if label == 'B' and u >> i & 1:
                        query = xor_all((mat[j][k] for j in range(n) if v >> j & 1), zero)
                        value ^= left(U, query)
                    solver.add(result[i][k] == value)
            return result
        new_qry = apply(qry, 'B')
        mat, qry = apply(mat, 'A'), new_qry
    basis = [z3.BitVec(f'basis_{i}', n) for i in range(n)]
    inverse = [z3.BitVec(f'inverse_{i}', n) for i in range(n)]
    def choose(selector, rows):
        return xor_all((z3.If(bit(selector, j), rows[j], z3.BitVecVal(0, n)) for j in range(n)), z3.BitVecVal(0, n))
    groups = [power(U, e) for e in exponents]
    for i in range(n):
        solver.add(choose(basis[i], inverse) == 1 << i)
        for g in range(24):
            selector = sum((h == g) << j for j, h in enumerate(groups))
            row = choose(basis[i] & selector, inverse)
            for j in range(n):
                solver.add(bit(mat[i][j], g) == bit(row, j))
        solver.add([mat[i][j] == qry[i][j] for j in range(n)])
    start = time.monotonic()
    status = str(solver.check())
    result = {'status': status, 'solve_seconds': time.monotonic()-start,
              'scope': 'Finite binary group-algebra endpoint only; no characteristic-zero or native caller claim',
              'parameters': vars(args) | {'output': str(args.output)}, 'query_frame': frame,
              'group_elements': ELEMENTS, 'advancer': ELEMENTS[U]}
    if status == 'unknown':
        result['reason'] = solver.reason_unknown()
    if status == 'sat':
        model = solver.model()
        aa = [[model.eval(x).as_long() for x in a] for a, _ in variables]
        bb = [[model.eval(x).as_long() for x in b] for _, b in variables]
        coefficients = replay(frame, aa, bb, exponents)
        result.update(a=aa, b=bb, coefficients=coefficients, coefficient_ranks=list(map(s.rank, coefficients)))
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('coefficients', 'group_elements')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--exponents', default='3,3,2,2')
    parser.add_argument('--timeout', type=int, default=60000)
    parser.add_argument('--seed', type=int, default=5786)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.n < 2 or args.n % 2:
        parser.error('Use a positive even field count')
    search(args)
