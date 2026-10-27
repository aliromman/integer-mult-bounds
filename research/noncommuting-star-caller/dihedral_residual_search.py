"""Binary full-frame search over an ordinary swap of two axis blocks.

U=R H0 has order four. U^2 is central. The target is U^2 S K(R), with
K invertible, plus an ordinary C(R) that retained input copies can remove.
This is a finite algebraic screen, not a characteristic-zero caller.
"""
import argparse
import json
from pathlib import Path
import time

import z3
import screen as binary


ELEMENTS = [(h, r) for h in range(4) for r in range(2)]
INDEX = {g: i for i, g in enumerate(ELEMENTS)}
def product(a, b):
    h, r = a
    k, t = b
    return h ^ (k if not r else ((k << 1) | (k >> 1)) & 3), r ^ t
TABLE = [[INDEX[product(a, b)] for b in ELEMENTS] for a in ELEMENTS]
U = INDEX[2, 1]


def xor(xs, initial=0):
    for x in xs:
        initial ^= x
    return initial


def left(g, x):
    return xor(1 << TABLE[g][j] for j in range(8) if x >> j & 1)


def mul(a, b):
    return xor(left(i, b) for i in range(8) if a >> i & 1)


def replay(result):
    n = 4
    a = [[int(i == j) for j in range(n)] for i in range(n)]
    b = [[0]*n for _ in range(n)]
    for u, v, x, y in result['frames']:
        for row, col, wanted in ((y, u, 1), (y, x, 0), (v, u, 0), (v, x, 1)):
            if xor(mul(i, j) for i, j in zip(row, col)) != wanted:
                raise ValueError('The ordinary frame is not dual')
        def query(matrix):
            va = [xor(mul(v[j], matrix[j][k]) for j in range(n)) for k in range(n)]
            return [[mul(u[i], left(U, va[k])) for k in range(n)] for i in range(n)]
        def apply(matrix):
            ya = [xor(mul(y[j], matrix[j][k]) for j in range(n)) for k in range(n)]
            va = [xor(mul(v[j], matrix[j][k]) for j in range(n)) for k in range(n)]
            qa = query(matrix)
            return [[matrix[i][k] ^ mul(u[i], left(U, ya[k]) ^ ya[k])
                     ^ mul(x[i], left(U, va[k]) ^ va[k]) ^ qa[i][k]
                     for k in range(n)] for i in range(n)]
        next_a, next_b, source = apply(a), apply(b), query(a)
        a, b = next_a, [[x ^ y for x, y in zip(row, other)] for row, other in zip(next_b, source)]
    allowed = 192 if result['parameters']['known_control'] else 3
    if any((a[i][j] ^ b[i][j]) & ~allowed or b[i][j] & ~192 for i in range(n) for j in range(n)):
        raise ValueError('The full endpoint has an extra group coefficient')
    augmentation = [sum((((b[i][j] >> 6) ^ (b[i][j] >> 7)) & 1) << j for j in range(n)) for i in range(n)]
    if binary.rank(augmentation) != n:
        raise ValueError('The Gaussian target basis is singular')
    return a, b


def search(args):
    n = 4
    solver = z3.Solver()
    solver.set(timeout=args.timeout, random_seed=5786)
    zero = z3.BitVecVal(0, 8)
    def bit(x, j):
        return z3.Extract(j, j, x) == 1
    def mono(g, x):
        inverse = [TABLE[g].index(j) for j in range(8)]
        return z3.Concat(*(z3.Extract(inverse[j], inverse[j], x) for j in range(7, -1, -1)))
    def poly(p, x):
        return z3.If(bit(p, 0), x, zero) ^ z3.If(bit(p, 1), mono(1, x), zero)
    def poly_mul(p, q):
        return z3.If(bit(p, 0), q, z3.BitVecVal(0, 2)) ^ z3.If(bit(p, 1), z3.RotateLeft(q, 1), z3.BitVecVal(0, 2))
    a = [[z3.BitVecVal(i == j, 8) for j in range(n)] for i in range(n)]
    b = [[zero]*n for _ in range(n)]
    variables = []
    for stage in range(4):
        if stage == 0 or args.known_control:
            first, second = [(0, 1), (1, 0), (2, 3), (3, 2)][stage]
            u, v, x, y = [[z3.BitVecVal(int(i == j), 2) for i in range(n)] for j in (first, second, second, first)]
        else:
            u, v, x, y = [[z3.BitVec(f'{name}_{stage}_{i}', 2) for i in range(n)] for name in ('u', 'v', 'x', 'y')]
            if args.fixed_augmentation:
                first, second = [(0, 1), (1, 0), (2, 3), (3, 2)][stage]
                for row, pivot in zip((u, v, x, y), (first, second, second, first)):
                    solver.add([z3.Xor(bit(p, 0), bit(p, 1)) == (i == pivot) for i, p in enumerate(row)])
            for row, col, wanted in ((y, u, 1), (y, x, 0), (v, u, 0), (v, x, 1)):
                solver.add(xor((poly_mul(i, j) for i, j in zip(row, col)), z3.BitVecVal(0, 2)) == wanted)
        variables.append((u, v, x, y))
        def query(matrix):
            va = [xor((poly(v[j], matrix[j][k]) for j in range(n)), zero) for k in range(n)]
            return [[poly(u[i], mono(U, va[k])) for k in range(n)] for i in range(n)]
        def apply(matrix, label):
            ya = [xor((poly(y[j], matrix[j][k]) for j in range(n)), zero) for k in range(n)]
            va = [xor((poly(v[j], matrix[j][k]) for j in range(n)), zero) for k in range(n)]
            q = query(matrix)
            extra = query(a) if label == 'B' else [[zero]*n for _ in range(n)]
            result = [[z3.BitVec(f'{label}_{stage}_{i}_{k}', 8) for k in range(n)] for i in range(n)]
            solver.add([result[i][k] == (matrix[i][k] ^ poly(u[i], mono(U, ya[k]) ^ ya[k])
                        ^ poly(x[i], mono(U, va[k]) ^ va[k]) ^ q[i][k] ^ extra[i][k])
                        for i in range(n) for k in range(n)])
            return result
        next_b = apply(b, 'B')
        a, b = apply(a, 'A'), next_b
    allowed = 192 if args.known_control else 3
    solver.add([((a[i][j] ^ b[i][j]) & (255 ^ allowed)) == 0 for i in range(n) for j in range(n)])
    solver.add([(b[i][j] & 63) == 0 for i in range(n) for j in range(n)])
    augmentation = [z3.Concat(*(z3.Extract(6, 6, b[i][j]) ^ z3.Extract(7, 7, b[i][j]) for j in range(n-1, -1, -1))) for i in range(n)]
    inverse = [z3.BitVec(f'inverse_{i}', n) for i in range(n)]
    for i in range(n):
        solver.add(xor((z3.If(bit(augmentation[i], j), inverse[j], z3.BitVecVal(0, n)) for j in range(n)), z3.BitVecVal(0, n)) == 1 << i)
    start = time.monotonic()
    status = str(solver.check())
    result = {'status': status, 'solve_seconds': time.monotonic()-start,
              'parameters': vars(args) | {'output': str(args.output)}, 'group_elements': ELEMENTS,
              'scope': 'Binary dihedral group algebra with arbitrary ordinary swap-dependent full frames; no lift or kappa'}
    if status == 'unknown':
        result['reason'] = solver.reason_unknown()
    if status == 'sat':
        model = solver.model()
        result['frames'] = [[[model.eval(x).as_long() for x in row] for row in frame] for frame in variables]
        aa, bb = replay(result)
        result.update(a=aa, b=bb, verification='PASS_LITERAL_DIHEDRAL_RESIDUAL_REPLAY')
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('frames', 'a', 'b', 'group_elements')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--timeout', type=int, default=180000)
    p.add_argument('--known-control', action='store_true')
    p.add_argument('--fixed-augmentation', action='store_true', help='Fix the R=1 frame to two two-call words; keep all swap-dependent terms')
    p.add_argument('--output', type=Path, required=True)
    search(p.parse_args())
