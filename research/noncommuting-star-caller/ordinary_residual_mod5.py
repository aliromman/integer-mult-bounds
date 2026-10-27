"""Check the ordinary-residual contract over F5.

Dyadic Gaussian scalars map to F5 with i mapped to 2. A modular positive
word still needs a lift. A negative result applies only to the stated
four-field contract and to bases that remain invertible modulo five.
"""
import argparse
import json
from pathlib import Path
import time

import z3


def literal_check(result):
    n, length = 4, 5
    def mm(a, b):
        return [[sum(x*y for x, y in zip(row, col)) % 5 for col in zip(*b)] for row in a]
    def add(a, b, sign=1):
        return [[(x+sign*y) % 5 for x, y in zip(ar, br)] for ar, br in zip(a, b)]
    def rank(a):
        rows, used = [row[:] for row in a], 0
        for j in range(n):
            pivot = next((i for i in range(used, n) if rows[i][j]), None)
            if pivot is None:
                continue
            rows[used], rows[pivot] = rows[pivot], rows[used]
            inv = pow(rows[used][j], -1, 5)
            rows[used] = [x*inv % 5 for x in rows[used]]
            for i in range(n):
                if i != used:
                    rows[i] = [(x-rows[i][j]*y) % 5 for x, y in zip(rows[i], rows[used])]
            used += 1
        return used
    ident = [[int(i == j) for j in range(n)] for i in range(n)]
    zero = [[0]*n for _ in range(n)]
    channels = [[ident]+[zero]*(length-1) for _ in range(2)]
    for stage, (u, v, x, y) in enumerate(result['frames']):
        p = [[(u[i]*y[j]+x[i]*v[j]) % 5 for j in range(n)] for i in range(n)]
        q = [[u[i]*v[j] % 5 for j in range(n)] for i in range(n)]
        if rank(p) != 2 or rank(q) != 1 or mm(p, p) != p or mm(q, q) != zero or mm(p, q) != q or mm(q, p) != q:
            raise ValueError('The modular native frame fails')
        nxt = []
        for channel, a in enumerate(channels):
            pq = add(p, q, 1 if channel == 0 else -1)
            c, l = (pq, add(ident, p, -1)) if result['parameters']['inverse_mask'] >> stage & 1 else (add(ident, p, -1), pq)
            out = [zero for _ in range(length)]
            for d in range(stage+1):
                out[d] = add(out[d], mm(c, a[d]))
                out[d+1] = add(out[d+1], mm(l, a[d]))
            nxt.append(out)
        channels = nxt
    left, right = result['left'], result['right']
    if rank(left) != n or rank(right) != n or channels != result['channels']:
        raise ValueError('The target basis or polynomial replay fails')
    for d in range(length):
        middle = [[left[i][j] if result['normalized_target_degrees'][j] == d else 0 for j in range(n)] for i in range(n)]
        target = mm(middle, right)
        if add(channels[0][d], channels[1][d], -1) != add(target, target):
            raise ValueError('The target coefficient fails')
        if d != result['residual_degree'] and add(channels[0][d], channels[1][d]) != zero:
            raise ValueError('The residual is at the wrong degree')


def search(args):
    n, calls, length = 4, 4, 5
    shift = args.inverse_mask.bit_count()
    degrees = [2]*n if not shift else [0, 0, 4, 4]
    residual_degree = shift if args.residual_degree is None else args.residual_degree
    solver = z3.SolverFor('QF_BV')
    solver.set(timeout=args.timeout, random_seed=5786)
    def val(x):
        return z3.BitVecVal(x % 5, 3)
    zero, one = val(0), val(1)
    def add(a, b):
        c = z3.ZeroExt(1, a)+z3.ZeroExt(1, b)
        return z3.Extract(2, 0, z3.If(z3.UGE(c, 5), c-5, c))
    def neg(a):
        return z3.If(a == 0, zero, 5-a)
    def mul(a, b):
        twice = add(a, a)
        return z3.If(b == 0, zero, z3.If(b == 1, a, z3.If(b == 2, twice, z3.If(b == 3, neg(twice), neg(a)))))
    def total(xs):
        result = zero
        for x in xs:
            result = add(result, x)
        return result
    def fresh(label, rows=n, cols=n):
        a = [[z3.BitVec(f'{label}_{i}_{j}', 3) for j in range(cols)] for i in range(rows)]
        solver.add([z3.ULT(x, 5) for row in a for x in row])
        return a
    def mm(a, b):
        return [[total(mul(a[i][k], b[k][j]) for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
    def plus(a, b):
        return [[add(x, y) for x, y in zip(ar, br)] for ar, br in zip(a, b)]
    def minus(a, b):
        return [[add(x, neg(y)) for x, y in zip(ar, br)] for ar, br in zip(a, b)]
    def equal(a, b):
        solver.add([x == y for ar, br in zip(a, b) for x, y in zip(ar, br)])
    ident = [[val(i == j) for j in range(n)] for i in range(n)]
    blank = [[zero]*n for _ in range(n)]
    channels = [[ident]+[blank]*(length-1) for _ in range(2)]
    frames = []
    for stage in range(calls):
        if stage == 0 or args.known_control:
            first, second = [(0, 1), (1, 0), (2, 3), (3, 2)][stage]
            u, v, x, y = [[val(i == k) for i in range(n)] for k in (first, second, second, first)]
        else:
            u, v, x, y = fresh(f'frame_{stage}', 4, n)
            solver.add(total(mul(a, b) for a, b in zip(u, v)) == 0,
                       total(mul(a, b) for a, b in zip(v, x)) == 1,
                       total(mul(a, b) for a, b in zip(y, u)) == 1,
                       total(mul(a, b) for a, b in zip(y, x)) == 0)
        frames.append((u, v, x, y))
        p, q = fresh(f'P_{stage}'), fresh(f'N_{stage}')
        equal(p, [[add(mul(u[i], y[j]), mul(x[i], v[j])) for j in range(n)] for i in range(n)])
        equal(q, [[mul(u[i], v[j]) for j in range(n)] for i in range(n)])
        next_channels = []
        backwards = args.inverse_mask >> stage & 1
        for channel, a in enumerate(channels):
            pq = plus(p, q) if channel == 0 else minus(p, q)
            const, linear = (pq, minus(ident, p)) if backwards else (minus(ident, p), pq)
            result = [blank for _ in range(length)]
            for d in range(stage+1):
                result[d] = plus(result[d], mm(const, a[d]))
                result[d+1] = plus(result[d+1], mm(linear, a[d]))
            explicit = [fresh(f'word_{stage}_{channel}_{d}') for d in range(length)]
            for arow, brow in zip(explicit, result):
                equal(arow, brow)
            next_channels.append(explicit)
        channels = next_channels
    left, right = fresh('left'), fresh('right')
    for label, a in (('left', left), ('right', right)):
        inv = fresh(label+'_inverse')
        equal(mm(a, inv), ident)
    for d in range(length):
        middle = [[left[i][j] if degrees[j] == d else zero for j in range(n)] for i in range(n)]
        target = mm(middle, right)
        equal(minus(channels[0][d], channels[1][d]), plus(target, target))
        if d != residual_degree:
            equal(plus(channels[0][d], channels[1][d]), blank)
    start = time.monotonic()
    status = str(solver.check())
    result = {'status': status, 'solve_seconds': time.monotonic()-start,
              'parameters': vars(args) | {'output': str(args.output)},
              'scope': 'Four-field F5 Laurent screen with a fixed first native frame; no kappa claim',
              'residual_degree': residual_degree, 'normalized_target_degrees': degrees}
    if status == 'unknown':
        result['reason'] = solver.reason_unknown()
    if status == 'sat':
        model = solver.model()
        def evaluate(a):
            return [[model.eval(x).as_long() for x in row] for row in a]
        result.update(frames=[evaluate(a) for a in frames], left=evaluate(left), right=evaluate(right),
                      channels=[[evaluate(a) for a in channel] for channel in channels],
                      verification='PASS_LITERAL_MOD5_POLYNOMIAL_REPLAY')
        literal_check(result)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('frames', 'left', 'right', 'channels')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inverse-mask', type=int, default=0)
    p.add_argument('--residual-degree', type=int, help='Control only: allow the residual at another normalized degree')
    p.add_argument('--timeout', type=int, default=60000)
    p.add_argument('--known-control', action='store_true', help='Fix two independent two-call words; residual is not ordinary')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.inverse_mask.bit_count() not in (0, 2) or args.inverse_mask >= 16:
        p.error('Use four forward calls or two forward and two inverse calls')
    if args.residual_degree is not None and not 0 <= args.residual_degree < 5:
        p.error('The residual degree must be from zero through four')
    if args.known_control and (args.inverse_mask or args.residual_degree != 2):
        p.error('The known control uses forward calls and residual degree two')
    search(args)
