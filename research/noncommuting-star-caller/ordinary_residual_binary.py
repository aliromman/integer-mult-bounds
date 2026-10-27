"""Exact binary screen for four full calls and an ordinary input residual.

The query is S=1+e with e^2=0. Pure and query coefficients are kept at
every Laurent degree. A positive result needs a characteristic-zero lift.
"""
import argparse
import json
from pathlib import Path
import time

import z3
import screen as s


def run(args):
    n, calls, period = 4, 4, args.period
    shift = args.inverse_mask.bit_count()
    degrees = [2]*n if not shift else [-2, -2, 2, 2]
    degrees = [(d+shift) % period for d in degrees]
    solver = z3.Solver()
    solver.set(timeout=args.timeout, random_seed=5786)
    zero = z3.BitVecVal(0, n)
    def bit(x, j):
        return z3.Extract(j, j, x) == 1
    def xor(xs, base=zero):
        for x in xs:
            base ^= x
        return base
    def select(rows, mask):
        return xor(z3.If(bit(mask, j), row, zero) for j, row in enumerate(rows))
    def parity(x):
        value = bit(x, 0)
        for j in range(1, n):
            value = z3.Xor(value, bit(x, j))
        return value
    a = [[z3.BitVecVal(1 << i, n) if d == 0 else zero for i in range(n)] for d in range(period)]
    b = [[zero]*n for _ in range(period)]
    variables = []
    for stage in range(calls):
        if stage == 0:
            u, v, x, y = [z3.BitVecVal(k, n) for k in (1, 2, 2, 1)]
        else:
            u, v, x, y = [z3.BitVec(f'{name}_{stage}', n) for name in ('u', 'v', 'x', 'y')]
            solver.add(z3.Not(parity(u & v)), parity(v & x), parity(y & u), z3.Not(parity(y & x)))
        variables.append((u, v, x, y))
        backwards = args.inverse_mask >> stage & 1
        aa = [[zero]*n for _ in range(period)]
        bb = [[zero]*n for _ in range(period)]
        for d in range(stage+1):
            va, vb = select(a[d], v), select(b[d], v)
            ya, yb = select(a[d], y), select(b[d], y)
            for i in range(n):
                pa = z3.If(bit(u, i), ya, zero) ^ z3.If(bit(x, i), va, zero)
                pb = z3.If(bit(u, i), yb, zero) ^ z3.If(bit(x, i), vb, zero)
                na = z3.If(bit(u, i), va, zero)
                nb = z3.If(bit(u, i), vb, zero)
                if backwards:
                    aa[d][i] ^= pa ^ na
                    bb[d][i] ^= pb ^ nb ^ na
                    aa[(d+1) % period][i] ^= a[d][i] ^ pa
                    bb[(d+1) % period][i] ^= b[d][i] ^ pb
                else:
                    aa[d][i] ^= a[d][i] ^ pa
                    bb[d][i] ^= b[d][i] ^ pb
                    aa[(d+1) % period][i] ^= pa ^ na
                    bb[(d+1) % period][i] ^= pb ^ nb ^ na
        anew = [[z3.BitVec(f'A_{stage}_{d}_{i}', n) for i in range(n)] for d in range(period)]
        bnew = [[z3.BitVec(f'B_{stage}_{d}_{i}', n) for i in range(n)] for d in range(period)]
        solver.add([anew[d][i] == aa[d][i] for d in range(period) for i in range(n)])
        solver.add([bnew[d][i] == bb[d][i] for d in range(period) for i in range(n)])
        a, b = anew, bnew
    left, right = [[z3.BitVec(f'{label}_{i}', n) for i in range(n)] for label in ('left', 'right')]
    for label, matrix in (('left', left), ('right', right)):
        inverse = [z3.BitVec(f'{label}_inverse_{i}', n) for i in range(n)]
        solver.add([select(inverse, matrix[i]) == 1 << i for i in range(n)])
    for d in range(period):
        selector = sum((e == d) << j for j, e in enumerate(degrees))
        solver.add([b[d][i] == select(right, left[i] & selector) for i in range(n)])
        if d != shift:
            solver.add([a[d][i] == b[d][i] for i in range(n)])
    start = time.monotonic()
    status = str(solver.check())
    result = {'status': status, 'solve_seconds': time.monotonic()-start,
              'parameters': vars(args) | {'output': str(args.output)},
              'scope': 'Four-field binary Laurent screen; first full frame fixed by conjugation; no native or kappa claim'}
    if status == 'unknown':
        result['reason'] = solver.reason_unknown()
    if status == 'sat':
        model = solver.model()
        witness = [[model.eval(x).as_long() for x in row] for row in variables]
        aa = [s.eye(n)]+[(0,)*n]*(period-1)
        bb = [(0,)*n]*period
        for stage, (u, v, x, y) in enumerate(witness):
            p = tuple((y if u >> i & 1 else 0) ^ (v if x >> i & 1 else 0) for i in range(n))
            q = tuple(v if u >> i & 1 else 0 for i in range(n))
            if s.mm(p, p) != p or s.rank(p) != 2 or s.rank(q) != 1 or any(s.mm(q, q)) or s.mm(p, q) != q or s.mm(q, p) != q:
                raise ValueError('The native projector or query is invalid')
            backwards = args.inverse_mask >> stage & 1
            const = s.mxor(p, q) if backwards else s.mxor(s.eye(n), p)
            linear = s.mxor(s.eye(n), p) if backwards else s.mxor(p, q)
            ca, cb = [(0,)*n]*period, [(0,)*n]*period
            for d in range(stage+1):
                ca[d] = s.mxor(ca[d], s.mm(const, aa[d]))
                cb[d] = s.mxor(cb[d], s.mm(const, bb[d]))
                ca[(d+1) % period] = s.mxor(ca[(d+1) % period], s.mm(linear, aa[d]))
                cb[(d+1) % period] = s.mxor(cb[(d+1) % period], s.mm(linear, bb[d]))
                qdegree = d if backwards else (d+1) % period
                cb[qdegree] = s.mxor(cb[qdegree], s.mm(q, aa[d]))
            aa, bb = ca, cb
        l, r = [tuple(model.eval(x).as_long() for x in row) for row in (left, right)]
        if s.rank(l) != n or s.rank(r) != n:
            raise ValueError('A target basis is singular')
        for d in range(period):
            diagonal = tuple(1 << i if degrees[i] == d else 0 for i in range(n))
            if bb[d] != s.mm(s.mm(l, diagonal), r) or (d != shift and aa[d] != bb[d]):
                raise ValueError('The literal target or residual fails')
        result.update(frames=witness, left=l, right=r, a_coefficients=aa, b_coefficients=bb,
                      ordinary_residual=s.mxor(aa[shift], bb[shift]), target_degrees=degrees,
                      normalization_degree=shift)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('a_coefficients', 'b_coefficients')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inverse-mask', type=int, default=0)
    p.add_argument('--period', type=int, choices=(4, 5), default=5, help='Four is cyclic; five retains every formal degree')
    p.add_argument('--timeout', type=int, default=60000)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.inverse_mask.bit_count() not in (0, 2) or args.inverse_mask >= 16:
        p.error('Use four forward calls or two forward and two inverse calls')
    run(args)
