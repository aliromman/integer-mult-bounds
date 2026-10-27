"""Search a shortest binary query frame with a finite-order stage.

Cinna Davis supplies the simplex and rank-two projector parametrization.
The earlier formal-stage search is credited in SOURCE-PINS.json. This test
uses the relation D^period=1. It does not establish a native caller.
"""
import argparse
import json
from pathlib import Path
import time

import z3
import screen as s


def mixed_frame(types):
    """Factor a corank-two Gram matrix in six-type coordinates.

    The six-type model is retained from the released 2026-10-23 controller
    investigation. Here it supplies a common-stage polynomial search.
    """
    pairs = [(a, b) for a in range(1, 4) for b in range(1, 4)
             if (a & b).bit_count() % 2]
    aa, bb = zip(*(pairs[t] for t in types))
    gram = tuple((1 << i) ^ sum(((a & b).bit_count() % 2) << j
                                for j, b in enumerate(bb)) for i, a in enumerate(aa))
    n, size = len(types)-2, len(types)
    if s.mm(gram, gram) != gram or s.rank(gram) != n:
        raise ValueError('The type sequence does not give a corank-two projector')
    columns = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(gram))
                    for j in range(size))
    basis = []
    for col in columns:
        if s.rank(basis+[col]) > len(basis):
            basis.append(col)
    coords = {0: 0}
    for j, col in enumerate(basis):
        coords.update({key ^ col: value | (1 << j) for key, value in list(coords.items())})
    u = [coords[col] for col in columns]
    v = [sum(((col >> i) & 1) << j for j, col in enumerate(basis)) for i in range(size)]
    total = (0,)*n
    for col, row in zip(u, v):
        if not col or not row or (col & row).bit_count() % 2:
            raise ValueError('A query is zero or is not square-zero')
        total = s.mxor(total, tuple(row if col >> i & 1 else 0 for i in range(n)))
    if total != s.eye(n):
        raise ValueError('The full query sum is not the identity')
    return list(zip(u, v))


def search(args):
    blocks = [args.n] if args.blocks is None else [int(x) for x in args.blocks.split(',')]
    if any(x < 2 or x % 2 for x in blocks):
        raise ValueError('Each simplex block must have positive even dimension')
    n = sum(blocks)
    query_blocks, offset = [], 0
    for size in blocks:
        mask = ((1 << size)-1) << offset
        query_blocks.append([(1 << (offset+i), mask ^ (1 << (offset+i))) for i in range(size)]+[(mask, mask)])
        offset += size
    if args.order == 'alternating':
        query_frame = [block[i] for i in range(max(map(len, query_blocks))) for block in query_blocks if i < len(block)]
    else:
        query_frame = [q for block in query_blocks for q in block]
    if args.types:
        query_frame = mixed_frame([int(t) for t in args.types.split(',')])
        n = len(query_frame)-2
    calls = len(query_frame)
    period = args.period or calls+1
    inverse_count = args.inverse_mask.bit_count()
    degrees = [int(x)+inverse_count for x in args.degrees.split(',')] if args.degrees else None
    if degrees and (len(degrees) != n or min(degrees) < 0 or max(degrees) >= period):
        raise ValueError('The declared normalized target degrees do not fit')
    target_rank = 2*(calls-n)
    solver = z3.Solver()
    solver.set(timeout=args.timeout, random_seed=args.seed)
    zero = z3.BitVecVal(0, n)
    def bit(value, i):
        return z3.Extract(i, i, value) == 1
    def parity(value):
        out = bit(value, 0)
        for j in range(1, n):
            out = z3.Xor(out, bit(value, j))
        return out
    def selected_sum(rows, selector):
        out = zero
        for j, row in enumerate(rows):
            out = out ^ z3.If(bit(selector, j), row, zero)
        return out
    def fresh(label, stage):
        return [[z3.BitVec(f'{label}_{stage}_{d}_{i}', n) for i in range(n)] for d in range(period)]
    aa = [[z3.BitVecVal(1 << i, n) if d == 0 else zero for i in range(n)] for d in range(period)]
    bb = [[zero for _ in range(n)] for _ in range(period)]
    projectors = []
    for stage, (u, v) in enumerate(query_frame):
        backwards = bool(args.inverse_mask >> stage & 1)
        a, b = z3.BitVec(f'a_{stage}', n), z3.BitVec(f'b_{stage}', n)
        solver.add(parity(a & v), parity(b & u), z3.Not(parity(a & b)))
        solver.add(z3.Not(bit(a, (u & -u).bit_length()-1)))
        projectors.append((u, v, a, b))
        next_a, next_b = fresh('A', stage), fresh('B', stage)
        ca = [[zero for _ in range(n)] for _ in range(period)]
        cb = [[zero for _ in range(n)] for _ in range(period)]
        for d in range(period):
            va, vb = selected_sum(aa[d], z3.BitVecVal(v, n)), selected_sum(bb[d], z3.BitVecVal(v, n))
            ba, bbb = selected_sum(aa[d], b), selected_sum(bb[d], b)
            for i in range(n):
                pa = (ba if u >> i & 1 else zero) ^ z3.If(bit(a, i), va, zero)
                pb = (bbb if u >> i & 1 else zero) ^ z3.If(bit(a, i), vb, zero)
                qa = va if u >> i & 1 else zero
                if args.shifted:
                    ca[d][i] = ca[d][i] ^ aa[d][i]
                    cb[d][i] = cb[d][i] ^ bb[d][i] ^ qa
                else:
                    ca[d][i] = ca[d][i] ^ (pa if backwards else aa[d][i] ^ pa)
                    cb[d][i] = cb[d][i] ^ (pb ^ qa if backwards else bb[d][i] ^ pb)
                if d+1 < period or not args.shifted:
                    ca[(d+1) % period][i] = ca[(d+1) % period][i] ^ (aa[d][i] ^ pa if backwards else pa)
                    cb[(d+1) % period][i] = cb[(d+1) % period][i] ^ (bb[d][i] ^ pb if backwards else pb ^ (va if u >> i & 1 else zero))
        solver.add([next_a[d][i] == ca[d][i] for d in range(period) for i in range(n)])
        solver.add([next_b[d][i] == cb[d][i] for d in range(period) for i in range(n)])
        aa, bb = next_a, next_b
    solver.add([aa[d][i] == bb[d][i] for d in range(period) for i in range(n)])
    if degrees:
        left = [z3.BitVec(f'left_{i}', n) for i in range(n)]
        right = [z3.BitVec(f'right_{i}', n) for i in range(n)]
        for i in range(n):
            solver.add(selected_sum(right, left[i]) == (1 << i))
            for d in range(period):
                selector = sum(((d & degree) == d if args.shifted else degree == d) << j
                               for j, degree in enumerate(degrees))
                solver.add(aa[d][i] == selected_sum(right, left[i] & selector))
    if args.balanced:
        if period < 4 or target_rank > n:
            raise ValueError('The balanced target does not fit this frame')
        left = [z3.BitVec(f'left_{i}', target_rank) for i in range(n)]
        right = [z3.BitVec(f'right_{i}', n) for i in range(target_rank)]
        for i in range(target_rank):
            for j in range(target_rank):
                terms = [z3.And(bit(right[i], k), bit(left[k], j)) for k in range(n)]
                value = terms[0]
                for term in terms[1:]:
                    value = z3.Xor(value, term)
                solver.add(value == (i == j))
        for i in range(n):
            p = zero
            for j in range(target_rank):
                p = p ^ z3.If(bit(left[i], j), right[j], zero)
            for d in range(period):
                want = (p if d == 3 else (p ^ (1 << i) if d == 2 else zero))
                if args.shifted:
                    want = p if d in (1, 3) else ((1 << i) if d in (0, 2) else zero)
                solver.add(aa[d][i] == want)
    start = time.monotonic()
    answer = str(solver.check())
    result = {'parameters': vars(args) | {'output': str(args.output)},
              'status': answer, 'solve_seconds': time.monotonic()-start,
              'scope': 'Binary multi-simplex query frame; no characteristic-zero lift or complete price',
              'query_frame': query_frame, 'effective_period': period,
              'formal_stage': args.period == 0, 'normalizing_degree': inverse_count}
    if answer == 'unknown':
        result['reason'] = solver.reason_unknown()
    if answer == 'sat':
        model = solver.model()
        ps = [tuple((model.eval(b).as_long() if u >> i & 1 else 0)
                    ^ (v if model.eval(a).as_long() >> i & 1 else 0) for i in range(n))
              for u, v, a, b in projectors]
        # This replay uses literal polynomial matrices, not solver expressions.
        a = [s.eye(n)]+[(0,)*n for _ in range(period-1)]
        b = [(0,)*n for _ in range(period)]
        for stage, p in enumerate(ps):
            backwards = bool(args.inverse_mask >> stage & 1)
            u, v = projectors[stage][:2]
            query = tuple(v if u >> i & 1 else 0 for i in range(n))
            if s.mm(p, p) != p or s.rank(p) != 2 or s.mm(p, query) != query or s.mm(query, p) != query:
                raise ValueError('Invalid native projector')
            ca, cb = [(0,)*n for _ in range(period)], [(0,)*n for _ in range(period)]
            for d in range(period):
                pa, pb = s.mm(p, a[d]), s.mm(p, b[d])
                qa = s.mm(query, a[d])
                ca[d] = s.mxor(ca[d], pa if backwards else s.mxor(a[d], pa))
                ca[(d+1) % period] = s.mxor(ca[(d+1) % period], s.mxor(a[d], pa) if backwards else pa)
                cb[d] = s.mxor(cb[d], s.mxor(pb, qa) if backwards else s.mxor(b[d], pb))
                cb[(d+1) % period] = s.mxor(cb[(d+1) % period], s.mxor(b[d], pb) if backwards else s.mxor(pb, qa))
            a, b = ca, cb
        if a != b:
            raise ValueError('Query endpoint fails the literal replay')
        if degrees:
            total = (0,)*n
            for d, block in enumerate(a):
                if s.rank(block) != degrees.count(d) or s.mm(block, block) != block:
                    raise ValueError('The target coefficient rank or projector fails')
                total = s.mxor(total, block)
                for e, other in enumerate(a):
                    if d != e and any(s.mm(block, other)):
                        raise ValueError('Target coefficient images are not independent')
            if total != s.eye(n):
                raise ValueError('Target coefficients do not sum to the identity')
        if args.balanced:
            if (any(any(a[d]) for d in range(period) if d not in (2, 3))
                    or s.mxor(a[2], a[3]) != s.eye(n) or s.mm(a[3], a[3]) != a[3]
                    or s.rank(a[3]) != target_rank):
                raise ValueError('Balanced target fails the literal replay')
        result.update(projectors=ps, coefficients=a, coefficient_ranks=list(map(s.rank, a)))
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['projectors', 'coefficients']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--period', type=int, default=4)
    p.add_argument('--blocks', help='Comma-separated even dimensions for a direct sum of query simplexes')
    p.add_argument('--types', help='Comma-separated six-type labels for a corank-two query frame')
    p.add_argument('--shifted', action='store_true', help='Solve in z=D+1 coordinates; the replay uses D')
    p.add_argument('--inverse-mask', type=int, default=0, help='Positions with inverse full stages')
    p.add_argument('--degrees', help='Comma-separated Laurent degrees of independent main targets')
    p.add_argument('--order', choices=('alternating', 'serial'), default='alternating')
    p.add_argument('--balanced', action='store_true')
    p.add_argument('--timeout', type=int, default=60000)
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 2 or args.n % 2 or args.period == 1 or args.period < 0:
        p.error('Use a positive even field count; period zero means a formal stage')
    if args.shifted and args.period not in (0, 4):
        p.error('Shifted coordinates require period zero or four')
    if args.inverse_mask and args.period:
        p.error('Inverse-stage targets currently use formal polynomial coordinates')
    if args.degrees and args.balanced:
        p.error('Use either explicit degrees or the balanced template')
    search(args)
