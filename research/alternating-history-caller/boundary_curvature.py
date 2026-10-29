"""Test boundary projectors against the second monomial coefficient.

The target stage powers are not fixed. Any constant monomial endpoint
satisfies W_2=(H^2-H)/2, with H=W_1 in the variable z-1. This is only a
necessary condition. Ordered variables permit arbitrary coupled widths.
"""
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'initialized-query-work'))
import exact_replay as e
import free_replay as nc
from discovery import construct


def fixed(n):
    out = []
    for site in range(1, n):
        if site % 2 == 0:
            out.append(e.zero(n))
            continue
        u = [(-1)**j if j < site else 0 for j in range(n)]
        v = [1 if j < site else (-1 if j > site else 0) for j in range(n)]
        out.append(e.numeric([[u[i]*v[j]+int(i == site)*int(j >= site)
                               for j in range(n)] for i in range(n)]))
    return out


def check(n):
    started = time.monotonic()
    names, ps, _, _ = construct(n, 0, tuple(range(1, n)))
    ps = [ps[0]]+fixed(n)+[ps[-1]]
    h, second = e.zero(n), e.zero(n)
    for p in ps:
        second = e.ma(second, nc.mm(p, h))
        h = e.ma(h, p)
    curvature = e.ma(e.ma(nc.mm(h, h), h, -1),
                     [[e.scale(p, 2) for p in row] for row in second], -1)
    system = [p for row in curvature for p in row]
    system += [p for end in [ps[0], ps[-1]]
               for row in e.ma(nc.mm(end, end), end, -1) for p in row]
    steps = []
    while True:
        reduced = e.row_reduce(system)
        affine = [p for p in reduced if all(len(m) <= 1 for m in p)]
        row = {'variables': len({v for p in system for m in p for v in m}),
               'rank': len(reduced), 'affine_consequences': len(affine)}
        steps.append(row)
        print(json.dumps(row), flush=True)
        if any(set(p) == {()} for p in affine):
            status = 'EXACT_BOUNDARY_CURVATURE_CONTRADICTION'
            break
        if not affine:
            status = 'UNRESOLVED_BOUNDARY_CURVATURE_SYSTEM' if any(system) else 'SECOND_COEFFICIENT_CONSISTENT'
            break
        values = e.linear_solution(affine)
        system = [nc.substitute(p, values) for p in system]
    return {'n': n, 'status': status, 'steps': steps,
            'seconds': time.monotonic()-started, 'new_kappa': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, choices=(4, 6, 8, 10, 12, 16), default=6)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = check(args.n)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
