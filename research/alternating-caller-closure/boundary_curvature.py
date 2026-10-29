"""Build the complete boundary charts and the fixed interior.

The chart is retained from the initialized-query investigation. Only the
zero-work, thin-interior specialization is used by this package.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).with_name('inputs')))
import exact_replay as e
import free_replay as nc


def construct(r, work, thin=()):
    n = r + work
    oldus, oldvs, _ = e.canonical(r)
    us = [u + [0]*work for u in oldus]
    vs = [v + [0]*work for v in oldvs]
    names, ps, qs, ts = [], [], [], []

    def variable(name):
        index = len(names)
        names.append(name)
        return {(index,): F(1)}

    for i, (u, v) in enumerate(zip(us, vs)):
        ts.append(e.numeric([[int(a == b)-2*u[a]*v[b] for b in range(n)]
                             for a in range(n)]))
        if i in thin:
            ps.append(e.zero(n))
            continue
        pivot = next(j for j in range(n) if v[j])
        drop = next(j for j in range(n) if u[j] and j != pivot)
        h = [F(int(j == pivot), v[pivot]) for j in range(n)]
        columns = [u, h] + [[int(k == j)-v[j]*h[k] for k in range(n)]
                            for j in range(n) if j not in (pivot, drop)]
        basis = list(map(list, zip(*columns)))
        av = [variable(f'a{i}_{j}') for j in range(n-2)]
        bv = [variable(f'b{i}_{j}') for j in range(n-2)]
        c = variable(f'c{i}')
        core = e.zero(n)
        core[0][0] = core[1][1] = e.const(1)
        core[0][1] = c
        for j in range(n-2):
            core[j+2][1], core[0][j+2] = av[j], bv[j]
        ps.append(e.mm(e.mm(e.numeric(basis), core), e.numeric(e.inverse(basis))))
        qs.append(e.add(c, *(e.mul(a, b) for a, b in zip(av, bv))))
    return names, ps, qs, ts


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

