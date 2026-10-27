"""Enumerate ordinary-monomial times one-Hadamard advancing cells.

The rank of the upper-right symplectic block is the Gaussian rank.
This finite census is a component screen, not a complete native caller.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path


def rank(rows):
    pivots = {}
    for x in rows:
        while x:
            k = x.bit_length()-1
            if k not in pivots:
                pivots[k] = x
                break
            x ^= pivots[k]
    return len(pivots)


def mm(a, b):
    out = []
    for row in a:
        value = 0
        for j, other in enumerate(b):
            if row >> j & 1:
                value ^= other
        out.append(value)
    return out


def inverse(a):
    n = len(a)
    rows = [row | 1 << (n+i) for i, row in enumerate(a)]
    for j in range(n):
        k = next(i for i in range(j, n) if rows[i] >> j & 1)
        rows[k], rows[j] = rows[j], rows[k]
        for i in range(n):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return [row >> n for row in rows]


def transpose(a):
    return [sum(((row >> j) & 1) << i for i, row in enumerate(a)) for j in range(len(a))]


def cell(a, c):
    n = len(a)
    ai = transpose(inverse(a))
    ac = mm(ai, c)
    r = a+[x | y << n for x, y in zip(ac, ai)]
    return [x ^ (((x ^ (x >> n)) & 1)*((1 << n) | 1)) for x in r]


def profile(u):
    n = len(u)//2
    ident = [1 << i for i in range(2*n)]
    p, ranks = ident, []
    for step in range(1, 121):
        p = mm(u, p)
        ranks.append(rank([x >> n for x in p[:n]]))
        if p == ident:
            return ranks
    raise ValueError('The symplectic order exceeds the declared bound')


def phase_power(a, c, order):
    """Multiply the literal unnormalized complex cell using integer pairs."""
    n, size = len(a), 1 << len(a)
    def mul(z, w):
        return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]
    phases = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    u = [[(0, 0) for _ in range(size)] for _ in range(size)]
    for j in range(size):
        for bit in range(2):
            x = (j & ~1) | bit
            exponent = sum(((c[k] >> k) & 1)*((x >> k) & 1) for k in range(n))
            exponent += 2*sum(((c[k] >> ell) & 1)*((x >> k) & 1)*((x >> ell) & 1)
                              for k in range(n) for ell in range(k+1, n))
            y = sum(((row & x).bit_count() % 2) << k for k, row in enumerate(a))
            z = phases[exponent % 4]
            u[y][j] = tuple(-v if (j & 1) and bit else v for v in z)
    p = [[(int(i == j), 0) for j in range(size)] for i in range(size)]
    support_ranks = []
    for _ in range(order):
        q = [[(0, 0) for _ in range(size)] for _ in range(size)]
        for i, j, k in itertools.product(range(size), repeat=3):
            z = mul(u[i][k], p[k][j])
            q[i][j] = q[i][j][0]+z[0], q[i][j][1]+z[1]
        p = q
        support_sizes = {sum(z != (0, 0) for z in row) for row in p}
        if len(support_sizes) != 1:
            raise ValueError('The literal power has unequal row support sizes')
        count = support_sizes.pop()
        if not count or count & (count-1):
            raise ValueError('The literal power does not have a power-of-two support')
        support_ranks.append(count.bit_length()-1)
    # The symplectic identity can still leave a Pauli. Keep that distinction.
    support = [[j for j, z in enumerate(row) if z != (0, 0)] for row in p]
    scalar = all(support[i] == [i] and p[i][i] == p[0][0] for i in range(size))
    eighth = (1, 0)
    if scalar:
        for _ in range(8):
            eighth = mul(eighth, p[0][0])
    return {'scalar': scalar, 'unnormalized_diagonal': p[0][0] if scalar else None,
            'denominator_sqrt2_power': order, 'matrix': p if not scalar else None,
            'literal_support_ranks': support_ranks,
            'eighth_tensor_power_has_phase_one': scalar and eighth == (1 << (4*order), 0)}


def census(n):
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    counts, examples = Counter(), {}
    for a in itertools.product(range(1 << n), repeat=n):
        if rank(a) != n:
            continue
        for mask in range(1 << len(pairs)):
            c = [0]*n
            for k, (i, j) in enumerate(pairs):
                if mask >> k & 1:
                    c[i] |= 1 << j
                    c[j] |= 1 << i
            u = cell(list(a), c)
            ranks = profile(u)
            key = tuple(ranks)
            counts[key] += 1
            examples.setdefault(key, (list(a), c, u))
    rows = []
    for key in sorted(counts, key=lambda k: (len(k), k)):
        a, c, u = examples[key]
        phase_check = phase_power(a, c, len(key))
        if phase_check['literal_support_ranks'] != list(key):
            raise ValueError('The literal matrix and symplectic rank profiles disagree')
        rows.append({'order': len(key), 'ranks': key, 'count': counts[key],
                     'address_matrix': a, 'phase_matrix': c, 'symplectic': u,
                     'phase_check': phase_check})
    return {'coordinates': n, 'cell_count': sum(counts.values()), 'profiles': rows}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coordinates', type=int, choices=(1, 2, 3), default=3)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = census(args.coordinates)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cell_count': result['cell_count'], 'profiles': [
        {k: x[k] for k in ('order', 'ranks', 'count')} for x in result['profiles']]}))
