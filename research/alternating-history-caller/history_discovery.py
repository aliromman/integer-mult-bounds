"""Find pure native monomial rows in the complete saved prefix span.

Use coefficients of z, not coefficients of z-1. Saved values are correlated
copies. Only independent original-input rows can earn productive credit.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
from exact_replay import canonical


def projectors(n):
    out = []
    for site in range(1, n):
        if site % 2 == 0:
            out.append(s.zeros(n))
            continue
        u = s.Matrix([(-1)**j if j < site else 0 for j in range(n)])
        v = s.Matrix([[1 if j < site else (-1 if j > site else 0) for j in range(n)]])
        unit = s.eye(n)[:, site]
        tail = s.Matrix([[int(j >= site) for j in range(n)]])
        out.append(u*v+unit*tail)
    return [out[0]] + out + [out[-1]]


def advance(word, p, t):
    old = [t*x for x in word]
    complement = s.eye(p.rows)-p
    return [complement*old[0]] + [complement*old[i]+p*old[i-1]
                                 for i in range(1, len(old))] + [p*old[-1]]


def search(args):
    started = time.monotonic()
    n = args.n
    ps = projectors(n)
    if args.ends in ('thin', 'last'):
        ps[0] = s.zeros(n)
    if args.ends in ('thin', 'first'):
        ps[-1] = s.zeros(n)
    ts = [s.Matrix(t) for t in canonical(n)[2]]
    ident = s.eye(n)
    plus, minus = [ident], [ident]
    records = [(plus.copy(), minus.copy())]
    for p, t in zip(ps, ts):
        if p*p != p or p*t != t*p:
            raise ValueError('The projector control fails')
        plus, minus = advance(plus, p, ident), advance(minus, p, t)
        records.append((plus, minus))
    if any(a+b != s.zeros(n) for a, b in zip(plus, minus)):
        raise ValueError('The full query endpoint fails')
    degree = args.period-1 if args.period else max(i for i, a in enumerate(plus) if a != s.zeros(n))
    rows = []
    for plus, minus in records:
        both = []
        for word in [plus, minus]:
            folded = [s.zeros(n) for _ in range(degree+1)]
            for index, value in enumerate(word):
                if index <= degree or args.period:
                    folded[index % (degree+1)] += value
                elif value != s.zeros(n):
                    raise ValueError('The formal history degree is too large')
            both.append(folded)
        for row in range(n):
            rows.append([x/2 for sign in [1, -1] for i in range(degree+1)
                         for x in (both[0][i]+sign*both[1][i]).row(row)])
    histories = s.Matrix(s.Matrix(rows).rowspace())
    spaces = []
    for query in [True, False]:
        for exponent in range(degree+1):
            offset = (degree+1 if query else 0)+exponent
            good = list(range(offset*n, (offset+1)*n))
            bad = [j for j in range(histories.cols) if j not in good]
            reduced, _ = histories[:, bad+good].rref()
            vectors = [row[-n:] for row in reduced.tolist()
                       if not any(row[:-n]) and any(row[-n:])]
            spaces.append({'central_query': query, 'degree': exponent,
                           'dimension': len(vectors),
                           'basis': [[str(x) for x in row] for row in vectors]})
            print(json.dumps({k: v for k, v in spaces[-1].items() if k != 'basis'}), flush=True)
    full = sum(p != s.zeros(n) for p in ps)
    result = {'n': n, 'period': args.period, 'boundary_stages': args.ends,
              'full_calls': full, 'thin_calls': n+1-full,
              'coefficient_variable': 'z', 'saved_history_rows': histories.rows,
              'spaces': spaces, 'seconds': time.monotonic()-started, 'new_kappa': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'spaces'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=6)
    parser.add_argument('--period', type=int, default=0)
    parser.add_argument('--ends', choices=('full', 'thin', 'first', 'last'), default='full')
    parser.add_argument('--output', type=Path, required=True)
    search(parser.parse_args())
