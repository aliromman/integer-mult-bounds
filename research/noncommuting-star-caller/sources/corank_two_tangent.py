#!/usr/bin/env python3
"""Screen corank-two binary controller tangents. Apache-2.0; Codex assistance."""
import argparse
from collections import Counter
from itertools import product
import json
from random import Random


def parity(x):
    return x.bit_count() & 1


def rank(rows):
    basis = {}
    for row in rows:
        while row:
            k = row.bit_length() - 1
            if k not in basis:
                basis[k] = row
                break
            row ^= basis[k]
    return len(basis)


def mul(a, b):
    out = []
    for row in a:
        x = 0
        while row:
            bit = row & -row
            x ^= b[bit.bit_length() - 1]
            row ^= bit
        out.append(x)
    return tuple(out)


def outer(col, row, size):
    return tuple(row if col >> i & 1 else 0 for i in range(size))


def plus(*mats):
    return tuple(xor(items) for items in zip(*mats))


def xor(values):
    out = 0
    for value in values:
        out ^= value
    return out


TYPES = [(a, b) for a in range(1, 4) for b in range(1, 4) if parity(a & b)]
SMALL = [tuple(a if b >> j & 1 else 0 for j in range(2)) for a, b in TYPES]


def valid(word):
    return plus(*(SMALL[t] for t in word)) == (1, 2)


def data(word):
    size = len(word)
    aa, bb = zip(*(TYPES[t] for t in word))
    gram = tuple((1 << i) ^ sum(parity(a & bb[j]) << j for j in range(size)) for i, a in enumerate(aa))
    if mul(gram, gram) != gram or rank(gram) != size - 2 or any(row >> i & 1 for i, row in enumerate(gram)):
        raise ValueError('The corank-two Gram conditions failed.')
    suffix = (0,) * size
    candidates, ambiguous, excluded = {}, [], []
    for i in reversed(range(size)):
        u = sum((gram[j] >> i & 1) << j for j in range(size))
        v = gram[i]
        ri = outer(u, v, size)
        w = sum(parity(row & u) << j for j, row in enumerate(suffix))
        z = mul((v,), suffix)[0]
        if parity(v & w):
            e = plus(outer(u, z, size), outer(w, v, size))
            if parity(z & w):
                e = plus(e, ri)
            if rank(e) != 2 or mul(e, e) != e or mul(e, ri) != ri or mul(ri, e) != ri:
                raise ValueError('The rank-two projector formula failed.')
            if mul(e, suffix) == mul(suffix, e):
                candidates[i] = e
            else:
                excluded.append(i)
        elif w in (0, u) and z in (0, v):
            ambiguous.append(i)
        else:
            excluded.append(i)
        suffix = plus(suffix, ri)
    if suffix != gram:
        raise ValueError('The query tangent sum failed.')
    return gram, candidates, ambiguous, excluded


def clique(candidates):
    keys = list(candidates)
    adj = [0] * len(keys)
    for i, j in product(range(len(keys)), repeat=2):
        if i != j and mul(candidates[keys[i]], candidates[keys[j]]) == mul(candidates[keys[j]], candidates[keys[i]]):
            adj[i] |= 1 << j
    best = []
    def visit(chosen, available):
        nonlocal best
        if len(chosen) + available.bit_count() <= len(best):
            return
        if not available:
            best = chosen
            return
        while available:
            bit = available & -available
            available ^= bit
            i = bit.bit_length() - 1
            visit(chosen + [i], available & adj[i])
            if len(chosen) + available.bit_count() <= len(best):
                return
    visit([], (1 << len(keys)) - 1)
    return [keys[i] for i in best]


def run(length, trials, seed):
    rng = Random(seed)
    seen = set()
    counts = Counter()
    best = None
    for trial in range(trials):
        word = tuple(rng.randrange(6) for _ in range(length))
        if word in seen or not valid(word):
            continue
        seen.add(word)
        _, candidates, ambiguous, excluded = data(word)
        fixed = clique(candidates)
        upper = len(fixed) + len(ambiguous)
        counts[upper] += 1
        row = dict(word=word, upper=upper, fixed_clique=fixed, ambiguous=ambiguous, excluded=excluded)
        if best is None or upper > best['upper']:
            best = row
    return dict(length=length, trials=trials, seed=seed, valid_distinct=len(seen), upper_histogram=dict(sorted(counts.items())), best=best)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=10000)
    ap.add_argument('--length', type=int, default=18)
    ap.add_argument('--seed', type=int, default=20261023)
    args = ap.parse_args()
    print(json.dumps(run(args.length, args.trials, args.seed), indent=2))
