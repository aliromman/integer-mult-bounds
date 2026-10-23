"""Construct the full sharp Schur normal form over a prime field.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
The written theorem applies over every field.
"""


def need(ok, message):
    if not ok:
        raise ValueError(message)


def columns(a):
    return [list(x) for x in zip(*a)]


def matrix(cols, rows):
    return [list(x) for x in zip(*cols)] if cols else [[] for _ in range(rows)]


def image(a, f, p):
    return [columns(a)[j] for j in f.rr(a, p)[1]]


def basis_extend(cols, candidates, f, p):
    out = [v[:] for v in cols]
    for v in candidates:
        if f.rank(matrix(out + [v], len(v)), p) > len(out):
            out.append(v[:])
    return out


def solve(a, v, f, p):
    r = len(a[0])
    reduced, pivots = f.rr([row + [b] for row, b in zip(a, v)], p)
    need(r not in pivots, 'consistent factor lift')
    out = [0] * r
    for row, col in enumerate(pivots):
        out[col] = reduced[row][-1]
    return out


def complement(a, f, p):
    d = len(a)
    cols = image(a, f, p)
    return basis_extend(cols, f.eye(d), f, p)[len(cols):]


def full_factor_basis(x, y, lifts, f, p):
    d, r = len(x), len(y)
    need(f.rank(x, p) == d and f.rank(y, p) == d, 'sharp factor ranks')
    cols = columns(y) + [solve(x, v, f, p) for v in lifts]
    cols = basis_extend(cols, f.ker(x, p), f, p)
    need(len(cols) == r, 'complete factor basis')
    return matrix(cols, r)


def construct(h, factors, f, p):
    d = len(h)
    sizes = [len(y) for x, y in factors]
    w = f.add(f.eye(d), h, p, -1)
    for (x, y), target in zip(factors, ([[(-v) % p for v in row] for row in w], w, w, h)):
        need(all(sum(x[i][k] * y[k][j] for k in range(len(y))) % p == target[i][j] % p
                 for i in range(d) for j in range(d)), 'sector product')
    wd = f.power(w, d, p)
    star, one = image(wd, f, p), f.ker(wd, p)
    ds, du = len(star), len(one)
    q = matrix(star + one, d)
    qi = f.inv(q, p)
    hh = f.mm(f.mm(qi, h, p), q, p)
    ww = f.add(f.eye(d), hh, p, -1)
    need(all(hh[i][j] == 0 for i in range(d) for j in range(d) if (i < ds) != (j < ds)), 'Fitting block split')
    hs = [row[:ds] for row in hh[:ds]]
    ws = [row[:ds] for row in ww[:ds]]
    hu = [row[ds:] for row in hh[ds:]]
    wu = [row[ds:] for row in ww[ds:]]
    n0 = complement(hs, f, p) if ds else []
    n1 = complement(wu, f, p) if du else []
    a0, a1 = len(n0), len(n1)
    s00 = sizes[0] - ds
    pads = [sizes[1] - d - a1, sizes[2] - d - a1, sizes[3] - d - a0]
    need(sizes[0] >= d - a1 and min(pads) >= 0, 'sharp sector capacities')
    changed = [(f.mm(qi, x, p), f.mm(y, q, p)) if y else ([[] for _ in range(d)], []) for x, y in factors]
    x0, y0 = changed[0]
    b0cols = columns(y0)[:ds] if y0 else []
    b0cols += f.ker(x0[:ds], p) if ds else f.eye(sizes[0])
    b0 = matrix(b0cols, sizes[0])
    need(len(b0cols) == sizes[0], 'sector00 split dimension')
    bases = [b0]
    for (x, y), lifts in zip(changed[1:], ([([0] * ds + v) for v in n1], [([0] * ds + v) for v in n1], [(v + [0] * du) for v in n0])):
        bases.append(full_factor_basis(x, y, lifts, f, p))
    x00 = f.mm(x0, b0, p) if sizes[0] else [[] for _ in range(d)]
    x1 = [row[ds:] for row in x00[ds:]]
    wis = f.inv(ws, p) if ds else []
    hui = f.inv(hu, p) if du else []
    cc = f.inv(f.add(f.eye(du), wu, p), p) if du else []
    n0m, n1m = matrix(n0, ds), matrix(n1, du)
    core0, core1 = 3 * ds + a0, 2 * (du + a1)
    m = sum(sizes) - d

    def plus(*vs):
        return [sum(x) % p for x in zip(*vs)]

    def neg(v):
        return [-x % p for x in v]

    def mv(a, v):
        return [sum(x * y for x, y in zip(row, v)) % p for row in a]

    def split(v):
        at = 0
        def take(k):
            nonlocal at
            out = v[at:at + k]
            at += k
            return out
        result = [take(k) for k in (ds, ds, ds, a0, du, a1, du, a1, s00, *pads)]
        need(at == m, 'full normal dimension')
        return result

    def embed(v):
        xs, ys, zs, t0, xb, bb, yb, cb, t, j10, j01, j11 = split(v)
        v00s = plus(xs, ys, mv(wis, plus(mv(hs, zs), mv(n0m, t0))))
        xt = mv(x1, t)
        shift = mv(cc, xt)
        xu, yu = plus(xb, neg(shift)), plus(yb, neg(shift))
        zu = neg(mv(hui, plus(xt, mv(wu, xu), mv(n1m, bb), mv(wu, yu), mv(n1m, cb))))
        parts = (v00s + t, xs + xu + bb + j10, ys + yu + cb + j01, zs + zu + t0 + j11)
        return sum((mv(g, part) for g, part in zip(bases, parts)), [])

    def action(v, which):
        xs, ys, zs, t0, xb, bb, yb, cb, t, j10, j01, j11 = split(v)
        nt = mv(n0m, t0)
        chosen = xs if which == 0 else ys
        first = plus(mv(hs, plus(chosen, neg(zs))), neg(nt))
        znew = plus(mv(ws, plus(zs, neg(chosen))), neg(nt))
        zero = [0] * ds
        core_star = (first + zero if which == 0 else zero + first) + znew + t0
        xnew = plus(xb, mv(wu, yb), mv(n1m, cb))
        ynew = plus(yb, mv(wu, xb), mv(n1m, bb))
        core_one = (xnew + bb + [0] * (du + a1) if which == 0 else [0] * (du + a1) + ynew + cb)
        return core_star + core_one + [0] * s00 + (j10 if which == 0 else [0] * len(j10)) + (j01 if which else [0] * len(j01)) + j11

    embedding = matrix([embed(v) for v in f.eye(m)], sum(sizes))
    pair = [matrix([action(v, j) for v in f.eye(m)], m) for j in (0, 1)]
    return dict(pair=pair, embedding=embedding, fitting_basis=q, dimensions=[ds, du, a0, a1], scalar_counts=[s00, *pads])
