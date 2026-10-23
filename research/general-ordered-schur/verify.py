"""Compare the normal form with literal full-header Schur complements.

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import sys
import types

P = Path(__file__).resolve().parent
ROOT = P.parents[1]


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load(name, path):
    m = types.ModuleType(name)
    m.__file__ = str(path)
    sys.modules[name] = m
    exec(compile(path.read_bytes(), str(path), 'exec'), m.__dict__)
    return m


def rank_factor(a, f, p):
    d = len(a)
    piv = f.rr(a, p)[1]
    x = [[row[j] for j in piv] for row in a]
    if not piv:
        return x, []
    r, _ = f.rr([u + v for u, v in zip(x, a)], p)
    return x, [row[len(piv):] for row in r[:len(piv)]]


def factors(h, f, p, variant, padding):
    d = len(h)
    w = f.add(f.eye(d), h, p, -1)
    x, y = rank_factor([[-v % p for v in row] for row in w], f, p)
    if variant:
        x = [row + [int(variant in (1, 3) and i == 0), 0] for i, row in enumerate(x)]
        y += [[0] * d, [int(variant in (2, 3) and j == d - 1) for j in range(d)]]
    out = [(x, y)]
    for a, pad in zip((w, w, h), padding):
        nn = f.complement(a, p) or [[] for _ in range(d)]
        nullity = d - f.rank(a, p)
        out.append(([row + nn[i] + [0] * pad for i, row in enumerate(a)],
                    f.eye(d) + [[0] * d for _ in range(nullity + pad)]))
    return out


def inspect(h, f, normal, p, variant=0, padding=(0, 0, 0), twist=False):
    fs = factors(h, f, p, variant, padding)
    if twist:
        rng = random.Random(1000 * p + len(h) + variant)
        changed = []
        for x, y in fs:
            r = len(y)
            if not r:
                changed.append((x, y))
                continue
            g = f.eye(r)
            for _ in range(3 * r):
                if r > 1:
                    i, j = rng.sample(range(r), 2)
                    c = rng.randrange(1, p)
                    g[i] = [(a + c * b) % p for a, b in zip(g[i], g[j])]
            changed.append((f.mm(x, f.inv(g, p), p), f.mm(g, y, p)))
        fs = changed
    sizes = [len(y) for x, y in fs]
    d, n = len(h), sum(sizes)
    x = [sum((a[i] for a, b in fs), []) for i in range(d)]
    y = sum((b for a, b in fs), [])
    need(f.mm(x, y, p) == f.eye(d), 'actual split')
    kb = f.tr(f.ker(x, p))
    header = [a + b for a, b in zip(y, kb)]
    hi = f.inv(header, p)
    literal = []
    for labels in ((0, 1, 0, 1), (0, 0, 1, 1)):
        diagonal = [v for size, v in zip(sizes, labels) for _ in range(size)]
        ambient = [[int(i == j) * diagonal[i] for j in range(n)] for i in range(n)]
        a = f.mm(f.mm(hi, ambient, p), header, p)
        need([row[:d] for row in a[:d]] == f.eye(d), 'identity pivot')
        literal.append(f.add([row[d:] for row in a[d:]],
                             f.mm([row[:d] for row in a[d:]], [row[d:] for row in a[:d]], p), p, -1))
    normal_form = normal.construct(h, fs, f, p)
    embedding = normal_form['embedding']
    need(f.mm(x, embedding, p) == [[0] * (n - d) for _ in range(d)], 'explicit map lies in kernel')
    link = f.mm(hi, embedding, p)[d:]
    need(f.rank(link, p) == n - d, 'explicit residual map invertible')
    for a, b in zip(literal, normal_form['pair']):
        need(f.mm(a, link, p) == f.mm(link, b, p), 'whole ordered intertwiner')
        need(f.mm(a, a, p) == a, 'idempotent residual')
    comm = f.add(f.mm(literal[0], literal[1], p), f.mm(literal[1], literal[0], p), p, -1)
    need(f.rank(comm, p) == 2 * d, 'actual sharpness')
    return dict(p=p, h=h, dimensions=normal_form['dimensions'], scalars=normal_form['scalar_counts'], residual_dimension=n-d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write-result', action='store_true')
    args = ap.parse_args()
    for row in json.loads((P / 'SOURCES.json').read_text())['files']:
        if row['usage'] == 'input':
            need(hashlib.sha256((ROOT / row['path']).read_bytes()).hexdigest() == row['sha256'], 'source hash')
    f = load('retained_field', ROOT / 'research/unipotent-ordered-schur/verify.py')
    normal = load('normal', P / 'normal.py')
    digest, count, transported = hashlib.sha256(), 0, 0
    domains = []
    for p, d in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (5, 1), (5, 2)):
        local = 0
        for data in product(range(p), repeat=d*d):
            h = [list(data[i*d:(i+1)*d]) for i in range(d)]
            result = inspect(h, f, normal, p)
            digest.update(json.dumps(result, sort_keys=True).encode())
            count += 1
            local += 1
            if count % 41 == 0:
                for variant in (1, 2, 3):
                    inspect(h, f, normal, p, variant, (1, 2, 1), True)
                    transported += 1
        domains.append([p, d, local])
    mixed = []
    for p in (2, 3, 5, 7):
        for h in ([[0,0,0],[0,1,1],[0,0,1]],
                  [[1,1,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,1]],
                  [[0,1,0,0,0],[0,0,0,0,0],[0,0,1,1,0],[0,0,0,1,0],[0,0,0,0,2%p]]):
            mixed.append(inspect(h, f, normal, p, 3, (1, 1, 1), True))
    controls = []
    # Sector00 ranks differ, while the complete ordered representative is fixed.
    h = [[0,0],[0,1]]
    variants = []
    for v in (1, 2, 3):
        fs = factors(h, f, 3, v, (0, 0, 0))
        out = normal.construct(h, fs, f, 3)
        variants.append(out['pair'])
    need(variants[0] == variants[1] == variants[2], 'sector00 data disappear across mixed components')
    controls.append('inequivalent sector00 factor systems')
    fs = factors([[1]], f, 3, 1, (0, 0, 0))
    fs[0][0][0][0] = 2
    fs[0][1][0][0] = 1
    try:
        normal.construct([[1]], fs, f, 3)
    except ValueError:
        controls.append('wrong sector product')
    else:
        raise ValueError('wrong product accepted')
    # H=2 is legal: I+W is singular there and the one-primary formula cannot be used globally.
    boundary = inspect([[2]], f, normal, 3)
    need(boundary['dimensions'] == [1,0,0,0], 'H2 belongs to invertible-W component')
    controls.append('singular 2I-H outside one-primary component')
    result = dict(status='PASS_GENERAL_ORDERED_SCHUR', domains=domains, complete_H=count,
                  transported=transported, mixed=mixed, controls=controls, digest=digest.hexdigest(),
                  normal_sha256=hashlib.sha256((P/'normal.py').read_bytes()).hexdigest())
    if args.write_result:
        (P/'RESULT.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        need(json.dumps(result,sort_keys=True)==json.dumps(json.loads((P/'RESULT.json').read_text()),sort_keys=True),'typed complete receipt')
    print(json.dumps(dict(status=result['status'],complete_H=count,transported=transported,mixed=len(mixed),controls=len(controls))))


if __name__ == '__main__':
    main()
