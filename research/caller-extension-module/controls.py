"""Replay the module query and an exact five-call stage control."""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import importlib.util

spec = importlib.util.spec_from_file_location('extension_discovery', Path(__file__).with_name('discovery.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def product(a, b):
    out = [[F(0) for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for k in range(len(b)):
            for j in range(len(b[0])):
                out[i][j] += a[i][k]*b[k][j]
    return out


def literal(ps, ts):
    n = len(ps[0])
    word = [[[F(i == j) for j in range(n)] for i in range(n)]]
    for p, t in zip(ps, ts):
        old = [product(t, a) for a in word]
        word = [old[0]]
        for degree in range(1, len(old)+1):
            extra = product(p, old[degree-1])
            if degree < len(old):
                extra = [[x+y for x, y in zip(a, b)] for a, b in zip(extra, old[degree])]
            word.append(extra)
    return word


def verify():
    identity = [[F(i == j) for j in range(6)] for i in range(6)]
    outcomes = []
    for b, t in ((F(1), F(1)), (F(1), F(0)), (F(2), F(1))):
        pairs, queries, target = m.frame(b, t)
        word = identity
        for (u, v), query in zip(pairs, queries):
            if sum(x*y for x, y in zip(u, v)):
                raise ValueError('The literal query is not square zero')
            word = product(query, word)
        if word != target:
            raise ValueError('The literal complete query product fails')
        result = m.search(b, t, 90)
        if result['status'] != 'EXACT_FIRST_COEFFICIENT_CONTRADICTION':
            raise ValueError('A stated linear failure did not replay')
        outcomes.append({'b': str(b), 't': str(t), 'status': result['status']})
    ps, ts = [], []
    for out, source in ((0, 4), (1, 5)):
        projector = [[F(i == j and i in (out, source)) for j in range(6)] for i in range(6)]
        query = [row[:] for row in identity]
        query[out][source] = F(-2)
        ps.append(projector)
        ts.append(query)
    us, vs, _ = m.e.canonical(2)
    for u, v in zip(us, vs):
        ps.append([[F(i == j and i in (2, 3)) for j in range(6)] for i in range(6)])
        query = [row[:] for row in identity]
        for i in range(2):
            for j in range(2):
                query[i+2][j+2] -= 2*u[i]*v[j]
        ts.append(query)
    target = [row[:] for row in identity]
    target[0][4] = target[1][5] = F(-2)
    target[2][2] = target[3][3] = F(-1)
    for minus in (False, True):
        word = literal(ps, ts if minus else [identity]*len(ps))
        for degree, coefficient in enumerate(word):
            diagonal = [[F(int(i == j)*(comb(3 if i in (2, 3) else 1, degree)
                         if degree <= (3 if i in (2, 3) else 1) else 0))
                         for j in range(6)] for i in range(6)]
            expected = product(diagonal, target) if minus else diagonal
            if coefficient != expected:
                raise ValueError('The complete five-call stage control fails')
    changed = ts[:]
    changed[-1] = identity
    if literal(ps, changed) == literal(ps, ts):
        raise ValueError('The omitted-query control was not rejected')
    result = {'status': 'PASS_EXTENSION_MODULE_CONTROLS', 'query_frames': outcomes,
              'five_call_stage_control': 'PASS_ALL_COEFFICIENTS',
              'omitted_query_rejected': True, 'new_kappa': False}
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    verify()
