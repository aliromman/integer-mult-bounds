"""Replay the saved binary word with scalar polynomial sets.

This interpreter does not import a search engine or use prefix inverses.
It multiplies all full and thin operations in execution order. The domain
is F2[z0,...,z4,e]/(e^2). This is not a native characteristic-zero proof.
"""
from copy import deepcopy
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rank(rows):
    pivots = {}
    for value in rows:
        while value:
            low = value & -value
            if low not in pivots:
                pivots[low] = value
                break
            value ^= pivots[low]
    return len(pivots)


def apply(rows, vector):
    return sum(((row & vector).bit_count() % 2) << i for i, row in enumerate(rows))


def query_from_image(projector, packed):
    n = len(projector)
    columns = [apply(projector, 1 << j) for j in range(n)]
    basis = []
    for column in columns:
        if rank(basis+[column]) > len(basis):
            basis.append(column)
    d = len(basis)
    lookup = {0: 0}
    for i, column in enumerate(basis):
        lookup.update({key ^ column: value | (1 << i) for key, value in list(lookup.items())})
    middle = [(packed >> (i*d)) & ((1 << d)-1) for i in range(d)]
    result = []
    for column in columns:
        coordinates = apply(middle, lookup[column])
        value = 0
        for i, b in enumerate(basis):
            if coordinates >> i & 1:
                value ^= b
        result.append(value)
    return [sum(((result[j] >> i) & 1) << j for j in range(n)) for i in range(n)]


def replay(record, witness):
    g, d = len(record['projections']), len(record['projections'][0])
    n = (g+1)*d
    zero = (0,)*(g+1)
    state = [[{zero} if i == j else set() for j in range(n)] for i in range(n)]
    def linear(rows):
        result = []
        for row in rows:
            result.append([set() for _ in range(n)])
            for k in range(n):
                if row >> k & 1:
                    for j in range(n):
                        result[-1][j] ^= state[k][j]
        return result
    def shift(cell, axis):
        return {mon[:axis]+(mon[axis]+1,)+mon[axis+1:] for mon in cell
                if axis != g or mon[g] == 0}
    def query(rows):
        delta = linear(rows)
        for i in range(n):
            for j in range(n):
                state[i][j] ^= shift(delta[i][j], g)
    order = [(j, 0) for j in range(g)]+[(j, 1) for j in reversed(range(g))]
    thin = dict(zip(witness['gaps'], witness['thin']))
    for t in range(len(order)+1):
        if t in thin:
            query(thin[t])
        if t == len(order):
            break
        leaf, half = order[t]
        a = [row ^ ((1 << i) if half else 0) for i, row in enumerate(record['projections'][leaf])]
        p = [0]*n
        for offset in (0, (leaf+1)*d):
            for i, row in enumerate(a):
                p[offset+i] = row << offset
        q = query_from_image(p, witness['full_queries'][t])
        if rank(p) != d or rank(q) != d//2:
            raise ValueError('Wrong full map rank')
        for j in range(n):
            v = apply(q, 1 << j)
            if apply(q, v) or apply(p, v) != v or apply(q, apply(p, 1 << j)) != v:
                raise ValueError('Incompatible full query')
        query(q)
        delta = linear(p)
        for i in range(n):
            for j in range(n):
                state[i][j] ^= delta[i][j] ^ shift(delta[i][j], leaf)
    for i in range(n):
        exponents = (1,)*g if i < d else tuple(int(k == i//d-1) for k in range(g))
        for j in range(n):
            target = {exponents+(0,), exponents+(1,)} if i == j else set()
            if state[i][j] != target:
                return False
    return True


def run():
    data = json.loads((HERE/'variable-query-g5-d4-two-maps.json').read_text())
    record = data['trials'][0]
    witness = next(w for w in record['systems'] if w['status'] == 'sat')
    if not replay(record, witness):
        raise ValueError('The complete literal endpoint fails')
    ranks = [rank(rows) for rows in witness['thin']]
    if ranks != [5, 22]:
        raise ValueError('The thin ranks differ')
    for index in range(2):
        changed = deepcopy(witness)
        changed['thin'][index] = [0]*24
        if replay(record, changed):
            raise ValueError('An omitted thin map was accepted')
    changed = deepcopy(witness)
    changed['thin'][0][0] ^= 1
    if replay(record, changed):
        raise ValueError('A corrupted thin coefficient was accepted')
    margin = 208*41+1040*9-10*1720-52*sum(ranks)
    if margin != -716:
        raise ValueError('Wrong complete linear price')
    print(json.dumps({'status': 'PASS_LITERAL_BINARY_TWO_MAP_WORD',
                      'matrix_entries': 24*24, 'full_calls': 10,
                      'thin_ranks': ranks, 'mutation_controls': 3,
                      'optimistic_linear_margin': margin,
                      'scope': 'Exact binary polynomial endpoint; no native lift or kappa'}))


if __name__ == '__main__':
    run()
