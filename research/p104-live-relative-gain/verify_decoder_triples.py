#!/usr/bin/env python3
"""Check the minimum integer movement rank for three frozen BIT returns."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
import gzip
import hashlib
from itertools import combinations
import json
from math import gcd, isqrt, lcm
from operator import or_
from pathlib import Path

WORD_SHA = 'd28f0d28d9d87d5fe57b8ef19e12faf98014f0bbef781987ce22daaa82f21cb9'
PORTS = 111
LIVE = (1 << PORTS) - 1


def require(condition, message):
    if not condition:
        raise ValueError(message)


def indices(row):
    while row:
        bit = row & -row
        yield bit.bit_length() - 1
        row ^= bit


def reduce_binary(row, basis):
    while row and (pivot := row.bit_length() - 1) in basis:
        row ^= basis[pivot]
    return row


def binary_basis(rows):
    basis = {}
    for row in rows:
        remainder = reduce_binary(row, basis)
        if remainder:
            basis[remainder.bit_length() - 1] = remainder
    return basis


def run_tape(word, reverse):
    start, finish = (word[key] for key in (('final', 'start') if reverse else ('start', 'final')))
    charts = list(start) + [None] * word['owned_fields']
    values = [1 << i for i in range(len(start))] + [0] * word['owned_fields']
    paid = {}
    program = word['program']
    times = reversed(range(len(program))) if reverse else range(len(program))
    for time in times:
        op, target, *args = program[time]
        if op == 'ADD':
            source, chart = args
            require(charts[source] == charts[target] == chart, 'ADD chart')
            if reverse:
                target, source = source, target
            values[target] ^= values[source]
        elif op == 'SWAP':
            source, = args
            values[target], values[source] = values[source], values[target]
            charts[target], charts[source] = charts[source], charts[target]
        elif op == ('CLEAR' if reverse else 'INIT'):
            require(charts[target] is None, 'a private cell must be free before use')
            charts[target], = args
            values[target] = 0
        elif op == ('INIT' if reverse else 'CLEAR'):
            require(charts[target] == args[0], 'private cell retirement chart')
            charts[target] = None
            values[target] = 0
        else:
            require(op == 'MOVE', 'unknown instruction')
            source, destination = reversed(args) if reverse else args
            require(charts[target] == source, 'MOVE source chart')
            paid[time] = (source, destination, values[target])
            charts[target] = destination
    require(charts == list(finish) + [None] * word['owned_fields'], 'final chart inventory')
    for i in range(len(start)):
        expected = 1 << i
        if i < PORTS and reverse:
            expected ^= 1 << (i + PORTS)
        if PORTS <= i < 2 * PORTS and not reverse:
            expected ^= 1 << (i - PORTS)
        require(values[i] == expected, 'public row')
    return paid


def write_integral(row):
    fractions = [Fraction(x) for x in row]
    denominator = lcm(*(x.denominator for x in fractions))
    integers = [int(x * denominator) for x in fractions]
    divisor = gcd(*integers)
    require(divisor > 0, 'nonzero chart basis row')
    return [v // divisor for v in integers]


def residue(row):
    return sum((v & 1) << i for i, v in enumerate(row))


def enlarge(rows, parity):
    pivot = (parity & -parity).bit_length() - 1
    return [[2 * row[pivot] if i == pivot else value - ((parity >> i) & 1) * row[pivot]
             for i, value in enumerate(row)] for row in rows]


def lattice_steps(roots):
    scaled = []
    for row in roots:
        norm = sum(x * x for x in row)
        length = isqrt(norm)
        require(length * length == norm and all(2 * x % length == 0 for x in row), 'integral scaled root')
        scaled.append([2 * x // length for x in row])
    steps = []
    for i in range(len(scaled)):
        parity = residue(scaled[i])
        if parity:
            steps.append(parity)
            scaled = enlarge(scaled, parity)
    return steps


def saturated_image(rows, steps):
    for parity in steps:
        rows = enlarge(rows, parity)
    rows = [write_integral(row) for row in rows]
    while True:
        pivots = {}
        for i, row in enumerate(rows):
            value, witness = residue(row), 1 << i
            while value and (j := value.bit_length() - 1) in pivots:
                old, mask = pivots[j]
                value ^= old
                witness ^= mask
            if value:
                pivots[value.bit_length() - 1] = (value, witness)
                continue
            terms = list(indices(witness))
            total = [sum(rows[j][k] for j in terms) for k in range(8)]
            require(any(total) and all(x % 2 == 0 for x in total), 'integral saturation relation')
            rows[i] = write_integral(total)
            break
        else:
            return [residue(row) for row in rows]


def collect(word, frames, sources, adjoints):
    returns = [time for time, (_, b, _) in sources.items() if b == 0]
    selectors = {adjoints[time][2] for time in returns}
    require(len(returns) == len(selectors) == 9, 'nine distinct return groups')
    require(all(x and x & ~(LIVE << PORTS) == 0 for x in selectors), 'only late public outputs use the returns')
    removed = {time for time in sources if adjoints[time][2] in selectors}
    free = sum(1 << i for i, chart in enumerate(word['start']) if chart == 0)
    targets = [sources[time][2] & ~free for time in returns]
    require(all(not (x & ~LIVE) for x in targets) and len(binary_basis(targets)) == 9, 'return input support and rank')
    stock = defaultdict(list)
    for i, chart in enumerate(word['start']):
        stock[chart].append((1 << i) & LIVE)
    for time, (_, chart, value) in sources.items():
        if time not in removed and not (value & (LIVE << PORTS)):
            stock[chart].append(value & LIVE)
    require(len(removed) == 27 and len(sources) == 2874, 'paid transfer inventory')
    for time, target in zip(returns, targets):
        group = [sources[t][:2] for t in sorted(removed) if adjoints[t][2] == adjoints[time][2]]
        require(len(group) == 3 and group[-1][1] == 0, 'three movements to the common chart')
        require(all(x[1] == y[0] for x, y in zip(group, group[1:])), 'return chain')
        for a, b in group:
            require(all(sum(x * y for x, y in zip(u, v)) == 0 for u in frames[b]['B'] for v in frames[a]['A']),
                    'nested movement')
        require(sum(frames[a]['dim'] - frames[b]['dim'] for a, b in group) == 7, 'return chain cost')
        sensors = [v for a, _ in group for v in stock[a]]
        require(reduce_binary(target, binary_basis(sensors)) == 0, 'constructive rank-seven source path')
    return targets, stock


def verify(path):
    blob = path.read_bytes()
    require(hashlib.sha256(blob).hexdigest() == WORD_SHA, 'pinned #7786 compressed word')
    word = json.loads(gzip.decompress(blob))
    frames = {int(k): v for k, v in word['frames'].items()}
    sources, adjoints = run_tape(word, False), run_tape(word, True)
    targets, stock = collect(word, frames, sources, adjoints)
    roots = [write_integral(frames[c]['B'][0]) for c in word['start'][:PORTS]]
    steps = lattice_steps(roots)
    reduced = {c: saturated_image([write_integral(v) for v in f['B']], steps) for c, f in frames.items()}
    require(all(len(binary_basis(reduced[c])) == f['dim'] for c, f in frames.items()), 'reduced chart dimension')
    lines = [reduced[c][0] for c in word['start'][:PORTS]]
    require(len(set(lines)) == PORTS, 'distinct reduced input lines')
    support = {c: reduce(or_, rows, 0) for c, rows in stock.items()}
    for c, fields in support.items():
        for i in indices(fields):
            require(all(sum(x * y for x, y in zip(roots[i], v)) == 0 for v in frames[c]['A']), 'stock root incidence')
    small = [bin(v).count('1') for c, v in support.items() if frames[c]['dim'] <= 2]
    require(max(small) == 2, 'small-bank support limit')
    require(frames[0]['dim'] == 0 and support.get(0, 0) == 0, 'empty source has no live root')
    hist, uncovered = Counter(), []
    for triple in combinations(range(9), 3):
        goals = [targets[i] for i in triple]
        labels = [sum(((g >> i) & 1) << j for j, g in enumerate(goals)) for i in range(PORTS)]
        def tensor(z, u):
            return sum(z << (8 * j) for j in indices(u))
        required = binary_basis(tensor(line, u) for line, u in zip(lines, labels))
        hist[len(required)] += 1
        require(len(required) >= 20, 'tensor lower bound')
        if len(required) > 20:
            continue
        admissible = {u: {z for z in range(256) if not reduce_binary(tensor(z, u), required)} for u in range(1, 8)}
        large = [c for c in support if frames[c]['dim'] > 2 and any(all(z in space for z in reduced[c]) for space in admissible.values())]
        require(all(frames[c]['dim'] == 4 for c in large), 'no usable stocked chart above grade four')
        remaining = reduce(or_, goals) & ~reduce(or_, (support[c] for c in large), 0)
        uncovered.append(bin(remaining).count('1'))
        require(bin(remaining).count('1') > 20 * max(small), 'insufficient small receiving banks')
    require(hist == {20: 83, 21: 1} and (min(uncovered), max(uncovered)) == (58, 103), 'all 84 triples')
    return {'word_sha256': WORD_SHA, 'stock_charts': len(support), 'old_return_moves': 27,
            'tensor_rank_histogram': dict(hist), 'min_uncovered_after_free_large_banks': min(uncovered),
            'max_uncovered_after_free_large_banks': max(uncovered), 'small_bank_support_cap': max(small),
            'sharp_integer_rank_of_each_triple': 21}


if __name__ == '__main__':
    result = verify(Path(__file__).with_name('BIT-WORD.json.gz'))
    print(json.dumps(result, indent=2))
