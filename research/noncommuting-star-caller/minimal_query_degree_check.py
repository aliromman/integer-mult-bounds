"""Check exact controls for the minimum-query degree theorem."""
import argparse
from collections import Counter
from decimal import Decimal, localcontext, ROUND_FLOOR
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import random


def add(a, b):
    out = a.copy()
    for d, v in b.items():
        out[d] = out.get(d, 0)+v
        if not out[d]:
            del out[d]
    return out


def scale(a, c):
    return {d: c*v for d, v in a.items() if c*v}


def mul(a, b):
    out = {}
    for (z, t), x in a.items():
        for (zz, tt), y in b.items():
            out = add(out, {(z+zz, t+tt): x*y})
    return out


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def pmm(a, b):
    out = [[{} for _ in b[0]] for _ in a]
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            for j, y in enumerate(b[k]):
                out[i][j] = add(out[i][j], mul(x, y))
    return out


def literal(a):
    return [[{(0, 0): x} if x else {} for x in row] for row in a]


def det(a):
    out, n = {}, len(a)
    for perm in itertools.permutations(range(n)):
        term = {(0, 0): (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))}
        for i, j in enumerate(perm):
            term = mul(term, a[i][j])
        out = add(out, term)
    return out


def stage(p, q):
    return [[add(add({(0, 0): int(i == j)-p[i][j]}, {(1, 0): p[i][j]}),
                 {(1, 1): q[i][j]}) for j in range(len(p))] for i in range(len(p))]


def run():
    rng = random.Random(5840)
    transfer_cases = 0
    for n in range(1, 6):
        for _ in range(4):
            k = [[rng.randrange(-2, 3) if j < i else 0 for j in range(n)] for i in range(n)]
            u = [[rng.randrange(-2, 3) if j < i else int(i == j) for j in range(n)] for i in range(n)]
            v = [[rng.randrange(-2, 3) if j > i else (i+1)*int(i == j) for j in range(n)] for i in range(n)]
            powers = [eye(n)]
            for j in range(n):
                powers.append(mm(powers[-1], k))
            if any(any(row) for row in powers[-1]):
                raise ValueError('The chronological matrix is not nilpotent')
            h = [[sum(powers[j][a][b] for j in range(0, n, 2)) for b in range(n)] for a in range(n)]
            target = mm(mm(u, h), v)
            d = [[-x for x in row] for row in mm(mm(mm(u, k), h), v)]
            value = literal(d)
            for j in range(n):
                coefficient = mm(mm(u, powers[j]), v)
                value = [[add(value[a][b], {(0, j+1): coefficient[a][b]}) for b in range(n)] for a in range(n)]
            for t in (-1, 1):
                actual = [[sum(c*t**degree[1] for degree, c in x.items()) for x in row] for row in value]
                if actual != [[t*x for x in row] for row in target]:
                    raise ValueError('A query-character endpoint fails')
            wanted = {(z, t+n): c for (z, t), c in det(literal(target)).items()}
            if det(value) != wanted or det(literal(d)):
                raise ValueError('The full determinant identity fails')
            transfer_cases += 1

    excess_controls = []
    for layers, block_rank, excess in ((3, 2, 1), (2, 3, 1), (3, 2, 2)):
        total, n = layers*block_rank, layers*block_rank-excess
        k = [[rng.randrange(-2, 3) if j//block_rank < i//block_rank else 0
              for j in range(total)] for i in range(total)]
        powers = [eye(total)]
        for _ in range(layers):
            powers.append(mm(powers[-1], k))
        if any(any(row) for row in powers[-1]):
            raise ValueError('A query-layer matrix has the wrong nilpotency index')
        h = [[sum(powers[j][a][b] for j in range(0, layers, 2)) for b in range(total)] for a in range(total)]
        for _ in range(20):
            u = [[rng.randrange(-2, 3) for _ in range(total)] for _ in range(n)]
            v = [[rng.randrange(-2, 3) for _ in range(n)] for _ in range(total)]
            target = mm(mm(u, h), v)
            if det(literal(target)):
                break
        else:
            raise ValueError('No invertible rectangular transfer control was found')
        d = [[-x for x in row] for row in mm(mm(mm(u, k), h), v)]
        value = literal(d)
        for j in range(layers):
            coefficient = mm(mm(u, powers[j]), v)
            value = [[add(value[a][b], {(0, j+1): coefficient[a][b]}) for b in range(n)] for a in range(n)]
        polynomial = det(value)
        lower = max(0, total-layers*excess)
        if not polynomial or min(t for _, t in polynomial) < lower:
            raise ValueError('The extra-query divisibility bound fails')
        for t in (-1, 1):
            actual = [[sum(c*t**degree[1] for degree, c in x.items()) for x in row] for row in value]
            if actual != [[t*x for x in row] for row in target]:
                raise ValueError('A rectangular query-character endpoint fails')
        excess_controls.append({'layers': layers, 'query_rank_per_layer': block_rank,
                                'output_count': n, 'excess_query_rank': excess,
                                'proved_minimum_t_order': lower,
                                'actual_t_order': min(t for _, t in polynomial)})

    minor_cases = 0
    for m, n in ((4, 2), (5, 3), (5, 4)):
        for _ in range(3):
            a = [[rng.randrange(-2, 3) for _ in range(2)] for _ in range(m-2)]
            b = [[rng.randrange(-2, 3) for _ in range(m-2)] for _ in range(2)]
            ba = mm(b, a)
            x = eye(2)+a
            y = [[int(i == j)-ba[i][j] for j in range(2)]+b[i] for i in range(2)]
            p = mm(x, y)
            q = mm(mm(x, [[0, 1], [0, 0]]), y)
            f = stage(p, q)
            for rows in itertools.combinations(range(m), n):
                for cols in itertools.combinations(range(m), n):
                    minor = det([[f[i][j] for j in cols] for i in rows])
                    if any(z+t > 2 for z, t in minor):
                        raise ValueError('An exterior coefficient exceeds total degree two')
                    minor_cases += 1

    block_minor_cases = 0
    m, n, rank = 6, 3, 4
    a = [[rng.randrange(-2, 3) for _ in range(rank)] for _ in range(m-rank)]
    b = [[rng.randrange(-2, 3) for _ in range(m-rank)] for _ in range(rank)]
    ba = mm(b, a)
    x = eye(rank)+a
    y = [[int(i == j)-ba[i][j] for j in range(rank)]+b[i] for i in range(rank)]
    nilpotent = [[0]*rank for _ in range(rank)]
    nilpotent[0][1] = nilpotent[2][3] = 1
    f = stage(mm(x, y), mm(mm(x, nilpotent), y))
    for rows in itertools.combinations(range(m), n):
        for cols in itertools.combinations(range(m), n):
            minor = det([[f[i][j] for j in cols] for i in rows])
            if any(z+t > 4 for z, t in minor):
                raise ValueError('A rank-two query block exceeds total degree four')
            block_minor_cases += 1

    sharp = []
    for n in range(1, 5):
        m, word = 2*n, literal(eye(2*n))
        for i in range(n):
            p, q = [[0]*m for _ in range(m)], [[0]*m for _ in range(m)]
            p[i][i] = p[n+i][n+i] = q[n+i][i] = 1
            word = pmm(stage(p, q), word)
        output = [[word[n+i][j] for j in range(n)] for i in range(n)]
        if output != [[{(1, 1): 1} if i == j else {} for j in range(n)] for i in range(n)]:
            raise ValueError('The initialized-work sharpness control fails')
        sharp.append({'inputs': n, 'work_fields': n, 'calls': n, 'target_degree_sum': n})

    word = literal(eye(2))
    for q in ([[0, 1], [0, 0]], [[0, 0], [-1, 0]], [[0, 1], [0, 0]]):
        word = pmm(stage(eye(2), q), word)
    word = pmm(literal([[0, -1], [1, 0]]), word)
    if det(word) != {(6, 0): 1}:
        raise ValueError('The additional-query triangle control fails')
    for t in (-1, 1):
        evaluated = [[{(z, 0): sum(c*t**tt for (zz, tt), c in a.items() if zz == z)
                       for z in {d[0] for d in a} if sum(c*t**tt for (zz, tt), c in a.items() if zz == z)}
                      for a in row] for row in word]
        if evaluated != [[{(3, 0): t} if i == j else {} for j in range(2)] for i in range(2)]:
            raise ValueError('The triangle endpoint fails')
    invalid = det(stage(eye(2), [[1, 0], [0, 0]]))
    if not any(z+t > 2 for z, t in invalid):
        raise ValueError('The nonnilpotent-query control was accepted')
    root = Path(__file__).resolve().parent
    pins = json.loads((root/'SOURCE-PINS.json').read_text())
    pin = next(p for p in pins if p['local_path'] == 'sources/selected-word.json.gz')
    data = (root/pin['local_path']).read_bytes()
    if hashlib.sha256(data).hexdigest() != pin['sha256']:
        raise ValueError('The retained source bytes differ from their pin')
    source = json.loads(gzip.decompress(data))
    paid = Counter(abs(len(source['frames'][op[3]])-len(source['frames'][op[2]]))
                   for op in source['program'] if op[0] == 'MOVE')
    paid.pop(0, None)
    main = set(source['live']) | {r+112 for r in source['live']}
    helpers = Counter(len(source['frames'][source['final'][r]])-len(source['frames'][source['initial'][r]])
                      for r in range(source['public']) if r not in main)
    expected_paid = {1:2203, 2:459, 3:493, 4:92, 5:168, 6:18, 7:14}
    expected_helpers = {1:2, 2:5, 4:89, 5:11, 7:5, 8:16, 9:412}
    if dict(paid) != expected_paid or dict(helpers) != expected_helpers or len(source['live']) != 104:
        raise ValueError('The complete source profiles differ')
    ranks = sorted(set(paid) | set(helpers))
    net = {r: paid[r]-helpers[r] for r in ranks}
    linear = sum(r*c for r, c in net.items())
    hinges = {j: sum(c*max(r-j, 0) for r, c in net.items()) for j in range(1, max(ranks)+1)}
    if linear != 1720 or any(v > 0 for v in hinges.values()) or 9*104-linear != -784:
        raise ValueError('The complete concave tariff control fails')
    cutoffs = {h: h//(8*h-784)+1 for h in range(99, 105)}
    if cutoffs != {99:13, 100:7, 101:5, 102:4, 103:3, 104:3}:
        raise ValueError('The first possible positive linear cases differ')
    if any(8*h-784 > 0 or 7*h-784 > -98 for h in range(99)):
        raise ValueError('The unbounded full-call count exclusion fails')
    power_controls = []
    for excess, numerator, denominator in ((101, 997, 1000), (102, 4983, 5000)):
        unit, bounds = 10**9, {}
        for r in ranks+[17]:
            with localcontext() as context:
                context.prec = 50
                estimate = (Decimal(r).ln()*numerator/denominator).exp()*unit
                lo = int(estimate.to_integral_value(rounding=ROUND_FLOOR))
            target = pow(r, numerator)*pow(unit, denominator)
            if not pow(lo, denominator) <= target < pow(lo+1, denominator):
                raise ValueError('An integer-power interval does not enclose its target')
            bounds[r] = (lo, lo+1)
        coefficients = {r: -c for r, c in net.items()}
        coefficients[9] = coefficients.get(9, 0)+104-excess
        coefficients[17] = excess
        low = sum(c*bounds[r][0 if c >= 0 else 1] for r, c in coefficients.items())
        high = sum(c*bounds[r][1 if c >= 0 else 0] for r, c in coefficients.items())
        if high >= 0 or any(c > 0 and r < 8 or c < 0 and r >= 8 for r, c in coefficients.items()):
            raise ValueError('A power-tariff exclusion or its sign separation fails')
        power_controls.append({'excess_query_rank': excess, 'tariff_exponent': [numerator, denominator],
                               'interval_scale': unit, 'margin_numerator_interval': [low, high],
                               'power_numerator_bounds': bounds,
                               'excludes_every_full_call_count': True})
    return {'status': 'PASS_MINIMAL_QUERY_DEGREE_CONTROLS', 'new_kappa': False,
            'exact_transfer_cases': transfer_cases, 'exact_exterior_minors': minor_cases,
            'excess_query_controls': excess_controls,
            'rank_two_query_exterior_minors': block_minor_cases,
            'sharp_initialized_controls': sharp, 'additional_query_triangle': 'PASS',
            'nonnilpotent_query_mutant': 'REJECT',
            'source_word_sha256': pin['sha256'], 'source_commit': pin['commit'],
            'source_paid': dict(sorted(paid.items())), 'source_helpers': dict(sorted(helpers.items())),
            'net_linear_mass': linear, 'hinges': hinges, 'complete_margin_per_full_call_upper_bound': -784,
            'strict_exclusion_max_excess_rank': 98, 'first_possible_full_call_counts': cutoffs,
            'power_tariff_exclusions': power_controls,
            'scope': 'Exact integer controls for the written general proof; no new profitable caller'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
