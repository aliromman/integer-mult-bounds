"""Check the cyclic-stage rank budget and selected-source price ceiling.

The exact finite checks support the written proof. They do not search callers.
This program uses the standard library and does not change any file.
"""
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import random
import runpy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rank(a, prime):
    a = [list(row) for row in a]
    pivot = 0
    for col in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][col] % prime), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        inv = pow(a[pivot][col] % prime, -1, prime)
        a[pivot] = [x*inv % prime for x in a[pivot]]
        for i in range(len(a)):
            if i != pivot:
                factor = a[i][col] % prime
                a[i] = [(x-factor*y) % prime for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def multiply(a, b, prime):
    return [[sum(x*y for x, y in zip(row, col)) % prime for col in zip(*b)] for row in a]


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def add(a, b, scalar=1):
    return [[x+scalar*y for x, y in zip(r, s)] for r, s in zip(a, b)]


def roots(period):
    prime = period+1
    while any(prime % d == 0 for d in range(2, isqrt(prime)+1)) or (prime-1) % period:
        prime += 1
    root = next(g for g in range(2, prime) if pow(g, period, prime) == 1
                and all(pow(g, d, prime) != 1 for d in range(1, period)))
    return prime, [pow(root, j, prime) for j in range(period)]


def loewner(left, right, value, prime):
    lhs, rhs = {x: value(x) for x in left}, {y: value(y) for y in right}
    n = len(lhs[left[0]])
    return [[(lhs[x][i][j]-rhs[y][i][j])*pow((x-y) % prime, -1, prime) % prime
             for y in right for j in range(n)] for x in left for i in range(n)]


def scalar_checks():
    cases = []
    for period in range(2, 20):
        prime, points = roots(period)
        for shift in range(min(period, 4)):
            ordered = points[shift:]+points[:shift]
            left, right = ordered[:period//2], ordered[period//2:]
            for degree in range(period):
                actual = rank(loewner(left, right, lambda z: [[pow(z, degree, prime)]], prime), prime)
                require(actual == min(degree, period-degree), 'A scalar cyclic rank differs')
                cases.append((period, prime, shift, degree, actual))
    return cases


def product_checks():
    rng, results = random.Random(290629), []
    for n in (2, 3, 4):
        for length in range(1, 6):
            prime, points = roots(7)
            factors = []
            for _ in range(length):
                width = rng.randint(1, 2)
                a = [[rng.randrange(prime) for _ in range(n)] for _ in range(n)]
                x = [[rng.randrange(prime) for _ in range(width)] for _ in range(n)]
                y = [[rng.randrange(prime) for _ in range(n)] for _ in range(width)]
                factors.append((a, multiply(x, y, prime)))

            def word(z, lower=0, upper=None):
                out = identity(n)
                for a, b in factors[lower:upper]:
                    out = multiply(add(a, b, z), out, prime)
                return out

            left, right = points[:3], points[3:]
            direct = loewner(left, right, word, prime)
            actual, budget = rank(direct, prime), sum(rank(b, prime) for _, b in factors)
            require(actual <= budget, 'An affine word exceeds its slope-rank budget')
            for x in left:
                for y in right:
                    total = [[0]*n for _ in range(n)]
                    for site, (_, b) in enumerate(factors):
                        term = multiply(multiply(word(x, site+1), b, prime), word(y, 0, site), prime)
                        total = add(total, term)
                    expected = add(word(x), word(y), -1)
                    require(all((expected[i][j]-(x-y)*total[i][j]) % prime == 0
                                for i in range(n) for j in range(n)), 'The literal telescoping product differs')
            results.append({'fields': n, 'factors': length, 'rank': actual, 'budget': budget})
    return results


def target_checks():
    records = []
    for period, degrees, calls in ((7, [3]+[4]*5, 8), (7, [3]*5+[4]*3, 10),
                                    (6, [2]*8+[3]*4, 14)):
        require((sum(degrees)-2*calls) % period == 0, 'The determinant-congruence control differs')
        prime, points = roots(period)
        for coupled in (1, 2):
            expanded = [d for d in degrees for _ in range(coupled)]

            def value(z):
                return [[pow(z, d, prime) if i == j else 0 for j in range(len(expanded))]
                        for i, d in enumerate(expanded)]

            actual = rank(loewner(points[:period//2], points[period//2:], value, prime), prime)
            expected = coupled*sum(min(d, period-d) for d in degrees)
            require(actual == expected, 'A complete target rank differs')
            budget = 2*coupled*calls
            require((actual > budget) == (period == 7), 'A target gate changed')
            records.append({'period': period, 'coupled': coupled, 'target_rank': actual,
                            'paid_rank': budget, 'excluded': actual > budget})
    prime, points = roots(7)
    low = rank(loewner(points[:2], points[2:4], lambda z: [[pow(z, 3, prime)]], prime), prime)
    require(low == 2 and low != 3, 'The omitted-character control fails')
    require(min(4, 7-4) != 4, 'The false polynomial-degree control fails')
    return records


def inputs():
    manifest = json.loads((HERE/'SOURCES.json').read_text())
    digests = {}
    for item in manifest['files']:
        if item['usage'] != 'input':
            continue
        raw = (ROOT/item['path']).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        require(digest == item['sha256'], 'A scientific input hash differs')
        require(hashlib.sha256(raw+b' ').hexdigest() != digest, 'The changed-source control fails')
        digests[item['path']] = digest
    return digests


def price_check():
    module = runpy.run_path(str(ROOT/'research/alternating-caller-closure/price_check.py'), run_name='pinned_price_bounds')
    power_bounds, outward = module['power_bounds'], module['outward']
    source = json.loads((ROOT/'research/alternating-caller-closure/sources/selected6030-result.json').read_text())
    coefficients = {int(r): -v for r, v in source['local']['paid'].items()}
    for r, count in source['local']['helpers'].items():
        coefficients[int(r)] = coefficients.get(int(r), 0)+count
    coefficients[17] = coefficients.get(17, 0)+104
    require(sum(coefficients.values()) == -2787, 'The count endpoint differs')
    require(sum(r*v for r, v in coefficients.items()) == 48, 'The linear endpoint differs')
    require(coefficients[1] == -2205 and coefficients[4] == 7, 'The derivative guard differs')
    require(all(c < 0 if r < 6 else c > 0 for r, c in coefficients.items() if r not in (1, 4, 6)),
            'A derivative sign differs')
    records = []
    for value, positive in (('0.00511060055', True), ('0.00511060056', False)):
        prices = {r: power_bounds(r, F(value)) for r in coefficients}
        lower = sum(c*prices[r][0 if c >= 0 else 1] for r, c in coefficients.items())
        upper = sum(c*prices[r][1 if c >= 0 else 0] for r, c in coefficients.items())
        require(lower > 0 if positive else upper < 0, 'The ceiling bracket fails')
        records.append({'kappa': value, 'margin_interval': outward((lower, upper)), 'positive': positive})
    return {'bracket': records, 'attained': False, 'new_kappa': False}


def main():
    authenticated = inputs()
    scalars = scalar_checks()
    result = {'scope': 'Cyclic common-stage rank budget; favorable central-query source ceiling',
              'scalar_rank_cases': len(scalars),
              'scalar_cases_sha256': hashlib.sha256(json.dumps(scalars).encode()).hexdigest(),
              'affine_words': product_checks(), 'complete_targets': target_checks(),
              'price': price_check(), 'authenticated_inputs': authenticated, 'new_kappa': False}
    print(json.dumps(result, indent=2))
    print('PASS_CYCLIC_STAGE_RANK_BUDGET')


if __name__ == '__main__':
    main()
