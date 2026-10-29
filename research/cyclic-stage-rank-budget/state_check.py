"""Check the nilpotent state proof in several characteristics.

Literal polynomial products and the separate state formula must agree.
Every test is finite corroboration of FORWARD-TARIFF.md.
"""
import itertools
from pathlib import Path
import random
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify import add, identity, multiply, require


def zero(rows, columns):
    return [[0]*columns for _ in range(rows)]


def polynomial_product(left, right, prime):
    result = [zero(len(left[0]), len(right[0][0])) for _ in range(len(left)+len(right)-1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] = add(result[i+j], multiply(a, b, prime))
    return [[[x % prime for x in row] for row in a] for a in result]


def constant_product(factors, first, last, n, prime):
    result = identity(n)
    for a, _, _ in factors[first:last]:
        result = multiply(a, result, prime)
    return result


def realization(factors, n, prime):
    widths = [len(x[0]) for _, x, _ in factors]
    offsets = [sum(widths[:i]) for i in range(len(widths))]
    size = sum(widths)
    t, b, c = zero(size, size), zero(size, n), zero(n, size)
    for j, (_, xj, yj) in enumerate(factors):
        bj = multiply(yj, constant_product(factors, 0, j, n, prime), prime)
        cj = multiply(constant_product(factors, j+1, len(factors), n, prime), xj, prime)
        for a in range(widths[j]):
            b[offsets[j]+a] = bj[a]
            for k in range(n):
                c[k][offsets[j]+a] = cj[k][a]
        for i, (_, xi, _) in enumerate(factors[:j]):
            block = multiply(multiply(yj, constant_product(factors, i+1, j, n, prime), prime), xi, prime)
            for a in range(widths[j]):
                for k in range(widths[i]):
                    t[offsets[j]+a][offsets[i]+k] = block[a][k]
    return t, b, c


def coefficients(t, b, c, prime):
    result, current = [], b
    for _ in range(len(t)):
        result.append(multiply(c, current, prime))
        current = multiply(t, current, prime)
    require(not any(x for row in current for x in row), 'The state is not nilpotent')
    return [zero(len(c), len(b[0]))]+result


def convolution(a, b, prime):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = (out[i+j]+x*y) % prime
    return out


def determinant(poly, prime):
    n = len(poly[0])
    out = [0]*(n*(len(poly)-1)+1)
    for permutation in itertools.permutations(range(n)):
        sign = (-1)**sum(permutation[i] > permutation[j] for i in range(n) for j in range(i+1, n))
        term = [sign % prime]
        for i, j in enumerate(permutation):
            term = convolution(term, [a[i][j] for a in poly], prime)
        for k, x in enumerate(term):
            out[k] = (out[k]+x) % prime
    return out


def check_state():
    rng, count, nonzero_determinants = random.Random(290630), 0, 0
    for prime in (2, 3, 5, 7):
        for n in (1, 2, 3):
            for length in (2, 4, 6):
                factors = []
                for i in range(length):
                    width = 1+i % min(n, 2)
                    a = [[rng.randrange(prime) for _ in range(n)] for _ in range(n)]
                    x = [[rng.randrange(prime) for _ in range(width)] for _ in range(n)]
                    y = [[rng.randrange(prime) for _ in range(n)] for _ in range(width)]
                    factors.append((a, x, y))
                literal = [identity(n)]
                for a, x, y in factors:
                    literal = polynomial_product([a, multiply(x, y, prime)], literal, prime)
                t, b, c = realization(factors, n, prime)
                state = coefficients(t, b, c, prime)
                state[0] = constant_product(factors, 0, length, n, prime)
                padded = literal+[zero(n, n)]*(len(state)-len(literal))
                require(state == padded, 'The state formula differs from the literal word')
                for period in range(2, 8):
                    # Form B_star by adding T^(a p) B directly.
                    star, current = zero(len(t), n), b
                    for power in range(len(t)):
                        if power % period == 0:
                            star = add(star, current)
                        current = multiply(t, current, prime)
                    star = [[x % prime for x in row] for row in star]
                    lifted = coefficients(t, star, c, prime)
                    for residue in range(1, period):
                        expected = zero(n, n)
                        for index in range(residue, len(literal), period):
                            expected = add(expected, literal[index])
                        expected = [[x % prime for x in row] for row in expected]
                        actual = lifted[residue] if residue < len(lifted) else zero(n, n)
                        require(actual == expected, 'Cyclic coefficient removal differs')
                    det = determinant(lifted, prime)
                    require(not any(det[len(t)+1:]), 'The block determinant degree bound fails')
                    nonzero_determinants += int(any(det))
                    count += 1
    require(nonzero_determinants > 0, 'Every determinant control is zero')
    return {'state_and_alias_cases': count, 'nonzero_determinants': nonzero_determinants,
            'characteristics': [2, 3, 5, 7], 'cyclic_orders': [2, 3, 4, 5, 6, 7]}


def direct_targets():
    records = []
    for prime in (2, 3, 5):
        for period, powers in ((3, [2, 2, 2, 2, 0]), (4, [3, 1, 2, 0]),
                               (5, [4, 2, 1, 0]), (6, [5, 3, 0])):
            n = len(powers)
            word = [identity(n)]
            for target, power in enumerate(powers):
                p = zero(n, n)
                p[target][target] = 1
                for _ in range(power):
                    word = polynomial_product([add(identity(n), p, -1), p], word, prime)
            expected = [zero(n, n) for _ in range(max(powers)+1)]
            for j, d in enumerate(powers):
                expected[d][j][j] = 1
            require(word[:len(expected)] == expected and not any(x for a in word[len(expected):] for row in a for x in row),
                    'A direct attaining word differs on an original or dirty input')
            records.append({'characteristic': prime, 'period': period,
                            'target_residues': powers, 'calls': sum(powers)})
    # At z=1 the nonfaithful representation hides every positive power.
    require([sum(a[i][j] for a in word) % prime for i in range(n) for j in range(n)]
            == [int(i == j) for i in range(n) for j in range(n)], 'The nonfaithful control differs')
    return records


def check():
    return check_state() | {'attaining_words': direct_targets(), 'new_kappa': False}


if __name__ == '__main__':
    import json
    print(json.dumps(check(), indent=2))
