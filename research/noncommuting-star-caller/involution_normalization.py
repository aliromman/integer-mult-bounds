"""Check the exact central involution and eight-axis stage normalization.

The proof for every width is in ADVANCING-CONTRACT.md. These integer
checks cover complete kernels at the stated widths. They do not check
a complete caller, native instruction word, or recursive invoice.
"""
import argparse
import json
from pathlib import Path


def mul(a, b):
    x, y = a
    u, v = b
    return x*u-y*v, x*v+y*u


def power(a, n):
    result = (1, 0)
    for _ in range(n):
        result = mul(result, a)
    return result


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def kernel(f):
    size = 1 << f
    q = [mul(power((1, 1), f-x.bit_count()), power((1, -1), x.bit_count())) for x in range(size)]
    result = [add(mul(value, (1, -1)), mul(q[x ^ (size-1)], (1, 1))) for x, value in enumerate(q)]
    if any(y for _, y in result):
        raise ValueError('The normalized kernel is not real')
    return [x for x, _ in result], 1 << (f+1)


def check(max_width):
    rows = []
    for f in range(1, max_width+1):
        values, den = kernel(f)
        size = len(values)
        for z in range(size):
            square = sum(values[x]*values[x ^ z] for x in range(size))
            if square != (den*den if z == 0 else 0):
                raise ValueError('The square is not the identity')
            eigenvalue = sum((-1 if (x & z).bit_count() % 2 else 1)*values[x] for x in range(size))
            w = z.bit_count()
            if eigenvalue != den*((-1)**(w*(w-1)//2)):
                raise ValueError('The complete Fourier eigenvalue fails')
        rows.append({'width': f, 'addresses': size, 'checked_kernel_pairs': size*size,
                     'negative_eigenvalues': sum((z.bit_count()*(z.bit_count()-1)//2) % 2 for z in range(size))})
    for x in range(256):
        for y in range(256):
            d = (x ^ y).bit_count()
            entry = mul(power((1, 1), 8-d), power((1, -1), d))
            entry = mul(entry, power((0, 1), x.bit_count()+y.bit_count()))
            if entry != (16*(-1 if (x & y).bit_count() % 2 else 1), 0):
                raise ValueError('The eight-axis Hadamard phase fails')
    a, b = (1, -1), (1, 1)
    # G has denominator two. Check both entries of G G^-1.
    if add(mul(a, b), mul(b, a)) != (4, 0) or add(mul(a, a), mul(b, b)) != (0, 0):
        raise ValueError('The two-address inverse fails')
    one, den1 = kernel(1)
    two, den2 = kernel(2)
    if one != [den1, 0] or two == [den2, 0, 0, 0]:
        raise ValueError('The factorwise-normalization control fails')
    # Q G^-1 equals X at width one. It is not the selected Q G = I.
    wrong = [add(mul(q, (1, 1)), mul(r, (1, -1))) for q, r in (((1, 1), (1, -1)), ((1, -1), (1, 1)))]
    if wrong == [(den1, 0), (0, 0)]:
        raise ValueError('The inverse-substitution control does not differ')
    return {'status': 'PASS_CENTRAL_INVOLUTION_NORMALIZATION', 'new_kappa': False,
            'width_checks': rows, 'eight_axis_entries': 65536,
            'controls': ['factorwise G differs at width two', 'G inverse differs at width one'],
            'scope': 'Exact component identities; no complete native caller or price improvement'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-width', type=int, default=6)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 2 <= args.max_width <= 10:
        parser.error('Use a width from two through ten')
    result = check(args.max_width)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
