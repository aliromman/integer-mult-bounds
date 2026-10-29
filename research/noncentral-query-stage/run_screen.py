"""Check the finite run patterns in the rational S3 boundary chart.

This exact scalar screen covers the stated field count. The rank reduction
does not by itself make this a theorem for every field count.
"""
import argparse
import json
from pathlib import Path
import sympy as s


def patterns(n, x):
    length = n+1
    for a in range(1, length):
        t = [1]*a+[x]*(length-a) if a % 2 else [x]*a+[1]*(length-a)
        yield [a, length-a], list(map(s.sympify, t))
    for a in range(1, length-1):
        for b in range(1, length-a):
            c = length-a-b
            if a % 2 and b % 2 and c % 2:
                first, middle, last = 4*x, 4*x*x, x
            elif a % 2:
                first, middle, last = 1, x, s.Rational(1, 4)
            elif b % 2:
                first, middle, last = 4*x, 1, x
            else:
                first, middle, last = 4, x, 1
            yield [a, b, c], list(map(s.sympify, [first]*a+[middle]*b+[last]*c))


def check(n):
    x = s.Symbol('x')
    length = n+1
    a, b, c, d = s.symbols('a b c d')
    values = (a, b, c, d)
    matrix = s.Matrix([[values[i] if i == j else 2*values[max(i, j)]-values[min(i, j)]
                        for j in range(4)] for i in range(4)])
    if s.expand(matrix[:3, :3].det()+3*(a-b)*(a-4*c)*(b-c)) or \
       s.expand(matrix[:3, 1:4].det()+6*(a-b)*(b-c)*(c-d)):
        raise ValueError('A symbolic rank-minor identity fails')
    records = []
    for runs, t in patterns(n, x):
        gram = s.Matrix([[0 if i == j else (-1 if i < j else t[j]/t[i])
                          for j in range(length)] for i in range(length)])
        k = gram[:n, :n]
        u, v = k.inv()*gram[:n, :], gram[:, :n]
        h = s.Matrix([[2*(value-1) for value in t]])
        r = k.inv()*s.Matrix([2*(1-1/value) for value in t[:n]])
        w = h[:, :n]
        boundary = s.eye(n)+r*w/2
        other_u = (u+r*(s.ones(1, length)+h/2)).applyfunc(s.cancel)
        other_v = v-s.ones(length, 1)*w
        word = s.eye(n)
        for site in range(length):
            word = ((s.eye(n)-2*other_u[:, site]*other_v[site, :])*word).applyfunc(s.cancel)
        numerators = [s.cancel(value).as_numer_denom()[0] for value in word+boundary if value != 0]
        if any(p.has(s.Float) for p in numerators):
            raise ValueError('The rational screen contains a floating coefficient')
        gcd = s.Poly(s.gcd_list(numerators), x, domain=s.QQ).monic().as_expr() if numerators else s.Integer(0)
        roots = s.solve(gcd, x) if gcd != 0 else []
        nonzero = [root for root in roots if root != 0]
        boundary_identity = [all(s.cancel(value).subs(x, root) == 0
                                 for value in boundary-s.eye(n)) for root in nonzero]
        row = {'runs': runs, 't': list(map(str, t)), 'gcd': str(s.factor(gcd)),
               'roots': list(map(str, roots)), 'nonzero_roots_force_identity': boundary_identity}
        records.append(row)
        print(json.dumps(row), flush=True)
    return {'n': n, 'patterns': records, 'exact_coefficient_domain': 'QQ',
            'rank_minor_identities': 'PASS_SYMBOLIC',
            'sympy_version': s.__version__, 'new_kappa': False,
            'scope': 'All two-run and admissible three-run cases in the stated rational chart'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = check(args.n)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
