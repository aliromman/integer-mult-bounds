"""Construct a rational shortest-query chart for the S3 standard block.

The sign representation is a separate equation. Passing the standard
block does not establish a complete group-algebra query identity.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
from exact_replay import inverse


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def construct(t):
    length = len(t)
    n = length-1
    g = [[0 if i == j else (-1 if i < j else t[j]/t[i])
          for j in range(length)] for i in range(length)]
    inverse_k = inverse([row[:n] for row in g[:n]])
    v = [row[:n] for row in g]
    u = mm(inverse_k, g[:n])
    h = [2*(x-1) for x in t]
    target = [[2*(1-1/x)] for x in t]
    r = mm(inverse_k, target[:n])
    w = h[:n]
    if mm(v, r) != target or sum(w[i]*r[i][0] for i in range(n)):
        raise ValueError('The rank-one boundary conditions fail')
    other_u = [[u[i][j]+r[i][0]*(1+h[j]/2) for j in range(length)] for i in range(n)]
    other_v = [[v[i][j]-w[j] for j in range(n)] for i in range(length)]
    boundary = [[F(i == j)+r[i][0]*w[j]/2 for j in range(n)] for i in range(n)]
    return u, v, other_u, other_v, boundary


def replay(t):
    u, v, other_u, other_v, boundary = construct(t)
    n = len(u)
    ident = [[F(i == j) for j in range(n)] for i in range(n)]
    signed = ident
    standard = [[F(i == j) for j in range(2*n)] for i in range(2*n)]
    s = [[F(-1, 2), F(3, 2)], [F(1, 2), F(1, 2)]]
    for site in range(n+1):
        us, vs = (u, other_u), (v, other_v)
        if any(sum(us[a][j][site]*vs[a][site][j] for j in range(n)) for a in range(2)):
            raise ValueError('A sector query is not square zero')
        query = [[F(i == j)+(s[i//n][j//n]-int(i//n == j//n))*
                  us[i//n][i % n][site]*vs[j//n][site][j % n]
                  for j in range(2*n)] for i in range(2*n)]
        standard = mm(query, standard)
        signed = mm([[ident[i][j]-2*other_u[i][site]*other_v[site][j]
                      for j in range(n)] for i in range(n)], signed)
    target = [[s[i//n][j//n]*(ident if i < n else boundary)[i % n][j % n]
               for j in range(2*n)] for i in range(2*n)]
    if standard != target:
        raise ValueError('The complete standard-representation word fails')
    error = [[signed[i][j]+boundary[i][j] for j in range(n)] for i in range(n)]
    return {'n': n, 't': list(map(str, t)), 'standard': 'PASS',
            'sign': 'PASS' if not any(x for row in error for x in row) else 'FAIL',
            'sign_error': [[str(x) for x in row] for row in error], 'new_kappa': False}


def ideal(n, seconds):
    import sympy as sp
    start = time.monotonic()
    variables = sp.symbols(f't0:{n}')
    t = list(variables)+[sp.prod(variables[1::2])/sp.prod(variables[::2])]
    g = sp.Matrix([[0 if i == j else (-1 if i < j else t[j]/t[i])
                    for j in range(n+1)] for i in range(n+1)])
    k = g[:n, :n]
    u, v = k.inv()*g[:n, :], g[:, :n]
    h = sp.Matrix([[2*(x-1) for x in t]])
    r = k.inv()*sp.Matrix([2*(1-1/x) for x in t[:n]])
    w = h[:, :n]
    boundary = sp.eye(n)+r*w/2
    other_u = (u+r*(sp.ones(1, n+1)+h/2)).applyfunc(sp.cancel)
    other_v = v-sp.ones(n+1, 1)*w
    word = sp.eye(n)
    for site in range(n+1):
        word = ((sp.eye(n)-2*other_u[:, site]*other_v[site, :])*word).applyfunc(sp.cancel)
    equations = [sp.Poly(sp.cancel(x).as_numer_denom()[0], *variables).as_expr()
                 for x in word+boundary if x != 0]
    equations = list(dict.fromkeys(equations))
    inverse_variable = sp.Symbol('inverse_product')
    equations.append(inverse_variable*sp.prod(variables)-1)
    print(json.dumps({'n': n, 'equations': len(equations),
                      'max_degree': int(max(sp.total_degree(p) for p in equations)),
                      'term_count': sum(len(sp.Poly(p, *variables, inverse_variable).terms()) for p in equations)}), flush=True)

    def expire(signum, frame):
        raise TimeoutError('The boundary ideal time limit expired')

    signal.signal(signal.SIGALRM, expire)
    signal.alarm(seconds)
    result = {'n': n, 'method': 'Exact QQ Groebner basis with nonzero-parameter saturation',
              'sympy_version': sp.__version__, 'new_kappa': False}
    try:
        basis = sp.groebner(equations, inverse_variable, *variables, order='grevlex', method='f5b')
        numerators = [sp.cancel(x).as_numer_denom()[0] for x in boundary-sp.eye(n)]
        residues = [basis.reduce(x)[1] for x in numerators]
        result.update(status='EXACT_IDEAL_COMPUTED', basis=list(map(str, basis.polys)),
                      boundary_residues=list(map(str, residues)),
                      boundary_forced_identity=not any(residues))
    except TimeoutError:
        result['status'] = 'UNKNOWN_GROEBNER_TIMEOUT'
    finally:
        signal.alarm(0)
    result['seconds'] = time.monotonic()-start
    result['limit_seconds'] = seconds
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--ideal', action='store_true')
    parser.add_argument('--seconds', type=int, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.ideal:
        result = ideal(args.n, args.seconds)
    else:
        t = list(map(F, range(2, args.n+2)))
        t.append(math.prod(t[1::2])/math.prod(t[::2]))
        result = replay(t)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
