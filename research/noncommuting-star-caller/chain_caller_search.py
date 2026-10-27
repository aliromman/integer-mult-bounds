"""Search full-stage projectors for a fixed chain with one ordinary residue.

At stage value one, the n raw queries give S B+C. B is invertible and
C=(e0+e2+...)e0^T. This numerical search does not prove a full endpoint.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def search(args):
    n, ident = args.n, np.eye(args.n)
    queries, bases = [], []
    for i in range(n):
        u = ident[:, i]
        v = np.array([int(j == i-1)-int(j == i+1) for j in range(n)], dtype=float)
        queries.append(np.outer(u, v))
        pivot = i-1 if i else 1
        a = ident[:, pivot]/v[pivot]
        cols = [u, a]
        for j in range(n):
            if j not in (i, pivot):
                cols.append(ident[:, j]-v[j]*a)
        basis = np.array(cols).T
        bases.append((basis, np.linalg.inv(basis)))
    def raw(sign):
        value = ident.copy()
        for q in queries:
            value = (ident+sign*q) @ value
        return value
    b, c = (raw(1)-raw(-1))/2, (raw(1)+raw(-1))/2
    expected = np.zeros((n, n)); expected[::2, 0] = 1
    if not np.array_equal(c, expected) or abs(np.linalg.det(b)) != 1:
        raise ValueError('The fixed chain or its target basis fails')
    def unpack(x):
        if args.complex:
            a, z = np.split(x, 2)
            return a+1j*z
        return x
    def projectors(x):
        ps = []
        for row, (basis, inverse) in zip(unpack(x).reshape(n, 2*(n-2)), bases):
            a, z = row[:n-2], row[n-2:]
            p = np.zeros((n, n), dtype=complex if args.complex else float)
            p[0, 0] = p[1, 1] = 1
            p[2:, 1], p[0, 2:], p[0, 1] = a, z, -z @ a
            ps.append(basis @ p @ inverse)
        return ps
    def endpoint(ps, d, sign):
        value = ident.astype(complex)
        for stage, (p, q) in enumerate(zip(ps, queries)):
            step = 1/d if args.inverse_mask >> stage & 1 else d
            value = (ident+(step-1)*p+step*sign*q) @ value
        return value
    def residual(x):
        ps = projectors(x)
        values, e = [], None
        if args.signed_target:
            probe = np.exp(2j*np.pi/args.period)
            odd = (endpoint(ps, probe, 1)-endpoint(ps, probe, -1))/2
            e = (odd @ np.linalg.inv(b)-probe**-2*ident)/(probe**2-probe**-2)
            values.extend((e @ e-e, np.array([[np.trace(e)-n/2]])))
            if not args.complex:
                values.append(1j*e.imag)
        for j in range(1, args.period if args.complex else args.period//2+1):
            d = np.exp(2j*np.pi*j/args.period)
            target = d*d*b if e is None else (d*d*e+d**-2*(ident-e)) @ b
            values.extend((endpoint(ps, d, 1)-target-c,
                           endpoint(ps, d, -1)+target-c))
        return np.concatenate([part.ravel() for a in values for part in (a.real, a.imag)])
    rng = np.random.default_rng(args.seed)
    records, best = [], None
    start = time.monotonic()
    for trial in range(args.trials):
        width = n*2*(n-2)*(2 if args.complex else 1)
        x = np.zeros(width) if trial == 0 else rng.normal(0, args.scale, width)
        answer = least_squares(residual, x, max_nfev=args.steps,
                               ftol=1e-12, xtol=1e-12, gtol=1e-12)
        row = {'trial': trial, 'residual_norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(np.max(abs(answer.fun))), 'evaluations': answer.nfev,
               'max_parameter': float(np.max(abs(answer.x)))}
        records.append(row)
        if best is None or row['residual_norm'] < best['residual_norm']:
            best = row | {'parameters': answer.x.tolist()}
            print(json.dumps(row), flush=True)
        if best['max_residual'] < 1e-9:
            break
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'parameters': vars(args) | {'output': str(args.output)},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'query_matrices': [q.tolist() for q in queries], 'target_bank_matrix': b.tolist(),
              'ordinary_residual': c.tolist()}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=4)
    p.add_argument('--period', type=int, default=5)
    p.add_argument('--complex', action='store_true')
    p.add_argument('--inverse-mask', type=int, default=0)
    p.add_argument('--signed-target', action='store_true')
    p.add_argument('--trials', type=int, default=20)
    p.add_argument('--steps', type=int, default=500)
    p.add_argument('--scale', type=float, default=.5)
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 4 or args.n % 2 or args.period < 2:
        p.error('Use an even field count at least four and a period at least two')
    if args.inverse_mask >= 1 << args.n or (args.signed_target and 4 % args.period == 0):
        p.error('The inverse-stage mask or the signed target period is invalid')
    search(args)
