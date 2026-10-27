"""Search one projector orbit on the signed simplex query frame.

The query frame is from rational_cyclic_search.py. This search is numeric.
It does not establish a caller or a kappa bound.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from rational_cyclic_search import frame


def search(args):
    n = args.n
    ident = np.eye(n)
    queries, affine = frame(n)
    columns = []
    for q in queries:
        j = int(np.argmax(np.sum(abs(q)**2, axis=0)))
        u = q[:, j]
        # The frame uses coordinate unit vectors for the first n columns.
        if len(columns) < n:
            u = ident[:, len(columns)]
        else:
            u = np.array([(-1.)**(n-i+1) for i in range(n)])
        columns.append(u)
    t = np.array(columns[1:]).T
    orbit = [(np.linalg.matrix_power(t, j), np.linalg.matrix_power(t, -j))
             for j in range(n+1)]
    if not np.allclose(np.linalg.matrix_power(t, n+1), -ident):
        raise ValueError('The signed cyclic frame fails')
    if any(not np.allclose(a @ queries[0] @ b, q)
           for (a, b), q in zip(orbit, queries)):
        raise ValueError('The query orbit fails')
    q = queries[0]
    u = ident[:, 0]
    v = q[0]
    pivot = int(np.argmax(abs(v)))
    a = ident[:, pivot]/v[pivot]
    cols = [u, a]
    for k in range(n):
        if k not in (0, pivot):
            cols.append(ident[:, k]-v[k]*a)
    basis = np.array(cols).T
    inverse = np.linalg.inv(basis)

    def projector(x):
        if args.complex:
            a, b = np.split(x, 2)
            x = a+1j*b
        if args.symplectic:
            p, ds = affine[0]
            return p+np.einsum('k,kij->ij', x, ds)
        a, b = np.split(x, 2)
        p = np.zeros((n, n), dtype=complex if args.complex else float)
        p[0, 0] = p[1, 1] = 1
        p[2:, 1], p[0, 2:], p[0, 1] = a, b, -b @ a
        return basis @ p @ inverse

    def endpoints(p, d):
        a, b = ident.astype(complex), ident.astype(complex)
        for (t, ti), q in zip(orbit, queries):
            stage = ident+(d-1)*(t @ p @ ti)
            a, b = stage @ a, stage @ (ident-2*q) @ b
        return a, b

    def residual(x):
        p = projector(x)
        probe = np.exp(2j*np.pi/args.period)
        ap, _ = endpoints(p, probe)
        e = (ap-probe**2*ident)/(probe**3-probe**2)
        parts = [e @ e-e, np.array([[np.trace(e)-2]])]
        if not args.complex:
            parts.append(1j*e.imag)
        for j in range(1, args.period if args.complex else args.period//2+1):
            d = np.exp(2j*np.pi*j/args.period)
            a, b = endpoints(p, d)
            target = d*d*(ident-e)+d*d*d*e
            parts.extend((a-target, b+target))
        return np.concatenate([part.ravel() for a in parts for part in (a.real, a.imag)])

    rng, records, best = np.random.default_rng(args.seed), [], None
    start = time.monotonic()
    width = (n-2)*(1 if args.symplectic else 2)*(2 if args.complex else 1)
    for trial in range(args.trials):
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
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY',
              'parameters': vars(args) | {'output': str(args.output)},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'cyclic_map': t.tolist()}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=8)
    p.add_argument('--period', type=int, default=7)
    p.add_argument('--symplectic', action='store_true')
    p.add_argument('--complex', action='store_true')
    p.add_argument('--trials', type=int, default=20)
    p.add_argument('--steps', type=int, default=500)
    p.add_argument('--scale', type=float, default=.3)
    p.add_argument('--seed', type=int, default=5833)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 4 or args.n % 2 or args.period < 4:
        p.error('Use even dimension at least four and period at least four')
    search(args)
