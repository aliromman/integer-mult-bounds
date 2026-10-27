"""Search a rectangular caller with initialized and discarded work fields.

All factors have compatible rank-two projectors and square-zero queries.
This is numerical discovery. It does not prove a caller or a kappa bound.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def search(args):
    n, m, calls = args.n, args.n+args.work, args.calls
    ident, embed = np.eye(m), np.eye(m)[:, :n]
    degrees = [int(d) for d in args.degrees.split(',')] if args.degrees else [args.exponent]*n
    if len(degrees) != n or min(degrees) < 0 or max(degrees) > calls:
        raise ValueError('The target degrees do not fit the polynomial word')
    width = 4*(m-2)+2
    orders = []
    for i in range(calls):
        pair = [0, 1] if n == 2 else [i % n, (i+1) % n]
        orders.append(pair+[j for j in range(m) if j not in pair])

    def unpack(values):
        if args.complex:
            a, b = np.split(values, 2)
            values = a+1j*b
        stages, left = values[:calls*width].reshape(calls, width), values[calls*width:].reshape(n, m)
        factors = []
        for i, row in enumerate(stages):
            a = row[:2*(m-2)].reshape(m-2, 2)
            b = row[2*(m-2):4*(m-2)].reshape(2, m-2)
            c, t = row[-2:]
            x = np.concatenate((np.eye(2), a))
            y = np.concatenate((np.eye(2)-b @ a, b), axis=1)
            u, v = (np.array([1, t]), np.array([-t, 1])) if i % 2 == 0 else (np.array([t, 1]), np.array([1, -t]))
            p, q = x @ y, c*x @ np.outer(u, v) @ y
            order = np.argsort(orders[i])
            factors.append((p[np.ix_(order, order)], q[np.ix_(order, order)]))
        return factors, left

    points = [np.exp(2j*np.pi*j/(calls+1))
              for j in range(calls+1 if args.complex else (calls+1)//2+1)]

    def residual(values):
        factors, left = unpack(values)
        parts = []
        for d in points:
            for sign in (1, -1):
                value = embed.astype(complex)
                for p, q in factors:
                    value = (ident+(d-1)*p+d*sign*q) @ value
                parts.append(left @ value-sign*np.diag([d**degree for degree in degrees]))
        return np.concatenate([part.ravel() for a in parts for part in (a.real, a.imag)])

    rng, best, records = np.random.default_rng(args.seed), None, []
    start = time.monotonic()
    size = calls*width+n*m
    for trial in range(args.trials):
        x = np.zeros(size*(2 if args.complex else 1))
        for i in range(calls):
            x[i*width+width-2] = 1 if i % 2 == 0 else -1
        x[calls*width:size] = embed.T.ravel()
        if n == 2 and calls == 3:
            x[calls*width:size] = np.array([[0, -1]+[0]*(m-2), [1, 0]+[0]*(m-2)]).ravel()
        if trial:
            x += rng.normal(0, args.scale, len(x))
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
    factors, left = unpack(np.array(best['parameters']))
    def encode(a):
        return {'real': a.real.tolist(), 'imag': a.imag.tolist()}
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY',
              'parameters': vars(args) | {'output': str(args.output)},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'projectors': [encode(p) for p, _ in factors],
              'queries': [encode(q) for _, q in factors], 'output_map': encode(left),
              'polynomial_degree': calls, 'root_count': calls+1, 'target_degrees': degrees}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=4)
    p.add_argument('--work', type=int, default=2)
    p.add_argument('--calls', type=int, default=4)
    p.add_argument('--exponent', type=int, default=2)
    p.add_argument('--degrees', help='Set one target stage degree for each output')
    p.add_argument('--complex', action='store_true')
    p.add_argument('--trials', type=int, default=12)
    p.add_argument('--steps', type=int, default=600)
    p.add_argument('--scale', type=float, default=.2)
    p.add_argument('--seed', type=int, default=5838)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 2 or args.work < 0 or args.calls < args.exponent or args.exponent < 1:
        p.error('Use at least two inputs and valid work, call and exponent counts')
    search(args)
