"""Numerical search for full queries and projectors at the same time.

Each site has P^2=P, rank(P)=2, N^2=0, and PN=NP=N by construction.
The target is Q times a balanced polynomial in an order-four stage.
Numerical residuals do not certify an endpoint or a kappa value.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def solve(args):
    n, calls = args.n, args.calls or args.n+2
    width = 4*(n-2)+2
    target_rank = 0 if args.uniform else 2*(calls-n)
    if not 0 <= target_rank <= n:
        raise ValueError('The balanced endpoint does not fit')
    ident = np.eye(n)
    pairs = [(2*(s//3) % n, (2*(s//3)+1) % n) for s in range(calls)]
    permutations = [[a, b]+[k for k in range(n) if k not in (a, b)] for a, b in pairs]
    def local(row, order):
        a = row[:2*(n-2)].reshape(n-2, 2)
        b = row[2*(n-2):4*(n-2)].reshape(2, n-2)
        t, scale = row[-2:]
        u = np.vstack((np.eye(2), a))
        v = np.hstack((np.eye(2)-b @ a, b))
        j = scale*np.array([[-t, 1.], [-t*t, t]])
        p, q = u @ v, u @ j @ v
        result_p, result_q = np.zeros((n, n)), np.zeros((n, n))
        result_p[np.ix_(order, order)] = p
        result_q[np.ix_(order, order)] = q
        return result_p, result_q
    def maps(x):
        return [local(row, order) for row, order in zip(x.reshape(calls, width), permutations)]
    def endpoints(ps, value):
        a, b = ident.astype(complex), ident.astype(complex)
        for p, q in ps:
            stage = ident+(value-1)*p
            a, b = stage @ a, stage @ (ident-2*q) @ b
        return a, b
    def residual(x):
        ps = maps(x)
        a0, b0 = endpoints(ps, 1)
        am, bm = endpoints(ps, -1)
        ai, bi = endpoints(ps, 1j)
        e = (ident-am)/2
        values = [b0+ident, am+bm, ai+bi]
        if args.balanced:
            values.extend([e @ e-e, ai+ident-(1-1j)*e])
        result = np.concatenate([part.reshape(-1) for a in values for part in (a.real, a.imag)])
        if args.balanced:
            result = np.append(result, np.trace(e).real-target_rank)
        return result
    rng = np.random.default_rng(args.seed)
    best, records = None, []
    start = time.monotonic()
    for trial in range(args.trials):
        x = rng.normal(0, args.scale, (calls, width))
        x[:, -1] += 1
        answer = least_squares(residual, x.reshape(-1), max_nfev=args.steps,
                               ftol=1e-11, xtol=1e-11, gtol=1e-11)
        row = {'trial': trial, 'norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(np.max(abs(answer.fun))),
               'max_parameter': float(np.max(abs(answer.x))), 'evaluations': answer.nfev}
        records.append(row)
        if best is None or row['norm'] < best['norm']:
            best = row | {'parameters': answer.x.tolist()}
            print(json.dumps(row), flush=True)
        if best['max_residual'] < 1e-9:
            break
    ps = maps(np.array(best['parameters']))
    for p, q in ps:
        if (not np.allclose(p @ p, p) or not np.allclose(q @ q, 0)
                or not np.allclose(p @ q, q) or not np.allclose(q @ p, q)):
            raise ValueError('Numerical projector parametrization fails')
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'parameters': vars(args) | {'output': str(args.output)},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'permutations': permutations,
              'projectors': [p.tolist() for p, _ in ps], 'queries': [q.tolist() for _, q in ps]}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'seconds': result['seconds'], 'status': result['status']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--calls', type=int)
    p.add_argument('--balanced', action='store_true')
    p.add_argument('--uniform', action='store_true', help='Require D^2 Q on every field; uses D^4=1')
    p.add_argument('--trials', type=int, default=4)
    p.add_argument('--steps', type=int, default=250)
    p.add_argument('--scale', type=float, default=.2)
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 4 or args.n % 2:
        p.error('Use an even dimension of at least four')
    if args.uniform:
        args.balanced = True
    solve(args)
