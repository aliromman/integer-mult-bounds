"""Search a three-full-call replacement for one full call.

The proposed replacement adds two independent D^2 Q fields. A finite exact
word would extend a caller without changing its query overhead. This script
only records numerical discovery, including unsuccessful searches.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def run(args):
    identity = np.eye(4)
    us = [np.array(x, float) for x in [(1, 0), (0, 1), (-1, 1)]]
    vs = [np.array(x, float) for x in [(0, -1), (1, 0), (1, 1)]]
    bases, queries = [], []
    for u, v in zip(us, vs):
        col, row = np.r_[1, 0, u], np.r_[0, 1, v]
        basis = np.column_stack((col, identity[:, 1], identity[:, 2]-v[0]*identity[:, 1],
                                 identity[:, 3]-v[1]*identity[:, 1]))
        bases.append((basis, np.linalg.inv(basis)))
        queries.append(np.outer(col, row))
    old = np.diag([1., 1., 0., 0.])
    new = identity-old
    old_query = np.outer(identity[:, 0], identity[:, 1])
    def projectors(x):
        result = []
        for row, (basis, inverse) in zip(x.reshape(3, 4), bases):
            a, b = row[:2], row[2:]
            p = np.zeros((4, 4))
            p[0, 0] = p[1, 1] = 1
            p[2:, 1], p[0, 2:], p[0, 1] = a, b, -b @ a
            result.append(basis @ p @ inverse)
        return result
    def residual(x):
        ps = projectors(x)
        errors = []
        for d in (-1, 1j):
            a, b = identity.astype(complex), identity.astype(complex)
            for p, q in zip(ps, queries):
                stage = identity+(d-1)*p
                a, b = stage @ a, stage @ (identity-2*q) @ b
            errors.extend([a-d*old-d*d*new, b-d*(old-2*old_query)+d*d*new])
        return np.concatenate([part.reshape(-1) for a in errors for part in (a.real, a.imag)])
    rng = np.random.default_rng(args.seed)
    best, records = None, []
    start = time.monotonic()
    for trial in range(args.trials):
        x = np.zeros(12) if trial == 0 else rng.normal(0, 1, 12)
        answer = least_squares(residual, x, max_nfev=args.steps,
                               ftol=1e-12, xtol=1e-12, gtol=1e-12)
        row = {'trial': trial, 'norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(max(abs(answer.fun))),
               'max_parameter': float(max(abs(answer.x))), 'evaluations': answer.nfev}
        records.append(row)
        if best is None or row['norm'] < best['norm']:
            best = row | {'parameters': answer.x.tolist()}
            print(json.dumps(row), flush=True)
        if best['max_residual'] < 1e-10:
            break
    data = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'seconds': time.monotonic()-start,
            'parameters': vars(args) | {'output': str(args.output)}, 'best': best, 'records': records,
            'queries': [q.astype(int).tolist() for q in queries]}
    args.output.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({'seconds': data['seconds'], 'status': data['status']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trials', type=int, default=30)
    p.add_argument('--steps', type=int, default=1000)
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args())
