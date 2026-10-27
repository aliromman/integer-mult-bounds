"""Numerical raw-query caller with arbitrary swap-dependent full frames.

All five complex irreducible representations of the dihedral group are
checked. The target is U^2 S K(R)+C(R). No numerical receipt proves it.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def solve(args):
    n, sites, local_width = 4, 4, 12
    ident, j = np.eye(n), np.array([[0., 1.], [0., 0.]])
    orders = [[a, b]+[k for k in range(n) if k not in (a, b)] for a, b in [(0, 1), (1, 0), (2, 3), (3, 2)]]
    def frames(x, signs):
        out = [(np.eye(n)[:, :2], np.eye(n)[:2])]*2
        result = [out]
        rows = x[:72].reshape(3, 2, local_width)
        for stage, pair in enumerate(rows, 1):
            site = []
            for side, row in enumerate(pair):
                a, b = row[:4].reshape(2, 2), row[4:8].reshape(2, 2)
                c, d, log_s, log_t = row[8:]
                g = np.array([[1+c*d, c], [d, 1.]]) @ np.diag([np.exp(log_s), signs[stage-1, side]*np.exp(log_t)])
                u = np.vstack((np.eye(2), a)) @ g
                v = np.linalg.solve(g, np.hstack((np.eye(2)-b @ a, b)))
                uu, vv = np.zeros((n, 2)), np.zeros((2, n))
                uu[orders[stage]], vv[:, orders[stage]] = u, v
                site.append((uu, vv))
            result.append(site)
        return result
    def endpoints(fs):
        small = {}
        for side in range(2):
            for value in (1, -1):
                for sign in (1, -1):
                    word = ident.copy()
                    for pair in fs:
                        u, v = pair[side]
                        word = (ident+(value-1)*(u @ v)+sign*value*(u @ j @ v)) @ word
                    small[side, value, sign] = word
        large = {}
        advance = np.array([[0., 1.], [-1., 0.]])
        for sign in (1, -1):
            word = np.eye(2*n)
            for pair in fs:
                blocks = [[(ident if a == b else np.zeros((n, n)))
                           +(advance[a, b]-(a == b))*(pair[a][0] @ pair[b][1])
                           +sign*advance[a, b]*(pair[a][0] @ j @ pair[b][1])
                           for b in range(2)] for a in range(2)]
                word = np.block(blocks) @ word
            large[sign] = word
        return small, large
    def residual(x, signs):
        small, large = endpoints(frames(x, signs))
        kk, cc, values = [], [], []
        inverses = x[72:].reshape(2, n, n)
        for side in range(2):
            plus, minus = small[side, 1, 1], small[side, 1, -1]
            k, c = (plus-minus)/2, (plus+minus)/2
            kk.append(k)
            cc.append(c)
            values.extend((inverses[side] @ k-ident, small[side, -1, 1]-k-c,
                           small[side, -1, -1]+k-c))
        zeros = np.zeros((n, n))
        big_k, big_c = np.block([[kk[0], zeros], [zeros, kk[1]]]), np.block([[cc[0], zeros], [zeros, cc[1]]])
        csign = -1 if args.relaxed_control else 1
        values.extend((large[1]+big_k-csign*big_c, large[-1]-big_k-csign*big_c))
        return np.concatenate([a.reshape(-1) for a in values])
    rng = np.random.default_rng(args.seed)
    best, records = None, []
    start = time.monotonic()
    for trial in range(args.trials):
        signs = np.ones((3, 2)) if trial == 0 else rng.choice([-1., 1.], size=(3, 2))
        params = np.zeros(72) if trial == 0 else rng.normal(0, args.scale, 72)
        x = np.r_[params, ident.reshape(-1), ident.reshape(-1)]
        small, _ = endpoints(frames(x, signs))
        for side in range(2):
            k = (small[side, 1, 1]-small[side, 1, -1])/2
            if abs(np.linalg.det(k)) > 1e-8:
                x[72+side*n*n:72+(side+1)*n*n] = np.linalg.inv(k).reshape(-1)
        # Bound only the logarithmic scale coordinates to prevent overflow.
        lower, upper = np.full(len(x), -np.inf), np.full(len(x), np.inf)
        for row in range(6):
            lower[row*12+10:row*12+12], upper[row*12+10:row*12+12] = -8, 8
        answer = least_squares(lambda y: residual(y, signs), x, bounds=(lower, upper),
                               max_nfev=args.steps, ftol=1e-11, xtol=1e-11, gtol=1e-11)
        row = {'trial': trial, 'residual_norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(np.max(abs(answer.fun))), 'evaluations': answer.nfev,
               'max_parameter': float(np.max(abs(answer.x)))}
        records.append(row)
        if best is None or row['residual_norm'] < best['residual_norm']:
            best = row | {'parameters': answer.x.tolist(), 'gauge_signs': signs.tolist()}
            print(json.dumps(row), flush=True)
        if best['max_residual'] < 1e-9:
            break
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'parameters': vars(args) | {'output': str(args.output)},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'representations': 'four scalar characters and one two-dimensional representation',
              'scope': 'Real swap-dependent frames in stated open charts; no exact identity or price claim'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trials', type=int, default=12)
    p.add_argument('--steps', type=int, default=500)
    p.add_argument('--scale', type=float, default=.15)
    p.add_argument('--relaxed-control', action='store_true', help='Permit a Gaussian residual; this has no new target credit')
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    solve(p.parse_args())
