"""Search an S3 query word with a rank-one ordinary boundary shear.

The three irreducible representations are included. The trivial one is
automatic. A numerical zero still needs exact reconstruction and stages.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'five-call-projector-rigidity'))
from exact_replay import canonical


def model(n, beta=1):
    length = n+1
    count = length*4*n
    ident = np.eye(n, dtype=complex)
    s = np.array([[-.5, 1.5], [.5, .5]], dtype=complex)
    coefficient = s-np.eye(2)
    boundary = ident.copy()
    boundary[0, 1] = beta
    target = np.block([[s[0, 0]*ident, s[0, 1]*ident],
                       [s[1, 0]*boundary, s[1, 1]*boundary]])

    def calculate(x):
        parameters = (x[:count]+1j*x[count:]).reshape(length, 2, 2, n)
        word = np.eye(2*n, dtype=complex)
        derivative = np.zeros((count, 2*n, 2*n), dtype=complex)
        signed = ident.copy()
        dsigned = np.zeros((count, n, n), dtype=complex)
        nilpotents, dnilpotents = [], []
        for site in range(length):
            us = parameters[site, :, 0]
            vs = parameters[site, :, 1]
            t = np.eye(2*n, dtype=complex)
            dt = np.zeros_like(derivative)
            for a in range(2):
                for b in range(2):
                    rows, cols = slice(a*n, (a+1)*n), slice(b*n, (b+1)*n)
                    t[rows, cols] += coefficient[a, b]*np.outer(us[a], vs[b])
                    for j in range(n):
                        dt[(4*site+2*a)*n+j, rows, cols][j, :] += coefficient[a, b]*vs[b]
                        dt[(4*site+2*b+1)*n+j, rows, cols][:, j] += coefficient[a, b]*us[a]
            ts = ident-2*np.outer(us[1], vs[1])
            dts = np.zeros_like(dsigned)
            for j in range(n):
                dts[(4*site+2)*n+j, j, :] = -2*vs[1]
                dts[(4*site+3)*n+j, :, j] = -2*us[1]
            derivative = dt@word+t@derivative
            word = t@word
            dsigned = dts@signed+ts@dsigned
            signed = ts@signed
            for a in range(2):
                nilpotents.append(us[a]@vs[a])
                row = np.zeros(count, dtype=complex)
                row[(4*site+2*a)*n:(4*site+2*a+1)*n] = vs[a]
                row[(4*site+2*a+1)*n:(4*site+2*a+2)*n] = us[a]
                dnilpotents.append(row)
        values = np.r_[(word-target).ravel(), (signed+boundary).ravel(), nilpotents]
        jacobian = np.vstack([derivative.reshape(count, -1).T,
                              dsigned.reshape(count, -1).T, dnilpotents])
        return np.r_[values.real, values.imag], np.block([[jacobian.real, -jacobian.imag],
                                                        [jacobian.imag, jacobian.real]])

    us, vs, _ = canonical(n)
    initial = np.array([[list(zip(us, vs))[i], list(zip(us, vs))[i]]
                        for i in range(length)], dtype=float).ravel()
    initial = np.r_[initial, np.zeros(count)]
    return count, initial, calculate


def search(args):
    count, initial, calculate = model(args.n)
    _, control, positive = model(args.n, beta=0)
    if np.any(positive(control)[0]) or not np.any(calculate(initial)[0]):
        raise ValueError('The canonical control or the changed boundary fails')
    rng = np.random.default_rng(args.seed)
    probe = initial+rng.normal(0, .01, 2*count)
    _, jacobian = calculate(probe)
    for j in np.linspace(0, 2*count-1, 24, dtype=int):
        delta = np.zeros(2*count)
        delta[j] = 1e-6
        difference = (calculate(probe+delta)[0]-calculate(probe-delta)[0])/2e-6
        if np.max(abs(difference-jacobian[:, j])) > 1e-4:
            raise ValueError('The analytic derivative differs from finite differences')
    start = time.monotonic()
    records, best = [], None
    for trial in range(args.trials):
        x = initial if trial == 0 else initial+rng.normal(0, .15, 2*count)
        result = least_squares(lambda p: calculate(p)[0], x,
                               jac=lambda p: calculate(p)[1], max_nfev=args.steps,
                               ftol=1e-11, xtol=1e-11, gtol=1e-11)
        row = {'trial': trial, 'norm': float(np.linalg.norm(result.fun)),
               'max_residual': float(np.max(abs(result.fun))), 'evaluations': result.nfev}
        records.append(row)
        print(json.dumps(row), flush=True)
        if best is None or row['norm'] < best['norm']:
            best = row | {'parameters': result.x.tolist()}
        if row['max_residual'] < 1e-10:
            break
    return {'n': args.n, 'queries': args.n+1, 'boundary': 'I+E01 on J-minus; I on J-plus',
            'records': records, 'best': best, 'seconds': time.monotonic()-start,
            'canonical_zero_boundary_control': 'PASS', 'changed_boundary_rejected': True,
            'derivative_control': 'PASS_24_COLUMNS', 'new_kappa': False,
            'scope': 'Numerical query-only S3 group-algebra search; no stage or native caller'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--trials', type=int, default=4)
    parser.add_argument('--steps', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=290627)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = search(args)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('best', 'records')}))
