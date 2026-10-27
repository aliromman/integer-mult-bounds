"""Search a Gaussian caller with a removable ordinary input residual.

The target is U^exponent S B + C. B is an invertible constant bank map.
C is constant in U and S. Original input copies can remove C. This
numerical test does not construct that consumer or certify its cost.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def solve(args):
    n, calls = args.n, args.calls or args.n
    width = 3*(n-2)+1 if args.image_basis else 4*(n-2)+2
    pairs = [(s % n, (s+1) % n) for s in range(calls)]
    orders = [[a, b]+[k for k in range(n) if k not in (a, b)] for a, b in pairs]
    ident = np.eye(n)
    def unpack(x):
        if args.complex_parameters:
            a, b = np.split(x, 2)
            return a+1j*b
        return x
    def pack(x):
        return np.r_[x.real, x.imag] if args.complex_parameters else x
    def maps(x):
        result = []
        for stage, (row, order) in enumerate(zip(x[:calls*width].reshape(calls, width), orders)):
            if args.image_basis:
                i, j = stage, (stage+1) % n
                other = [k for k in range(n) if k not in (i, j)]
                u = ident[:, i]
                z, v, w = [np.zeros(n, dtype=x.dtype) for _ in range(3)]
                z[j], w[i] = 1, 1
                z[other] = row[:n-2]
                v[other] = row[n-2:2*(n-2)]
                w[other] = row[2*(n-2):3*(n-2)]
                v[j], w[j] = 1-v @ z, -w @ z
                result.append((np.outer(u, w)+np.outer(z, v), row[-1]*np.outer(u, v)))
                continue
            a, b = row[:2*(n-2)].reshape(n-2, 2), row[2*(n-2):4*(n-2)].reshape(2, n-2)
            t, scale = row[-2:]
            u = np.vstack((np.eye(2), a))
            v = np.hstack((np.eye(2)-b @ a, b))
            j = scale*np.array([[-t, 1.], [-t*t, t]])
            p, q = np.zeros((n, n), dtype=x.dtype), np.zeros((n, n), dtype=x.dtype)
            p[np.ix_(order, order)], q[np.ix_(order, order)] = u @ v, u @ j @ v
            result.append((p, q))
        return result
    def endpoint(ps, d, sign):
        value = ident.astype(complex)
        for stage, (p, q) in enumerate(ps):
            step = 1/d if args.inverse_mask >> stage & 1 else d
            value = (ident+(step-1)*p+step*sign*q) @ value
        return value
    def residual(x):
        x = unpack(x)
        ps = maps(x)
        positive, negative = endpoint(ps, 1, 1), endpoint(ps, 1, -1)
        b, c = (positive-negative)/2, (positive+negative)/2
        inverse = x[calls*width:].reshape(n, n)
        values = [inverse @ b-ident]
        e = None
        if args.signed_target:
            probe = np.exp(2j*np.pi/args.period)
            odd = (endpoint(ps, probe, 1)-endpoint(ps, probe, -1))/2
            e = (odd @ inverse-probe**(-args.exponent)*ident)/(probe**args.exponent-probe**(-args.exponent))
            values.extend((e @ e-e, np.array([[np.trace(e)-n/2]])))
            if not args.complex_parameters:
                values.append(1j*e.imag)
        for j in range(1, args.period if args.complex_parameters else args.period//2+1):
            d = np.exp(2j*np.pi*j/args.period)
            target = d**args.exponent*b if e is None else (d**args.exponent*e+d**(-args.exponent)*(ident-e)) @ b
            errors = [endpoint(ps, d, 1)-target-c, endpoint(ps, d, -1)+target-c]
            values.extend(error @ inverse if args.relative_residual else error for error in errors)
        return np.concatenate([part.reshape(-1) for value in values for part in (value.real, value.imag)])
    rng = np.random.default_rng(args.seed)
    records, best = [], None
    start = time.monotonic()
    for trial in range(args.trials):
        params = np.zeros((calls, width)) if trial == 0 else rng.normal(0, args.scale, (calls, width))
        if args.complex_parameters:
            params = params.astype(complex)
            if trial:
                params += 1j*rng.normal(0, args.scale, (calls, width))
        params[:, -1] += 1
        x = np.r_[params.reshape(-1), ident.reshape(-1)]
        ps = maps(x)
        b = (endpoint(ps, 1, 1)-endpoint(ps, 1, -1))/2
        if abs(np.linalg.det(b)) > 1e-8:
            x[calls*width:] = (np.linalg.inv(b) if args.complex_parameters else np.linalg.inv(b).real).reshape(-1)
        x = pack(x)
        if trial == 0 and args.start:
            x = np.array(json.loads(args.start.read_text())['best']['parameters'])
            if len(x) != (calls*width+n*n)*(2 if args.complex_parameters else 1):
                raise ValueError('The saved parameter count does not match')
        answer = least_squares(residual, x, max_nfev=args.steps,
                               ftol=1e-11, xtol=1e-11, gtol=1e-11)
        ps = maps(unpack(answer.x))
        b = (endpoint(ps, 1, 1)-endpoint(ps, 1, -1))/2
        row = {'trial': trial, 'residual_norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(np.max(abs(answer.fun))), 'evaluations': answer.nfev,
               'max_parameter': float(np.max(abs(answer.x))), 'target_condition': float(np.linalg.cond(b))}
        records.append(row)
        if best is None or row['residual_norm'] < best['residual_norm']:
            best = row | {'parameters': answer.x.tolist()}
            print(json.dumps(row), flush=True)
        if best['max_residual'] < 1e-9:
            break
    ps = maps(unpack(np.array(best['parameters'])))
    plus, minus = endpoint(ps, 1, 1), endpoint(ps, 1, -1)
    def encoded(a):
        return {'real': a.real.tolist(), 'imag': a.imag.tolist()} if args.complex_parameters else a.real.tolist()
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'parameters': vars(args) | {'output': str(args.output), 'start': str(args.start) if args.start else None},
              'seconds': time.monotonic()-start, 'best': best, 'records': records,
              'projectors': [encoded(p) for p, _ in ps], 'queries': [encoded(q) for _, q in ps],
              'target_bank_matrix': encoded((plus-minus)/2), 'ordinary_residual': encoded((plus+minus)/2)}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=4)
    p.add_argument('--calls', type=int)
    p.add_argument('--period', type=int, default=6)
    p.add_argument('--exponent', type=int, default=2)
    p.add_argument('--inverse-mask', type=int, default=0)
    p.add_argument('--signed-target', action='store_true', help='Split the target between positive and negative stage powers')
    p.add_argument('--complex-parameters', action='store_true')
    p.add_argument('--image-basis', action='store_true', help='Fix the independent full-query images to the coordinate basis')
    p.add_argument('--relative-residual', action='store_true', help='Measure endpoint errors after the target input adapter')
    p.add_argument('--trials', type=int, default=12)
    p.add_argument('--steps', type=int, default=400)
    p.add_argument('--scale', type=float, default=.2)
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--start', type=Path, help='Continue the best parameters from a previous receipt')
    args = p.parse_args()
    if args.n < 2 or args.period < 2:
        p.error('Use at least two fields and a stage order of at least two')
    if args.signed_target and (args.n % 2 or (2*args.exponent) % args.period == 0):
        p.error('The two target powers must be distinct, with an even field count')
    if args.image_basis and args.calls not in (None, args.n):
        p.error('The image basis needs one query per field')
    solve(args)
