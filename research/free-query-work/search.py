"""Search arbitrary full-query planes with initialized, discarded work.

Each call is I+U A(z,s) V, with VU=I_2. All parameters may be complex.
Numerical output is discovery only. No small residual certifies an endpoint.
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


def frames(r, work):
    n = r+work
    us, vs, _ = canonical(r)
    out = []
    for u, v in zip(us, vs):
        u, v = u+[0]*work, v+[0]*work
        pivot = next(i for i in range(n) if v[i])
        drop = next(i for i in range(n) if u[i] and i != pivot)
        h = [int(i == pivot)/v[pivot] for i in range(n)]
        columns = [u, h]+[[int(i == j)-v[j]*h[i] for i in range(n)]
                          for j in range(n) if j not in (pivot, drop)]
        basis = np.array(columns, dtype=complex).T
        out.append(basis)
    return out


def initial(r, work):
    out = []
    for basis in frames(r, work):
        u2 = basis[:, :2]
        v2 = np.linalg.inv(basis)[:2]
        out.extend(u2.ravel())
        out.extend(v2.ravel())
    return np.array(out)


def calculator(r, work, period=0, normalized=False):
    n, calls = r+work, r+1
    local, count = 4*(r+work), 4*(r+work)*(r+1)
    ident = np.eye(n, dtype=complex)
    e01 = np.array([[0, 1], [0, 0]], dtype=complex)
    exponents = np.array([2]*(r-2)+[3, 3])
    samples = period or calls+1
    roots = np.exp(2j*np.pi*np.arange(samples)/samples)
    du = np.zeros((local, n, 2), dtype=complex)
    dv = np.zeros((local, 2, n), dtype=complex)
    du[:2*n] = np.eye(2*n).reshape(2*n, n, 2)
    dv[2*n:] = np.eye(2*n).reshape(2*n, 2, n)

    def calculate(x):
        values = x[:count]+1j*x[count:]
        pairs = [(v[:2*n].reshape(n, 2), v[2*n:].reshape(2, n))
                 for v in values.reshape(calls, local)]
        residuals, derivatives = [], []
        pair_derivatives = []
        for site, (u, v) in enumerate(pairs):
            if normalized:
                inverse = np.linalg.inv(v@u)
                actual = inverse@v
                actual_derivative = inverse@(dv-(dv@u+v@du)@actual)
                pairs[site] = (u, actual)
                pair_derivatives.append(actual_derivative)
                continue
            pair_derivatives.append(dv)
            residuals.append((v@u-np.eye(2)).ravel())
            derivative = np.zeros((count, 2, 2), dtype=complex)
            derivative[site*local:(site+1)*local] = dv@u+v@du
            derivatives.append(derivative.reshape(count, -1).T)
        for z in roots:
            target = np.diag(z**exponents)
            for sign in [1, -1]:
                a = (z-1)*np.eye(2)+z*(sign-1)*e01
                word = ident[:, :r].copy()
                jac = np.zeros((count, n, r), dtype=complex)
                for site, (u, v) in enumerate(pairs):
                    factor = ident+u@a@v
                    derivative = du@a@v+u@a@pair_derivatives[site]
                    jac = factor@jac
                    jac[site*local:(site+1)*local] += derivative@word
                    word = factor@word
                residuals.append((word[:r]-sign*target).ravel())
                derivatives.append(jac[:, :r].reshape(count, -1).T)
        f = np.concatenate(residuals)
        j = np.vstack(derivatives)
        return np.r_[f.real, f.imag], np.block([[j.real, -j.imag], [j.imag, j.real]])
    return calculate


def chart_initial(r, work):
    values = np.zeros((r+1, 4*(r+work)-6), dtype=complex)
    values[:, -2] = 1
    return values.ravel()


def chart_calculator(r, work, period=0):
    """Use a polynomial projector chart with no matrix inverse in the fit."""
    n, calls = r+work, r+1
    local, count = 4*n-6, (4*n-6)*(r+1)
    bases = frames(r, work)
    inverses = [np.linalg.inv(b) for b in bases]
    samples = period or calls+1
    roots = np.exp(2j*np.pi*np.arange(samples)/samples)
    exponents = np.array([2]*(r-2)+[3, 3])
    ident = np.eye(n, dtype=complex)
    da = np.zeros((local, n-2, 2), dtype=complex)
    db = np.zeros((local, 2, n-2), dtype=complex)
    if n > 2:
        da[:2*(n-2)] = np.eye(2*(n-2)).reshape(-1, n-2, 2)
        db[2*(n-2):4*(n-2)] = np.eye(2*(n-2)).reshape(-1, 2, n-2)

    def calculate(parameters):
        values = parameters[:count]+1j*parameters[count:]
        pairs = []
        for site, value in enumerate(values.reshape(calls, local)):
            a = value[:2*(n-2)].reshape(n-2, 2)
            b = value[2*(n-2):4*(n-2)].reshape(2, n-2)
            x, y = value[-2:]
            u = bases[site]@np.vstack((np.eye(2), a))
            v = np.hstack((np.eye(2)-b@a, b))@inverses[site]
            du = bases[site]@np.concatenate((np.zeros((local, 2, 2)), da), axis=1)
            dv = np.concatenate((-db@a-b@da, db), axis=2)@inverses[site]
            nilpotent = np.array([[-x*y, x*x], [-y*y, x*y]])
            dn = np.zeros((local, 2, 2), dtype=complex)
            dn[-2] = [[-y, 2*x], [0, y]]
            dn[-1] = [[-x, 0], [-2*y, x]]
            pairs.append((u, v, du, dv, nilpotent, dn))
        residuals, derivatives = [], []
        for z in roots:
            target = np.diag(z**exponents)
            for sign in (1, -1):
                word = ident[:, :r].copy()
                jac = np.zeros((count, n, r), dtype=complex)
                for site, (u, v, du, dv, nilpotent, dn) in enumerate(pairs):
                    core = (z-1)*np.eye(2)+z*(sign-1)*nilpotent
                    factor = ident+u@core@v
                    derivative = du@core@v+u@core@dv+u@(z*(sign-1)*dn)@v
                    jac = factor@jac
                    jac[site*local:(site+1)*local] += derivative@word
                    word = factor@word
                residuals.append((word[:r]-sign*target).ravel())
                derivatives.append(jac[:, :r].reshape(count, -1).T)
        f, j = np.concatenate(residuals), np.vstack(derivatives)
        return np.r_[f.real, f.imag], np.block([[j.real, -j.imag], [j.imag, j.real]])
    return calculate


def search(args):
    started = time.monotonic()
    if args.chart:
        calculate = chart_calculator(args.r, args.work, args.period)
        seed = chart_initial(args.r, args.work)
    else:
        calculate = calculator(args.r, args.work, args.period, args.normalized)
        seed = initial(args.r, args.work)
    x0 = np.r_[seed.real, seed.imag]
    rng = np.random.default_rng(args.seed)
    test = x0+rng.normal(0, .01, len(x0))
    _, jac = calculate(test)
    for i in rng.choice(len(x0), min(24, len(x0)), replace=False):
        delta = np.zeros(len(x0))
        delta[i] = 1e-6
        fd = (calculate(test+delta)[0]-calculate(test-delta)[0])/2e-6
        if np.max(abs(fd-jac[:, i])) > 2e-5:
            raise ValueError('The analytic Jacobian control fails')
    control = chart_initial(2, 0) if args.chart else initial(2, 0)
    control_calculator = chart_calculator(2, 0, args.period) if args.chart else calculator(2, 0, args.period, args.normalized)
    if np.max(abs(control_calculator(np.r_[control.real, control.imag])[0])) > 1e-12:
        raise ValueError('The exact triangle positive control fails')
    records, best = [], None
    cached_x, cached_result = None, None

    def cached(x):
        nonlocal cached_x, cached_result
        if cached_x is None or not np.array_equal(x, cached_x):
            cached_x, cached_result = x.copy(), calculate(x)
        return cached_result

    for trial in range(args.trials):
        start = x0 if trial == 0 else x0+rng.normal(0, args.scale, len(x0))
        try:
            result = least_squares(lambda x: cached(x)[0], start,
                                   jac=lambda x: cached(x)[1], max_nfev=args.steps,
                                   ftol=1e-11, xtol=1e-11, gtol=1e-11)
            row = {'trial': trial, 'norm': float(np.linalg.norm(result.fun)),
                   'max_residual': float(np.max(abs(result.fun))),
                   'evaluations': result.nfev}
            if best is None or row['norm'] < best['norm']:
                best = row | {'parameters': result.x.tolist()}
        except np.linalg.LinAlgError as error:
            row = {'trial': trial, 'status': 'NUMERICAL_FAILURE_NOT_AN_OBSTRUCTION',
                   'error': str(error)}
        records.append(row)
        print(json.dumps(row), flush=True)
        receipt = {'r': args.r, 'initialized_work': args.work, 'full_calls': args.r+1,
                   'stage_period': args.period, 'normalized_planes': args.normalized,
                   'polynomial_chart': args.chart,
                   'status': 'NUMERICAL_DISCOVERY_ONLY', 'seed': args.seed,
                   'scale': args.scale, 'steps': args.steps, 'records': records,
                   'best': best, 'seconds': time.monotonic()-started,
                   'jacobian_control': 'PASS', 'triangle_control': 'PASS', 'new_kappa': False}
        args.output.write_text(json.dumps(receipt, indent=2)+'\n')
        if row.get('max_residual', float('inf')) < 1e-9:
            break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, choices=(2, 4, 6, 8, 10, 12), default=4)
    parser.add_argument('--work', type=int, choices=(0, 1, 2, 3), default=1)
    parser.add_argument('--trials', type=int, default=8)
    parser.add_argument('--steps', type=int, default=600)
    parser.add_argument('--scale', type=float, default=.15)
    parser.add_argument('--seed', type=int, default=6281)
    parser.add_argument('--period', type=int, default=0)
    parser.add_argument('--normalized', action='store_true')
    parser.add_argument('--chart', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    search(parser.parse_args())
