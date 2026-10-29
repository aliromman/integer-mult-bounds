"""Search both complete stage endpoints with variable query directions.

This is a numerical discovery program. A small residual is not a proof.
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


def triangle(offset, n):
    us, vs, _ = canonical(2)
    factors = []
    for u, v in zip(us, vs):
        u, v = np.array(u, float), np.array(v, float)
        h = v / (v @ v)
        base = np.column_stack([u, h])
        b, c = np.zeros((n, 2)), np.zeros((2, n))
        b[offset:offset+2] = base
        c[:, offset:offset+2] = np.linalg.inv(base)
        factors.append((b, c))
    return factors


def pack(factors):
    return np.concatenate([np.r_[b.ravel(), c.ravel()] for b, c in factors])


def unpack(x, n):
    return [(a[:2*n].reshape(n, 2), a[2*n:].reshape(2, n))
            for a in x.reshape(-1, 4*n)]


def model(n, degrees, length, exact_dual=False, complex_parameters=False, period=0):
    ident, e01 = np.eye(n), np.array([[0., 1.], [0., 0.]])
    # A degree-length polynomial is determined by length+1 distinct points.
    samples = period or length+1
    points = np.exp(2j*np.pi*np.arange(samples)/samples)
    z = np.tile(points, 2)
    s = np.repeat([1, -1], samples)
    kernels = (z-1)[:, None, None]*np.eye(2)+(s-1)[:, None, None]*z[:, None, None]*e01
    target = s[:, None, None]*np.array([np.diag(t**np.array(degrees)) for t in z])
    count = length*4*n
    weight = np.sqrt(len(z))

    def calculate(x, derivative=True):
        if complex_parameters:
            x = x[:count]+1j*x[count:]
        factors = unpack(x, n)
        raw = factors
        if exact_dual:
            factors = [(b, np.linalg.solve(c@b, c)) for b, c in factors]
        matrices = [ident+b@kernels@c for b, c in factors]
        prefixes = [np.broadcast_to(ident, (len(z), n, n))]
        for a in matrices:
            prefixes.append(a@prefixes[-1])
        residues = [((prefixes[-1]-target)/weight).ravel()]
        if not exact_dual:
            residues += [(c@b-np.eye(2)).ravel() for b, c in factors]
        f = np.concatenate(residues)
        if not derivative:
            return np.r_[f.real, f.imag]
        j = np.zeros((len(f), count), complex)
        after = np.broadcast_to(ident, (len(z), n, n))
        total = len(z)*n*n
        for site in range(length-1, -1, -1):
            b, c = factors[site]
            left, right = after@b@kernels, kernels@c@prefixes[site]
            if exact_dual:
                inverse = np.linalg.inv(raw[site][1]@b)
                left_c = left@inverse
                right_c = (ident-b@c)@prefixes[site]
                left_b = left@c
                right_b = c@prefixes[site]
            start, constraint = site*4*n, total+site*4
            for i in range(n):
                for a in range(2):
                    j[:total, start+2*i+a] = (after[:, :, i, None]*right[:, a, None, :]/weight).ravel()
                    if exact_dual:
                        j[:total, start+2*i+a] -= (left_b[:, :, i, None]*right_b[:, a, None, :]/weight).ravel()
                        j[:total, start+2*n+a*n+i] = (left_c[:, :, a, None]*right_c[:, i, None, :]/weight).ravel()
                    else:
                        j[:total, start+2*n+a*n+i] = (left[:, :, a, None]*prefixes[site][:, i, None, :]/weight).ravel()
                        dc, db = np.zeros((2, 2), complex), np.zeros((2, 2), complex)
                        db[:, a] = c[:, i]
                        dc[a, :] = b[i, :]
                        j[constraint:constraint+4, start+2*i+a] = db.ravel()
                        j[constraint:constraint+4, start+2*n+a*n+i] = dc.ravel()
            after = after@matrices[site]
        jac = np.block([[j.real, -j.imag], [j.imag, j.real]]) if complex_parameters else np.vstack([j.real, j.imag])
        return np.r_[f.real, f.imag], jac

    return calculate


def controls(exact_dual=False, complex_parameters=False, period=0):
    n, length = 4, 6
    x = pack(triangle(0, n)+triangle(2, n))
    if complex_parameters:
        x = np.r_[x, np.zeros_like(x)]
    test = model(n, [3]*n, length, exact_dual, complex_parameters, period)
    if np.max(abs(test(x, False))) > 1e-12:
        raise ValueError('The two-triangle positive control fails')
    if np.max(abs(model(n, [2]*n, length, exact_dual, complex_parameters, period)(x, False))) < .1:
        raise ValueError('The changed degree control does not reject')
    rng = np.random.default_rng(290610)
    x += rng.normal(0, .02, len(x))
    _, jac = test(x)
    for col in np.linspace(0, len(x)-1, 24, dtype=int):
        dx = np.zeros_like(x)
        dx[col] = 1e-6
        delta = (test(x+dx, False)-test(x-dx, False))/2e-6
        if np.max(abs(delta-jac[:, col])) > 1e-7:
            raise ValueError('The derivative control fails')
    return {'exact_triangle_float_replay': 'PASS', 'changed_degrees_reject': True,
            'finite_difference_columns': 24}


def initial(n, length, rng, trial, random_frame):
    if random_frame:
        factors = []
        for _ in range(length):
            b = np.linalg.qr(rng.normal(size=(n, 2)))[0]
            c = b.T+rng.normal(0, .2, (2, n))@(np.eye(n)-b@b.T)
            scale = np.diag([2., .5])
            factors.append((b@scale, np.linalg.inv(scale)@c))
        return pack(factors)
    factors = sum((triangle(i, n) for i in range(0, n, 2)), [])
    while len(factors) > length:
        del factors[int(rng.integers(len(factors)))]
    if len(factors) != length:
        raise ValueError('The initializer requires at most three calls per pair')
    x = pack(factors)
    return x if trial == 0 else x+rng.normal(0, .05, len(x))


def run(args):
    check = controls(args.exact_dual, args.complex, args.period)
    degrees = [2]*(args.n-4)+[3]*4
    if args.degrees:
        degrees = list(map(int, args.degrees.split(',')))
        if len(degrees) != args.n:
            raise ValueError('The degree list does not match the field count')
    length = args.n+2
    if args.period:
        required = sum(min(d % args.period, (-d) % args.period) for d in degrees)
        if required > 2*length and not args.calibrate_excluded:
            output = {'status': 'EXCLUDED_BY_CYCLIC_STAGE_RANK', 'n': args.n,
                      'full_calls': length, 'target_degrees': degrees,
                      'stage_period': args.period, 'required_rank': required,
                      'available_rank': 2*length, 'new_kappa': False}
            args.output.write_text(json.dumps(output, indent=2)+'\n')
            print(json.dumps(output))
            return
    calculate = model(args.n, degrees, length, args.exact_dual, args.complex, args.period)
    rng = np.random.default_rng(args.seed)
    records, best, start = [], None, time.monotonic()
    for trial in range(args.trials):
        x = initial(args.n, length, rng, trial, args.random_frame)
        if args.complex:
            x = np.r_[x, rng.normal(0, .05, len(x))]
        if args.resume:
            previous = json.loads(args.resume.read_text())
            if previous['n'] != args.n or previous['full_calls'] != length:
                raise ValueError('The saved search has a different shape')
            if previous.get('complex_parameters', False) != args.complex:
                raise ValueError('The saved coefficient field differs')
            x = np.array(previous['best']['parameters'])
            if trial:
                x += rng.normal(0, .02, len(x))
        bound = args.bound if args.bound else np.inf
        result = least_squares(lambda x: calculate(x, False), x,
                               jac=lambda x: calculate(x)[1], max_nfev=args.steps,
                               bounds=(-bound, bound),
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        row = {'trial': trial, 'norm': float(np.linalg.norm(result.fun)),
               'max_residual': float(max(abs(result.fun))), 'evaluations': result.nfev,
               'max_parameter': float(max(abs(result.x))),
               'elapsed_seconds': time.monotonic()-start}
        records.append(row)
        print(json.dumps(row), flush=True)
        if best is None or row['norm'] < best['norm']:
            best = row | {'parameters': result.x.tolist()}
        output = {'n': args.n, 'full_calls': length, 'target_degrees': degrees,
                  'controls': check, 'records': records, 'best': best,
                  'seed': args.seed, 'steps': args.steps,
                  'random_frame': args.random_frame, 'parameter_bound': args.bound,
                  'exact_dual': args.exact_dual, 'complex_parameters': args.complex,
                  'stage_period': args.period,
                  'scope': 'Scalar B and C; complete sampled polynomial endpoints and C B = I',
                  'status': 'NUMERICAL_DISCOVERY_ONLY', 'new_kappa': False}
        args.output.write_text(json.dumps(output, indent=2)+'\n')
        if row['max_residual'] < 1e-10:
            break


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=6)
    parser.add_argument('--trials', type=int, default=4)
    parser.add_argument('--steps', type=int, default=300)
    parser.add_argument('--seed', type=int, default=290610)
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--random-frame', action='store_true')
    parser.add_argument('--bound', type=float, default=0.)
    parser.add_argument('--exact-dual', action='store_true')
    parser.add_argument('--complex', action='store_true')
    parser.add_argument('--period', type=int, default=0)
    parser.add_argument('--degrees')
    parser.add_argument('--calibrate-excluded', action='store_true',
                        help='Replay a historical negative domain despite its exact rank cut')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.n < 4 or args.n % 2:
        parser.error('Use an even dimension of at least four')
    run(args)
