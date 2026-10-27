"""Search real polynomial frames for a three-block advancing caller.

All irreducible complex representations of (C2)^3 semidirect C3 are
included. The result is numerical discovery evidence only. The constant
query frame is the credited skew simplex in rational_cyclic_search.py.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares
from rational_cyclic_search import frame


def make_bases(n):
    queries, _ = frame(n)
    bases = []
    for query in queries:
        query = np.rint(query)
        pivot = int(np.argmax(np.sum(abs(query)**2, axis=0)))
        u = query[:, pivot]
        active = int(np.argmax(abs(u)))
        v = query[active]/u[active]
        j = int(np.argmax(abs(v)))
        a = np.eye(n)[:, j]/v[j]
        cols = [u, a]
        for k in range(n):
            candidate = np.eye(n)[:, k]-v[k]*a
            if np.linalg.matrix_rank(np.array(cols+[candidate]).T) > len(cols):
                cols.append(candidate)
        basis = np.array(cols).T
        bases.append((basis, np.linalg.inv(basis)))
    return bases


def representations():
    values = []
    for sign in (1, -1):
        for j in range(3):
            r = np.array([[np.exp(2j*np.pi*j/3)]])
            values.append((r, sign*r))
    r = np.roll(np.eye(3), 1, axis=0).astype(complex)
    for sign in (1, -1):
        values.append((r, r @ (sign*np.diag([-1, 1, 1]))))
    return values


def solve(args):
    n, calls = args.n, args.n+1
    bases = make_bases(n)
    reps = representations()
    # The constant term of c commutes with U and cancels from the stage.
    # Fix it to zero. This removes an exact freedom from the coordinates.
    width = 6*(n-2)+2
    prepared = []
    for r, u in reps:
        d = len(r)
        ident = np.eye(d)
        rp = np.array([ident, r, r @ r])
        prepared.append((d, ident, rp, u, np.linalg.matrix_power(u, 2), np.linalg.matrix_power(u, 3)))

    def endpoints(x):
        rows = x.reshape(calls, width)
        output = []
        for d, ident, rp, u, u2, u3 in prepared:
            identity = np.eye(n*d, dtype=complex)
            a, b = identity.copy(), identity.copy()
            for row, (basis, inverse) in zip(rows, bases):
                c = np.einsum('c,cij->ij', row[:2], rp[1:])
                xx = np.einsum('kc,cij->kij', row[2:3*n-4].reshape(n-2, 3), rp)
                yy = np.einsum('kc,cij->kij', row[3*n-4:].reshape(n-2, 3), rp)
                ac = [c, ident]+list(xx)
                bc = [ident, -c-sum((y @ x for x, y in zip(xx, yy)), np.zeros((d, d), complex))]+list(yy)
                av = [sum((basis[i, j]*ac[j] for j in range(n)), np.zeros((d, d), complex)) for i in range(n)]
                bv = [sum((bc[j]*inverse[j, i] for j in range(n)), np.zeros((d, d), complex)) for i in range(n)]
                u0, v1 = basis[:, 0], inverse[1]
                delta = u-ident
                stage = identity+np.block([[u0[i]*delta @ bv[j]+av[i] @ delta*v1[j]
                                           for j in range(n)] for i in range(n)])
                query = np.kron(np.outer(u0, v1), u)
                a, b = stage @ a, (stage-2*query) @ b
            output.append((a, b))
        return output

    def residual(x):
        outputs = endpoints(x)
        if args.query_only:
            values = [a+b for a, b in outputs]
            return np.concatenate([part.reshape(-1) for v in values for part in (v.real, v.imag)])
        e = (np.eye(n)-outputs[3][0])/2
        values = [e @ e-e]
        for (a, b), (_, _, _, _, u2, u3) in zip(outputs, prepared):
            target = np.kron(np.eye(n)-e, u2)+np.kron(e, u3)
            values.extend((a-target, b+target))
        residuals = np.concatenate([part.reshape(-1) for v in values for part in (v.real, v.imag)])
        return np.append(residuals, np.trace(e).real-2)

    rng = np.random.default_rng(args.seed)
    best, records = None, []
    start = time.monotonic()
    for trial in range(args.trials):
        x = np.zeros((calls, width)) if trial == 0 else rng.normal(0, args.scale, (calls, width))
        answer = least_squares(residual, x.reshape(-1), max_nfev=args.steps,
                               ftol=1e-11, xtol=1e-11, gtol=1e-11)
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
              'representations': 'six scalar characters and two three-dimensional representations',
              'bases': [basis.tolist() for basis, _ in bases]}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'seconds': result['seconds'], 'status': result['status']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=4)
    p.add_argument('--trials', type=int, default=8)
    p.add_argument('--steps', type=int, default=300)
    p.add_argument('--scale', type=float, default=.2)
    p.add_argument('--query-only', action='store_true', help='Do not require a Gaussian main-stage endpoint')
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 2 or args.n % 2:
        p.error('Use a positive even dimension')
    solve(args)
