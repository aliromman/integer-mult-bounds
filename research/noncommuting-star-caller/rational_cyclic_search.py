"""Search rational symplectic full projectors with a finite-order stage.

This is a numerical search. A small residual is not an exact endpoint.
The skew simplex follows from the minimum transvection Gram identity.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares


def frame(n, kind='simplex'):
    size = n+(1 if kind == 'simplex' else 2)
    if kind == 'simplex':
        gram = np.array([[0 if i == j else (1 if i > j else -1)
                          for j in range(size)] for i in range(size)], dtype=float)
    else:
        if kind == 'polygon':
            complement = np.cos(np.pi*(np.arange(size)[:, None]-np.arange(size))/size)
        else:
            angles = [i % 4 for i in range(size-1)]
            delta = sum((-1 if i % 2 == 0 else 1)*a for i, a in enumerate(angles))
            angles.append((2-delta) % 4)
            directions = np.array([(1, 0), (1, 1), (0, 1), (-1, 1)], dtype=float)[angles]
            complement = directions @ directions.T / np.sum(directions**2, axis=1)[:, None]
        gram = np.where(np.arange(size)[:, None] > np.arange(size), complement, -complement)
        np.fill_diagonal(gram, 0)
    if np.linalg.matrix_rank(gram, tol=1e-9) != n:
        raise ValueError('Wrong query Gram rank')
    selected = []
    for j in range(size):
        if np.linalg.matrix_rank(gram[:, selected+[j]], tol=1e-9) > len(selected):
            selected.append(j)
    factor = gram[:, selected]
    coord = np.linalg.lstsq(factor, gram, rcond=None)[0]
    columns = [coord[:, j] for j in range(size)]
    rows = factor
    queries = [np.outer(u, v) for u, v in zip(columns, rows)]
    value = np.eye(n)
    for query in queries:
        value = (np.eye(n)-2*query) @ value
    if not np.allclose(value, -np.eye(n), atol=1e-9):
        raise ValueError('The signed query frame fails')
    affine = []
    for i in range(size):
        choices = [np.outer(columns[i], rows[j])/gram[j, i]+np.outer(columns[j], rows[i])/gram[i, j]
                   for j in range(size) if abs(gram[j, i]) > 1e-9]
        # The affine space has dimension n-2. Choose independent differences.
        directions = []
        for p in choices[1:]:
            delta = p-choices[0]
            if np.linalg.matrix_rank(np.array(directions+[delta]).reshape(len(directions)+1, -1), tol=1e-9) > len(directions):
                directions.append(delta)
        if kind != 'dyadic' and len(directions) != n-2:
            raise ValueError('Wrong symplectic projector dimension')
        affine.append((choices[0], np.array(directions)))
    return queries, affine


def search(args):
    n = args.n
    queries, affine = frame(n, args.frame)
    ident = np.eye(n)
    sites = len(queries)
    signs = [-1 if args.inverse_mask >> i & 1 else 1 for i in range(sites)]
    bases = []
    for q in queries:
        pivot = int(np.argmax(np.sum(abs(q)**2, axis=0)))
        u = q[:, pivot]
        active = int(np.argmax(abs(u)))
        v = q[active]/u[active]
        j = int(np.argmax(abs(v)))
        a = ident[:, j]/v[j]
        cols = [u, a]
        for k in range(n):
            candidate = ident[:, k]-v[k]*a
            if np.linalg.matrix_rank(np.array(cols+[candidate]).T, tol=1e-9) > len(cols):
                cols.append(candidate)
        basis = np.array(cols).T
        bases.append((basis, np.linalg.inv(basis)))
    width = (n-2)*(2 if args.general else 1)
    tangent_origin, tangent_basis = None, None
    if args.formal_tangent:
        tail, tails = ident.copy(), []
        for q in reversed(queries):
            tails.append(tail.copy())
            tail = tail @ (ident-2*q)
        tails.reverse()
        constant = sum(s*(t @ p @ np.linalg.inv(t)-p)
                       for s, t, (p, _) in zip(signs, tails, affine))
        columns = [s*(t @ d @ np.linalg.inv(t)-d).ravel()
                   for s, t, (_, ds) in zip(signs, tails, affine) for d in ds]
        matrix = np.array(columns).T
        _, singular, right = np.linalg.svd(matrix, full_matrices=True)
        rank = int(sum(singular > 1e-8))
        tangent_origin = np.linalg.lstsq(matrix, -constant.ravel(), rcond=1e-8)[0]
        tangent_basis = right[rank:].T
        if np.max(abs(matrix @ tangent_origin+constant.ravel())) > 1e-8:
            raise ValueError('The first mixed equation has no numerical solution')
    def projectors(x):
        if args.complex:
            real, imag = np.split(x, 2)
            x = real+1j*imag
        if tangent_origin is not None:
            x = tangent_origin+tangent_basis @ x
        if args.general:
            result = []
            for row, (basis, inverse) in zip(x.reshape(sites, width), bases):
                a, b = row[:n-2], row[n-2:]
                p = np.zeros((n, n), dtype=complex if args.complex else float)
                p[0, 0] = p[1, 1] = 1
                p[2:, 1], p[0, 2:], p[0, 1] = a, b, -b @ a
                result.append(basis @ p @ inverse)
            return result
        return [p+np.einsum('k,kij->ij', row, ds)
                for row, (p, ds) in zip(x.reshape(sites, n-2), affine)]
    def endpoints(x, value):
        a = ident.astype(complex)
        b = a.copy()
        for p, q, sign in zip(projectors(x), queries, signs):
            stage = ident+(value**sign-1)*p
            a, b = stage @ a, stage @ (ident-2*q) @ b
        return a, b
    def residual(x):
        if args.balanced:
            probe = np.exp(2j*np.pi/args.period)
            ap, _ = endpoints(x, probe)
            lo, hi = args.lower_degree, args.upper_degree
            e = ((sum(s*p for s, p in zip(signs, projectors(x)))-lo*ident)/(hi-lo)
                 if args.formal_tangent else (ap-probe**lo*ident)/(probe**hi-probe**lo))
            values = [e @ e-e]
            if not args.complex:
                values.append(1j*e.imag)
            for j in range(1, args.period if args.complex else args.period//2+1):
                d = np.exp(2j*np.pi*j/args.period)
                a, b = endpoints(x, d)
                target = d**lo*(ident-e)+d**hi*e
                values.extend((a-target, b+target))
            packed = np.concatenate([part.reshape(-1) for a in values for part in (a.real, a.imag)])
            target_rank = args.upper_rank if args.upper_rank is not None else 2*(sites-n)
            return np.append(packed, np.trace(e).real-target_rank)
        values = [sum(endpoints(x, d)) for d in ((-1, 1j, -1j) if args.complex else (-1, 1j))]
        return np.concatenate([part.reshape(-1) for a in values for part in (a.real, a.imag)])
    rng = np.random.default_rng(args.seed)
    best = None
    start = time.monotonic()
    records = []
    for trial in range(args.trials):
        size = (tangent_basis.shape[1] if args.formal_tangent else sites*width)*(2 if args.complex else 1)
        x = np.zeros(size) if trial == 0 else rng.normal(.5, args.scale, size)
        answer = least_squares(residual, x, max_nfev=args.steps, ftol=1e-12,
                               xtol=1e-12, gtol=1e-12)
        row = {'trial': trial, 'residual_norm': float(np.linalg.norm(answer.fun)),
               'max_residual': float(np.max(np.abs(answer.fun))),
               'evaluations': answer.nfev, 'parameters': answer.x.tolist()}
        records.append({k: v for k, v in row.items() if k != 'parameters'})
        if best is None or row['residual_norm'] < best['residual_norm']:
            best = row
            print(json.dumps({k: v for k, v in row.items() if k != 'parameters'}), flush=True)
        if best['max_residual'] < 1e-9:
            break
    result = {'status': 'NUMERICAL_DISCOVERY_ONLY', 'n': n, 'frame': args.frame, 'general_projectors': args.general,
              'period': args.period, 'balanced': args.balanced,
              'formal_tangent': args.formal_tangent,
              'inverse_mask': args.inverse_mask,
              'target_degrees': [args.lower_degree, args.upper_degree],
              'upper_target_rank': args.upper_rank,
              'complex_parameters': args.complex, 'trials': len(records),
              'seed': args.seed, 'seconds': time.monotonic()-start,
              'best': best, 'records': records,
              'query_matrices': [q.tolist() for q in queries],
              'affine_projectors': [[p.tolist(), ds.tolist()] for p, ds in affine]}
    if args.general:
        result['bases'] = [basis.tolist() for basis, _ in bases]
    if args.formal_tangent:
        result['tangent_origin'] = tangent_origin.tolist()
        result['tangent_basis'] = tangent_basis.tolist()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('best', 'records', 'query_matrices', 'affine_projectors', 'bases', 'tangent_origin', 'tangent_basis')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--n', type=int, default=4)
    p.add_argument('--trials', type=int, default=30)
    p.add_argument('--steps', type=int, default=300)
    p.add_argument('--scale', type=float, default=.5)
    p.add_argument('--general', action='store_true', help='Use all compatible rank-two projectors')
    p.add_argument('--complex', action='store_true', help='Allow complex projector parameters')
    p.add_argument('--period', type=int, default=4)
    p.add_argument('--balanced', action='store_true', help='Require the balanced D^2/D^3 target modulo the stage order')
    p.add_argument('--inverse-mask', type=int, default=0)
    p.add_argument('--lower-degree', type=int, default=2)
    p.add_argument('--upper-degree', type=int, default=3)
    p.add_argument('--upper-rank', type=int)
    p.add_argument('--formal-tangent', action='store_true', help='Eliminate the first mixed equation for a formal symplectic target')
    p.add_argument('--frame', choices=('simplex', 'polygon', 'dyadic'), default='simplex')
    p.add_argument('--seed', type=int, default=5786)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.n < 4 or args.n % 2:
        p.error('The dimension must be even and at least four')
    if args.frame == 'dyadic' and not args.general:
        p.error('Use general projectors for the nonsymmetric dyadic Gram frame')
    if args.period < 4:
        p.error('Use a period of at least four')
    if args.lower_degree >= args.upper_degree or (args.upper_degree-args.lower_degree) % args.period == 0:
        p.error('Use distinct target powers at the primitive probe')
    if args.inverse_mask < 0 or args.inverse_mask >= 1 << (args.n+2):
        p.error('The inverse mask is out of range')
    if args.formal_tangent and (args.general or not args.balanced or args.frame != 'simplex' or args.period <= args.n+1):
        p.error('The tangent chart needs a balanced symplectic simplex and more roots than calls')
    search(args)
