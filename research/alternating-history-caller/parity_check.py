"""Check the stage-rank proof on literal rational block matrices.

These finite controls support the written argument. They do not replace it.
"""
import json
from pathlib import Path
import sys
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'five-call-projector-rigidity'))
from exact_replay import canonical


def full_factor(n, q, site):
    us, vs, _ = canonical(n)
    u, v = us[site], vs[site]
    pivot = next(i for i in range(n) if v[i])
    drop = next(i for i in range(n) if u[i] and i != pivot)
    h = [s.Rational(int(i == pivot), v[pivot]) for i in range(n)]
    columns = [u, h]+[[int(i == j)-v[j]*h[i] for i in range(n)]
                      for j in range(n) if j not in (pivot, drop)]
    frame = s.kronecker_product(s.Matrix.hstack(*map(s.Matrix, columns)), s.eye(q))
    a = s.Matrix((n-2)*q, q, lambda i,j: (i+2*j+site)%3-1)
    b = s.Matrix(q, (n-2)*q, lambda i,j: (2*i+j+2*site)%3-1)
    left = s.BlockMatrix([[s.eye(q), -b*a], [s.zeros(q), s.eye(q)],
                          [s.zeros((n-2)*q, q), a]]).as_explicit()
    right = s.BlockMatrix([[s.eye(q), s.zeros(q), b],
                           [s.zeros(q), s.eye(q), s.zeros(q, (n-2)*q)]]).as_explicit()
    return frame*left, right*frame.inv(), s.kronecker_product(s.Matrix(u), s.eye(q)), s.kronecker_product(s.Matrix([v]), s.eye(q))


def check(n, q, full_sites, idle):
    size = n*q
    stages = [full_factor(n, q, i) for i in full_sites]
    if idle:
        stages.insert(1, (s.eye(size)[:, :1], s.eye(size)[:1, :], None, None))
    a = s.Matrix.hstack(*(x[0] for x in stages))
    b = s.Matrix.vstack(*(x[1] for x in stages))
    offsets = [0]
    for left, right, _, _ in stages:
        if right*left != s.eye(left.cols): raise ValueError('Projector factor control fails')
        offsets.append(offsets[-1]+left.cols)
    g, dimension = b*a, offsets[-1]
    d = s.eye(dimension)
    for i in range(len(stages)):
        for j in range(i):
            d[offsets[i]:offsets[i+1], offsets[j]:offsets[j+1]] = 2*g[offsets[i]:offsets[i+1], offsets[j]:offsets[j+1]]
    z = d-g
    word = s.eye(size)
    for left, right, _, _ in stages: word = (s.eye(size)-2*left*right)*word
    if word != s.eye(size)-2*a*d.inv()*b: raise ValueError('Chronological factorization fails')
    full = len(full_sites)
    row, col = s.zeros(full*q, dimension), s.zeros(dimension, full*q)
    slot = 0
    for i, (left, right, u, v) in enumerate(stages):
        if u is None: continue
        row[slot*q:(slot+1)*q, offsets[i]:offsets[i+1]] = v*left
        col[offsets[i]:offsets[i+1], slot*q:(slot+1)*q] = right*u
        slot += 1
    expected = s.kronecker_product(s.ones(full)-s.eye(full), s.eye(q))
    if row*z*col != expected: raise ValueError('Query compression fails')
    k = size-(s.eye(size)+word).rank()
    if z.rank() != dimension-k: raise ValueError('Eigenspace rank identity fails')
    if full >= 2 and k > full*q+int(idle): raise ValueError('Stage parity bound fails')
    return {'n': n, 'q': q, 'full_sites': list(full_sites), 'idle_rank': int(idle),
            'minus_one_dimension': k, 'rank_z': z.rank(), 'compressed_rank': expected.rank()}


def verify():
    cases = [check(4,1,(0,1,3,4),False), check(6,1,(0,1,3,5,6),True),
             check(4,2,(0,1,3,4),True), check(8,1,(0,2,5,8),False)]
    exception = check(2,1,(0,),False)
    if exception['minus_one_dimension'] != 2: raise ValueError('The one-full exception was lost')
    return {'cases': cases, 'one_full_exception': exception, 'new_kappa': False}


if __name__ == '__main__':
    result=verify()
    output=Path(__file__).with_name('discovery-results')/'parity-controls.json'
    if '--write' in sys.argv: output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('PASS_STAGE_PARITY_CONTROLS')
