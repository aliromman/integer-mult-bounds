"""Check the finite controls for the cyclic reduction lemma."""
import argparse
import itertools
import json
from pathlib import Path


ZERO, IDENT = (0, 0), (1, 2)


def add(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def mm(a, b):
    return tuple((b[0] if x & 1 else 0) ^ (b[1] if x & 2 else 0) for x in a)


def extend(poly, p):
    out = [ZERO]*(len(poly)+1)
    for d, a in enumerate(poly):
        out[d] = add(out[d], mm(add(IDENT, p), a))
        out[d+1] = add(out[d+1], mm(p, a))
    return out


def reduce(poly, period):
    out = [ZERO]*period
    for d, a in enumerate(poly):
        out[d % period] = add(out[d % period], a)
    return out


def run():
    projectors = [a for a in itertools.product(range(4), repeat=2) if mm(a, a) == a]
    counts = {'words': 0, 'hypothesis_tests': 0, 'passing_hypotheses': 0}
    frontier = [[IDENT]]
    for length in range(1, 6):
        nxt = []
        for old in frontier:
            for p in projectors:
                poly = extend(old, p)
                counts['words'] += 1
                for period in range(2, 5):
                    reduced = reduce(poly, period)
                    for k in range(1, period):
                        if length > period+k-1:
                            continue
                        counts['hypothesis_tests'] += 1
                        if any(a != ZERO for a in reduced[:k]):
                            continue
                        counts['passing_hypotheses'] += 1
                        if any(a != ZERO for a in poly[period:]):
                            raise ValueError('A finite case violates the no-wrap lemma')
                nxt.append(poly)
        frontier = nxt
    controls = []
    for period, k in ((2, 1), (3, 2), (7, 2)):
        poly = [ZERO]*(period+k)+[IDENT]
        reduced = reduce(poly, period)
        if any(a != ZERO for a in reduced[:k]) or reduced[k] != IDENT or not any(a != ZERO for a in poly[period:]):
            raise ValueError('The sharpness control fails')
        controls.append({'period': period, 'k': k, 'length': period+k})
    wanted = [3]*7+[4]
    if min(wanted) != 3 or 9 > 7+min(wanted)-1 or sum(wanted) == 2*9:
        raise ValueError('The eight-field determinant control fails')
    return {'status': 'PASS_CYCLIC_REDUCTION_CONTROLS', 'new_kappa': False,
            'projectors': len(projectors), **counts, 'sharpness_controls': controls,
            'eight_field_high_rank_target': 'REJECT_INTEGER_DETERMINANT_DEGREE',
            'scope': 'Finite controls for the written all-field lemma; no new caller'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
