"""Test all stated redistributions of the old four stage degrees.

These endpoints change the old full-call contract. They are not direct
replacements for two unchanged full calls.
"""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import importlib.util

spec = importlib.util.spec_from_file_location('module_discovery', Path(__file__).with_name('discovery.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def screen():
    records = []
    diagonal_tests = []
    for b, t in ((F(1), F(1)), (F(1), F(0)), (F(2), F(1))):
        _, ps, ts, target = m.construct(b, t)
        ps = [[[{tuple(v+6 for v in word): c for word, c in entry.items()}
                for entry in row] for row in p] for p in ps]
        h = [[{(i,): F(1)} if i == j else {} for j in range(6)] for i in range(6)]
        matrices = (m.e.ma(m.coefficient(ps, [m.e.eye(6)]*4, 1), h, -1),
                    m.e.ma(m.coefficient(ps, ts, 1), m.e.mm(h, target), -1))
        reduced = m.e.row_reduce([entry for matrix in matrices for row in matrix for entry in row])
        contradiction = any(set(row) == {()} for row in reduced)
        if not contradiction:
            raise ValueError('The unrestricted diagonal-derivative test did not close')
        diagonal_tests.append({'b': str(b), 't': str(t), 'variables': 42,
                               'rank': len(reduced), 'status': 'EXACT_CONSTANT_CONTRADICTION',
                               'scope': 'Every scalar diagonal first stage derivative, before projector quadratics'})
        for old in product(range(3), repeat=4):
            if sum(old) != 4:
                continue
            powers = old[:2]+(2, 2)+old[2:]
            result = m.search(b, t, 90, powers)
            records.append(result)
            print(json.dumps({k: v for k, v in result.items() if k != 'chart'}), flush=True)
    return {'domain': 'Three rational query frames; all19 old-degree profiles in {0,1,2}^4 with sum4',
            'cases': records, 'all_diagonal_derivatives': diagonal_tests, 'new_kappa': False,
            'scope': 'Changed old stage endpoints; no composition, native or price claim'}


if __name__ == '__main__':
    result = screen()
    print(json.dumps(result))
