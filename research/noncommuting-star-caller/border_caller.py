"""Verify the two-call residual as an exact Laurent identity."""
import argparse
from fractions import Fraction
import json
from pathlib import Path


def term(coefficient=1, e=0, u=0, s=0):
    return {(e, u, s): Fraction(coefficient)} if coefficient else {}


def add(*parts):
    result = {}
    for part in parts:
        for key, value in part.items():
            result[key] = result.get(key, 0)+value
    return {k: v for k, v in result.items() if v}


def mul(a, b):
    result = {}
    for (e, u, s), x in a.items():
        for (f, v, t), y in b.items():
            key = e+f, u+v, (s+t) % 2
            result[key] = result.get(key, 0)+x*y
    return {k: v for k, v in result.items() if v}


def mm(a, b):
    return [[add(*(mul(a[i][k], b[k][j]) for k in range(2))) for j in range(2)] for i in range(2)]


def run():
    z, one, advance = {}, term(), term(u=1)
    first = [[advance, term(-1, u=1, s=1)], [z, advance]]
    second = [[advance, z], [term(u=1, s=1), advance]]
    before = [[one, z], [z, term(e=-1)]]
    swap = [[z, one], [term(-1), z]]
    after = [[one, z], [z, term(e=1)]]
    actual = mm(after, mm(swap, mm(second, mm(first, before))))
    target = [[term(u=2, s=1), z], [term(-1, e=1, u=2), term(u=2, s=1)]]
    if actual != target:
        raise ValueError('The complete Laurent word fails')
    for residual in ({}, term(1, e=1, u=2)):
        incorrect = [[term(u=2, s=1), z], [residual, term(u=2, s=1)]]
        if actual == incorrect:
            raise ValueError('The residual control was accepted')
    margin = 2*(104*17-1720)
    if margin != 96 or margin-104 != -8:
        raise ValueError('The complete linear repair price fails')
    return {'status': 'PASS_TWO_CALL_BORDER_IDENTITY', 'new_kappa': False,
            'full_calls': 2, 'exact_first_output': 'U^2 S x',
            'second_output': 'U^2 S y - e U^2 x',
            'input_scaling': 'e^-1 on y', 'ordinary_residual': False,
            'counterfactual_linear_margin_without_error_fee': margin,
            'linear_margin_with_one_exact_query_repair': margin-104,
            'scope': 'Exact identity with a retained error term; recursive precision and a profitable repair remain open'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    result = run()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
