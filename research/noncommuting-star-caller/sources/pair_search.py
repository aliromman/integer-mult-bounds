"""Search the pair constraints forced by a shortest query word.

Ali Romman, with substantial OpenAI Codex assistance. Apache-2.0.
"""
import argparse
import json
from pathlib import Path
import z3


def search(n, seconds, constraint='all'):
    solver = z3.Solver()
    solver.set(timeout=seconds*1000)
    zero = z3.BitVecVal(0, n)
    mask = (1 << n)-1
    us = [(1 << i) if i < n else mask for i in range(n+1)]
    vs = [mask ^ (1 << i) if i < n else mask for i in range(n+1)]

    def parity(x):
        out = z3.Extract(0, 0, x) == 1
        for k in range(1, n):
            out = z3.Xor(out, z3.Extract(k, k, x) == 1)
        return out

    av, bv = {}, {}
    for i in range(1, n):
        av[i], bv[i] = z3.BitVec(f'a{i}', n), z3.BitVec(f'b{i}', n)
        solver.assert_and_track(z3.And(parity(av[i] & vs[i]), parity(bv[i] & us[i]),
                   z3.Not(parity(av[i] & bv[i])),
                   z3.Extract(i, i, av[i]) == 0), f'projector_{i}')
    for transpose in (False, True):
        for first in range(n-1):
            coeff = [z3.BitVecVal(vs[first] if transpose else us[first], n)]
            for last in range(first+2, n+1):
                s = last-1
                out = [zero for _ in range(len(coeff)+1)]
                for d, x in enumerate(coeff):
                    if transpose:
                        px = z3.If(parity(x & us[s]), bv[s], zero) ^ z3.If(parity(x & av[s]), z3.BitVecVal(vs[s], n), zero)
                    else:
                        px = z3.If(parity(x & bv[s]), z3.BitVecVal(us[s], n), zero) ^ z3.If(parity(x & vs[s]), av[s], zero)
                    out[d] = out[d] ^ x ^ px
                    out[d+1] = out[d+1] ^ px
                coeff = [z3.BitVec(f'x{int(transpose)}_{first}_{last}_{d}', n) for d in range(len(out))]
                solver.add([x == y for x, y in zip(coeff, out)])
                target = us[last] if transpose else vs[last]
                conditions = []
                tangent = zero
                for d, x in enumerate(coeff):
                    if d % 2:
                        tangent = tangent ^ x
                    if constraint == 'tangent':
                        continue
                    if constraint == 'constant' and d != 0:
                        continue
                    if constraint == 'leading' and d != len(coeff)-1:
                        continue
                    conditions.append(parity(x & target) == (d == len(coeff)-1))
                if constraint == 'tangent':
                    conditions.append(parity(tangent & target) == bool((len(coeff)-1) % 2))
                solver.assert_and_track(z3.And(conditions), f'pair_{int(transpose)}_{first}_{last}')
    result = solver.check()
    receipt = {'n': n, 'full_calls': n+1, 'solver_result': str(result),
               'timeout_seconds': seconds, 'z3_version': z3.get_version_string(),
               'constraint': constraint,
               'scope': 'positive common-stage binary query identity; endpoint factorization not imposed'}
    if result == z3.sat:
        model = solver.model()
        receipt['interior_projectors'] = {
            str(s): [(model.eval(bv[s]).as_long() if (us[s] >> i) & 1 else 0)
                     ^ (vs[s] if (model.eval(av[s]).as_long() >> i) & 1 else 0)
                     for i in range(n)] for s in range(1, n)}
    elif result == z3.unknown:
        receipt['reason_unknown'] = solver.reason_unknown()
    else:
        receipt['core'] = [str(x) for x in solver.unsat_core()]
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--seconds', type=int, default=60)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--constraint', choices=('all', 'constant', 'leading', 'tangent'), default='all')
    args = parser.parse_args()
    if args.n < 2 or args.n % 2:
        parser.error('n must be positive and even')
    result = search(args.n, args.seconds, args.constraint)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, sort_keys=True))
