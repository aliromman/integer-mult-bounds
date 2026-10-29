"""Price hypothetical corank-two endpoints on the pinned local profile.

The root is a floating estimate. The two fixed-exponent margins use
rational intervals. No endpoint or omitted cost is supplied by this file.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / 'alternating-caller-closure'
spec = importlib.util.spec_from_file_location('selected_price_intervals', ROOT/'price_check.py')
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def screen():
    raw = (ROOT/'sources/selected6030-result.json').read_bytes()
    if hashlib.sha256(raw).hexdigest() != '84fcbe2a0b64307b5703c363408042bdd91555b0419509aa50daa43b3b1924ad':
        raise ValueError('The selected source differs')
    local = json.loads(raw)['local']
    paid = {int(k): v for k, v in local['paid'].items()}
    helpers = {int(k): v for k, v in local['helpers'].items()}

    def source_float(k):
        return sum(v*r**(1-k) for r, v in paid.items())-sum(v*r**(1-k) for r, v in helpers.items())

    records = []
    for n in (6, 8, 10, 12, 16):
        full = n+2
        targets = {17: n-4, 25: 4}

        def margin_float(k):
            return 104*sum(v*r**(1-k) for r, v in targets.items())-full*source_float(k)

        lo, hi = 0., .02
        if not margin_float(lo) > 0 > margin_float(hi):
            raise ValueError('The exploratory root bracket fails')
        for _ in range(70):
            mid = (lo+hi)/2
            if margin_float(mid) > 0:
                lo = mid
            else:
                hi = mid
        intervals = {}
        for kappa in (F('0.003125804'), F('0.0034')):
            prices = {r: p.power_bounds(r, kappa) for r in set(paid)|set(helpers)|set(targets)}
            low = 104*sum(v*prices[r][0] for r, v in targets.items()) \
                + full*sum(v*prices[r][0] for r, v in helpers.items()) \
                - full*sum(v*prices[r][1] for r, v in paid.items())
            high = 104*sum(v*prices[r][1] for r, v in targets.items()) \
                + full*sum(v*prices[r][1] for r, v in helpers.items()) \
                - full*sum(v*prices[r][0] for r, v in paid.items())
            intervals[str(kappa)] = p.outward((low, high))
        records.append({'n': n, 'full_calls': full, 'target_multiplicities_per_port': targets,
                        'raw_root_estimate': (lo+hi)/2, 'fixed_exponent_margin_intervals': intervals,
                        'unpriced': ['native source and target binding', 'ordinary maps and routes',
                                     'compiler and outer closure'], 'new_kappa': False})
    return records


if __name__ == '__main__':
    print(json.dumps(screen(), indent=2))
