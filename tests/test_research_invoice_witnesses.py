"""Replay tariff witnesses by direct paid and credited sums."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import explain_recursive_invoices as explain
sys.path.pop(0)


def invoice(paid=None, targets=None):
    return {'schema': 'integer-mult-recursive-invoice/v1',
            'contract': {'model': 'fixture', 'rank_semantics': 'rank', 'ambient_rank': 10, 'boundary': 'same'},
            'paid': {str(r): n for r, n in (paid or {}).items()},
            'targets': {str(r): n for r, n in (targets or {}).items()}}


class TariffWitnessTests(unittest.TestCase):
    def test_mixed_signs_have_two_explicit_costs(self):
        report = explain.explain_invoices(invoice({2: 2}), invoice({1: 1, 4: 1}))
        self.assertEqual(report['classification'], 'inconclusive')
        witnesses = report['concave_tariff_witnesses']
        self.assertEqual((witnesses['negative_change']['cap'], witnesses['negative_change']['new_minus_old']), (2, -1))
        self.assertEqual((witnesses['positive_change']['cap'], witnesses['positive_change']['new_minus_old']), (4, 1))
        self.assertEqual(witnesses['negative_change']['old_paid_minus_targets'], 4)
        self.assertEqual(witnesses['negative_change']['new_paid_minus_targets'], 3)

    def test_exhaustive_small_signed_invoices_replay_every_witness(self):
        for delta in itertools.product(range(-2, 3), repeat=4):
            old = invoice({r: -d for r, d in enumerate(delta, 1) if d < 0}, {3: 2})
            new = invoice({r: d for r, d in enumerate(delta, 1) if d > 0}, {3: 2})
            report = explain.explain_invoices(old, new)
            for key, witness in report['concave_tariff_witnesses'].items():
                cap = witness['cap']
                expected = sum(d * min(r, cap) for r, d in enumerate(delta, 1))
                self.assertEqual(witness['new_minus_old'], expected)
                self.assertTrue(expected < 0 if key == 'negative_change' else expected > 0)
            self.assertEqual(len(report['concave_tariff_witnesses']),
                             0 if not any(delta) else 2 if report['classification'] == 'inconclusive' else 1)

    def test_surrendered_target_is_debited(self):
        report = explain.explain_invoices(invoice({1: 1, 8: 1}, {9: 1}), invoice())
        witness = report['concave_tariff_witnesses']['negative_change']
        self.assertEqual(witness['old_paid_minus_targets'], 1)
        self.assertEqual(witness['new_paid_minus_targets'], 0)

    def test_cli_binds_inputs_under_isolated_optimized_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp) / name for name in ('old.json', 'new.json')]
            for path, value in zip(paths, (invoice({2: 2}), invoice({1: 1, 4: 1}))):
                path.write_text(json.dumps(value))
            result = subprocess.run([sys.executable, '-I', '-O', '-B', str(ROOT / 'scripts/explain_recursive_invoices.py'),
                                     *map(str, paths)], text=True, capture_output=True, check=True)
            report = json.loads(result.stdout)
            self.assertEqual(report['inputs']['old']['sha256'], hashlib.sha256(paths[0].read_bytes()).hexdigest())
            self.assertEqual(len(report['concave_tariff_witnesses']), 2)

    def test_mismatched_contract_is_not_explained(self):
        old, new = invoice(), invoice(); new['contract']['boundary'] = 'different'
        with self.assertRaisesRegex(ValueError, 'contracts differ'):
            explain.explain_invoices(old, new)


if __name__ == '__main__':
    unittest.main()
