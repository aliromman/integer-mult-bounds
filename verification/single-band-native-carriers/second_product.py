"""Verify a dense signed product at both coprime square axes."""
from collections import Counter
from hashlib import shake_256, sha256
import json
from pathlib import Path
import sys
from time import perf_counter

pkg = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(pkg))
from band_product import PackProduct, Product
from check_product import one
from ordinary_carrier import ordinary_swap

def operand(label, size):
    return int.from_bytes(shake_256(label).digest(size), "little") | (1 << (8 * size - 1))

def verify():
    left = -operand(b"single-band-native-carriers left field", 40)
    right = operand(b"single-band-native-carriers right field", 39)
    calls = Counter()
    def route(source, target, bits, before, gap, after):
        calls[bits] += 1
        ordinary_swap(source, target, bits, before, gap, after)
    start = perf_counter()
    previous = sys.argv[2:] == ["--baseline"]
    assert previous or len(sys.argv) == 2
    machine = PackProduct((9, 25), 8, 136) if previous else Product((9, 25), 8, 136, route)
    setup = perf_counter() - start
    row = one(machine, left, right)
    assert row["bits"] == [320, 312] and row["H_calls"] == 1270464
    assert machine.total == 512 and machine.size == 225
    assert previous or calls
    print(json.dumps(dict(row, signed_operands=[str(left), str(right)], true_sides=[9, 25],
                          checked_coordinates=machine.size, padded_sites=machine.total,
                          logical_site_H_pairs=machine.h_calls*machine.total,
                          setup_seconds=round(setup, 4), physical_shuttles_by_width=dict(calls),
                          path=type(machine).__module__,
                          product_sha256=sha256(str(left*right).encode()).hexdigest())),flush=True)
    print("PASS_LINUX_TWO_SQUARE_PRODUCT", flush=True)

if __name__ == "__main__":
    verify()
