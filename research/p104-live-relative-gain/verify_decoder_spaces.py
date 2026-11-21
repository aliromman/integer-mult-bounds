#!/usr/bin/env python3
"""Check every binary three-space in the span of the nine decoder returns."""

import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_decoder_triples import (PORTS, WORD_SHA, collect, enlarge, indices,
                    lattice_steps, require, residue, run_tape, write_integral)


def parity(value):
    return bin(value).count('1') % 2


def index_two_coordinates(row, character):
    pivot = (character & -character).bit_length() - 1
    result = [2 * x for x in row]
    result[pivot] = sum(row[j] for j in indices(character))
    return write_integral(result)


def tables(path):
    blob = path.read_bytes()
    require(hashlib.sha256(blob).hexdigest() == WORD_SHA, 'pinned compressed word')
    word = json.loads(gzip.decompress(blob))
    frames = {int(k): v for k, v in word['frames'].items()}
    goals, stocks = collect(word, frames, run_tape(word, False), run_tape(word, True))
    roots = [write_integral(frames[c]['B'][0]) for c in word['start'][:PORTS]]
    for c, rows in stocks.items():
        for row in rows:
            for i in indices(row):
                require(all(sum(a * b for a, b in zip(roots[i], v)) == 0 for v in frames[c]['A']),
                        'every input represented at a bank lies in its chart')
    lattice_rows = roots
    for character in lattice_steps(roots):
        lattice_rows = enlarge(lattice_rows, character)
    lattice_rows = [write_integral(row) for row in lattice_rows]
    reductions = [[residue(row) for row in lattice_rows]]
    for character in range(1, 256):
        reductions.append([residue(index_two_coordinates(row, character)) for row in lattice_rows])
    require(all(all(row) for row in reductions), 'each line has a nonzero reduction')
    evaluation = [sum(((g >> i) & 1) << j for j, g in enumerate(goals)) for i in range(PORTS)]
    affine = []
    for j in range(9):
        options = [a for a in range(512) if all(parity(a & (z | 256)) == (x >> j) & 1
                                                for x, z in zip(evaluation, reductions[0]))]
        require(len(options) == 1, 'each return has unique affine coordinates in the base lattice')
        affine.append(options[0])
    return evaluation, reductions, affine


def main():
    package = Path(__file__).resolve().parent
    signatures, lattices, affine = tables(package / 'BIT-WORD.json.gz')
    rows = [signatures, *lattices, affine]
    data = ''.join(' '.join(map(str, row)) + '\n' for row in rows)
    with tempfile.TemporaryDirectory() as temporary:
        binary = Path(temporary) / 'decoder_spaces'
        subprocess.run(
            ['c++', '-std=c++20', '-O3', str(package / 'decoder_spaces.cpp'), '-o', str(binary)],
            check=True,
        )
        subprocess.run([str(binary)], input=data, text=True, check=True)


if __name__ == '__main__':
    main()
