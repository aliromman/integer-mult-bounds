#!/usr/bin/env python3
"""Verify the scalar index rank of each listed six-dimensional chart complement."""
import gzip
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_decoder_triples import (PORTS, WORD_SHA, binary_basis, collect, reduce_binary,
                                    require, run_tape, write_integral)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def check_case(entry, roots, blocked):
    clique = entry['clique']
    require(clique == sorted(set(clique)) and len(clique) == 8, 'eight distinct clique vertices')
    require(all(i in range(PORTS) and i not in blocked for i in clique), 'clique outside the chart')
    require(all(dot(roots[i], roots[j]) for i, j in combinations(clique, 2)), 'nonorthogonal clique')

    def allowed(i):
        require(i in range(PORTS) and i not in blocked and i not in clique, 'valid list vertex')
        colors = [j for j in clique if not dot(roots[i], roots[j])]
        require(0 < len(colors) <= 2, 'one or two remaining basis coordinates')
        return colors if len(colors) == 2 else colors * 2

    def edge(before, after):
        i, a = divmod(before, 2)
        j, b = divmod(after, 2)
        left, right = allowed(i), allowed(j)
        if i == j:
            require(left[0] == left[1] and (a, b) == (1, 0), 'canonical label for a singleton')
            return
        require(dot(roots[i], roots[j]) != 0, 'the two vertices require zero cross pairing')
        common = set(left) & set(right)
        require(common == {left[a]} and left[a] == right[1 - b], 'forced opposite choice')

    forward, reverse = entry['forward'], entry['reverse']
    require(len(forward) > 1 and len(reverse) > 1, 'two nonempty implication paths')
    require(forward[-1] == (forward[0] ^ 1) and reverse[0] == forward[-1]
            and reverse[-1] == forward[0], 'both choices force their negation')
    for path in (forward, reverse):
        for before, after in zip(path, path[1:]):
            edge(before, after)
    return len(forward) + len(reverse) - 2


def verify(path, receipts):
    data = path.read_bytes()
    require(hashlib.sha256(data).hexdigest() == WORD_SHA, 'pinned #7786 word')
    word = json.loads(gzip.decompress(data))
    frames = {int(k): frame for k, frame in word['frames'].items()}
    roots = [write_integral(frames[c]['B'][0]) for c in word['start'][:PORTS]]
    require(len(roots) == 111 and frames[0]['dim'] == 0, 'original line set and empty chart')
    summaries, _ = collect(word, frames, run_tape(word, False), run_tape(word, True))
    for i, root in enumerate(roots):
        unknown = sum(1 << j for j, other in enumerate(roots) if dot(root, other))
        require(not reduce_binary(1 << i, binary_basis(row & unknown for row in summaries)),
                'nine binary summaries recover each input with orthogonal side information')
    expected = {0} | {c for c, f in frames.items() if f['dim'] == 6}
    proofs = {int(c): row for c, row in json.loads(gzip.decompress(receipts.read_bytes())).items()}
    require(len(expected) == 181 and set(proofs) == expected, 'all 180 six-charts and the full graph')
    lengths = []
    for c, entry in proofs.items():
        excluded = {i for i, v in enumerate(roots) if all(not dot(v, u) for u in frames[c]['A'])}
        lengths.append(check_case(entry, roots, excluded))
    return dict(word_sha256=WORD_SHA, listed_six_charts=len(expected) - 1, whole_graph=True,
                scalar_rank_lower=9, maximum_implication_steps=max(lengths))


if __name__ == '__main__':
    root = Path(__file__).resolve().parent
    print(json.dumps(verify(root / 'BIT-WORD.json.gz', root / 'INDEX-OBSTRUCTIONS.json.gz'), indent=2))
