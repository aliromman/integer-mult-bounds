"""Replay one binary impulse through two levels of the two-contact caller."""

import gzip
import json
from fractions import Fraction
from functools import cache
from hashlib import sha256
from pathlib import Path

WORD_SHA256 = "d28f0d28d9d87d5fe57b8ef19e12faf98014f0bbef781987ce22daaa82f21cb9"
PORT = (0, 111)
AXES = (tuple(range(8)), (0, *range(8, 15)))
SOURCE = tuple(739 * j * j + 11 * j + 29 for j in range(15))
DESTINATION = (-81429, -779, -3007, -6713, -11897, -18559, -26699, -36317,
               -27762, -46248, -25971, 78786, 123585, 111011, 152854)


def require(test, message):
    if not test:
        raise ValueError(message)


def invert(rows):
    size = len(rows)
    work = [list(row) + [int(i == j) for j in range(size)] for i, row in enumerate(rows)]
    for col in range(size):
        pivot = next(i for i in range(col, size) if work[i][col])
        work[col], work[pivot] = work[pivot], work[col]
        factor = Fraction(1, work[col][col])
        work[col] = [v * factor for v in work[col]]
        for i, row in enumerate(work):
            if i != col and row[col]:
                work[i] = [u - row[col] * v for u, v in zip(row, work[col], strict=True)]
    return [row[size:] for row in work]


def reflect(frame):
    axes = [tuple(map(Fraction, row)) for row in frame["B"] + frame["A"]]
    inverse = invert(axes)
    signed = [[-x if i < frame["dim"] else x for x in row] for i, row in enumerate(axes)]
    return tuple(tuple(sum(a * b for a, b in zip(row, col, strict=True))
                       for col in zip(*signed, strict=True)) for row in inverse)


def caller(word, include_defect):
    @cache
    def matrix(chart):
        return reflect(word["frames"][str(chart)])

    @cache
    def change(address, chart, stage):
        positions = AXES[stage]
        data = [address[i] for i in positions]
        target = list(address)
        for i, row in zip(positions, matrix(chart), strict=True):
            target[i] = sum(a * b for a, b in zip(row, data, strict=True))
        return tuple(target)

    def contact(values, stage):
        """Apply one full outer native with one inner source and its mixed component."""
        size = word["public"] + word["owned_fields"]
        top = [set() for _ in range(size)]
        bottom = [set() for _ in range(size)]
        tags = word["start"].copy() + [None] * word["owned_fields"]
        # Store each address in its current native chart.
        for port, (x, error) in zip(PORT, values, strict=True):
            top[port] = {change(v, tags[port], stage) for v in x}
            bottom[port] = {change(v, tags[port], stage) for v in error}
        for op, i, *args in word["program"]:
            if op == "MOVE":
                before, after = args
                require(tags[i] == before, "move source chart")
                if include_defect and top[i]:
                    # The head defect sends the identity and the relative reflection.
                    moved = {change(change(v, before, stage), after, stage) for v in top[i]}
                    bottom[i] ^= top[i] ^ moved
                tags[i] = after
            elif op == "ADD":
                j, chart = args
                require(tags[i] == tags[j] == chart, "add source charts")
                top[i] ^= top[j]
                bottom[i] ^= bottom[j]
            elif op == "SWAP":
                j, = args
                for state in (top, bottom, tags):
                    state[i], state[j] = state[j], state[i]
            elif op == "INIT":
                chart, = args
                require(tags[i] is None, "owned field already in use")
                tags[i] = chart
            elif op == "CLEAR":
                chart, = args
                require(tags[i] == chart, "owned return chart")
                tags[i] = None
                top[i].clear()
                bottom[i].clear()
            else:
                raise ValueError(f"unknown instruction {op}")
        require(tags[:word["public"]] == word["final"], "main and helper end charts")
        require(all(t is None for t in tags[word["public"]:]), "owned field at return")
        return [({change(v, tags[i], stage) for v in top[i]},
                 {change(v, tags[i], stage) for v in bottom[i]}) for i in PORT]

    def add(x, y):
        return x[0] ^ y[0], x[1] ^ y[1]

    x, y = ({SOURCE}, set()), (set(), set())
    x, y = y, x
    x = add(x, y)
    x, y = contact((x, y), 0)
    y = add(y, x)
    y, x = contact((y, x), 1)
    return y, x


def main():
    packed = Path(__file__).with_name("BIT-WORD.json.gz").read_bytes()
    require(sha256(packed).hexdigest() == WORD_SHA256, "binary word digest")
    word = json.loads(gzip.decompress(packed))
    require(word["frames"][str(word["start"][0])]["B"] == [[1, 0, 0, 0, 0, 0, 0, 0]],
            "coordinate head at port zero")
    desired = {tuple(-v for v in SOURCE)}
    clean, wrong = caller(word, False)
    require(clean == (desired, set()) and wrong == (set(), set()), "diagonal control")
    clean, wrong = caller(word, True)
    require(clean[0] == desired and not wrong[0], "diagonal result with defect")
    require(len(wrong[1]) == 50 and DESTINATION in wrong[1], "mixed output witness")
    print("PASS. The diagonal control gives zero mixed addresses. The nested caller gives 50.")


if __name__ == "__main__":
    main()
