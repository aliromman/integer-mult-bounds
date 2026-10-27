"""Finite binary subspaces for the local frame experiments."""

def span(rows):
    out = {0}
    for row in rows:
        out |= {x ^ row for x in out}
    return frozenset(out)


def basis(space):
    out, current = [], {0}
    for row in sorted(space):
        if row not in current:
            out.append(row)
            current |= {x ^ row for x in current}
    return out


def join(a, b):
    return span(basis(a)+basis(b))


def all_between(lo, hi):
    spaces = {lo}
    level = {lo}
    while level:
        following = set()
        for current in level:
            for row in hi-current:
                new = current | {x ^ row for x in current}
                if new not in spaces:
                    following.add(new)
        spaces |= following
        level = following
    return sorted(spaces, key=lambda s: (len(s), sorted(s)))


