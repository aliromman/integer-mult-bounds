"""Rebuild the chart 117 circuit and check its missing rational direction."""
from collections import Counter
from fractions import Fraction
from graphlib import TopologicalSorter
from math import gcd, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_component import canonical, check, expected_columns, factor, frame_span, margin, read


def reconstruct(parent, delta):
    check(canonical(parent) == delta["parent"] == read("RATIONAL-RECEIPT.json")["word"], "parent word")
    forms = {k: dict(row) for k, row in parent["forms"]}
    moves = {k: (a, b) for k, a, b in parent["moves"]}
    for k in delta["retire"]:
        del moves[k], forms[k]
    for k, a, b in delta["add"]:
        check(k not in forms, "new transfer has a fresh label")
        moves[k], forms[k] = (a, b), {}
    for k, changes in delta["coefficients"]:
        for j, coefficient in changes:
            if coefficient is None:
                del forms[k][j]
            else:
                forms[k][j] = coefficient
    ends = list(parent["final"])
    for i, chart in delta["final"]:
        ends[i] = chart
    graph = dict(initial=parent["initial"], final=ends, first_output=parent["first_output"],
                 moves=[[k, *moves[k]] for k in sorted(moves)],
                 forms=[[k, sorted(row.items())] for k, row in sorted(forms.items())])
    check(canonical(graph) == delta["child"], "reconstructed circuit")
    return graph


def replay(graph, layout, delta, paid, coform_size):
    starts, ends = graph["initial"], graph["final"]
    count, first = len(starts), graph["first_output"]
    moves = {k: (a, b) for k, a, b in graph["moves"]}
    forms = dict(graph["forms"])
    check(len(starts) == len(ends) == layout["public"] == 1289, "public fields")
    check(len(moves) == paid and forms.keys() == moves.keys() | set(range(first, first + count)),
          "paid and terminal labels")
    cols, main = expected_columns(layout, count)
    rows = [{} for _ in starts]
    for i, column in enumerate(cols):
        for j, pair in column.items():
            rows[j][i] = Fraction(*pair)
    h = {i: Fraction(s) for i, s in delta["h"]}
    q = {i: Fraction(s) for i, s in delta["q"]}
    check(len(h) == 122 and all(starts[i] == 3 and h[i] for i in h), "source direction")
    check(len(q) == coform_size and all(i not in main and starts[i] == 3 and ends[i] == 127 and q[i]
                                        for i in q), "lost helper coordinates")
    check(sum(v * h.get(i, 0) for i, v in q.items()) == 1, "projection normalization")
    controls = {}
    if "repair" in delta:
        transfer, source, receiver = delta["repair"]["transfer"]
        check(transfer not in forms and transfer >= 0 and (source, receiver) == (3, 117),
              "fresh repair transfer from the chart that holds q")
        controls = {k: Fraction(v) for k, v in delta["repair"]["controls"]}
        check(len(controls) == len(delta["repair"]["controls"]) == 15 and controls.keys() <= forms.keys(),
              "repair coefficients")
    correction = {}
    needed = Counter(j for row in forms.values() for j, _ in row)
    stock = {~i: (chart, {i: 1}, 1) for i, chart in enumerate(starts)}
    dag = TopologicalSorter({k: [j for j, _ in row if j >= 0] for k, row in forms.items()})
    products, grid = 0, 1
    for k in dag.static_order():
        row = forms[k]
        chart = moves[k][0] if k in moves else ends[k - first]
        if controls:
            check(k not in controls or chart in (source, receiver), f"repair signal available at {k}")
            amplitude = controls.get(k, 0) + sum(Fraction(s) * correction[j]
                                                 for j, s in row if j in correction)
            if amplitude and k in moves:
                correction[k] = amplitude
        common = lcm(*(factor(s)[1] * stock[j][2] for j, s in row))
        result = {}
        for j, s in row:
            source, values, divisor = stock[j]
            check(source == chart, f"donor chart at {k}")
            num, den = factor(s)
            scale = num * (common // (den * divisor))
            for i, value in values.items():
                result[i] = result.get(i, 0) + value * scale
            products += len(values)
            needed[j] -= 1
            if not needed[j]:
                del stock[j]
        divisor = gcd(common, *result.values())
        den = common // divisor
        result = {i: v // divisor for i, v in result.items() if v}
        grid = lcm(grid, den)
        if k in moves:
            stock[k] = (moves[k][1], result, den)
            continue
        want = rows[k - first].copy()
        response = sum(v * h.get(i, 0) for i, v in want.items())
        if controls:
            check(amplitude == response, f"repaired row {k - first}")
        if response:
            for i, v in q.items():
                want[i] = want.get(i, 0) - response * v
            want = {i: v for i, v in want.items() if v}
        check(result.keys() == want.keys() and all(result[i] * v.denominator == den * v.numerator
                                                  for i, v in want.items()), f"projected row {k - first}")
    check(products == delta["products"] and str(grid) == delta["grid"], "integer replay census")
    return main


def repair_budget(graph, layout, main, repair=None):
    caller = read("CALLER-CENSUS.json")
    calls = caller["calls"]
    spans = list(map(frame_span, layout["frames"]))

    def grade(a, b):
        check(spans[a] < spans[b] or spans[b] < spans[a], "strict chart transition")
        return abs(len(layout["frames"][a]) - len(layout["frames"][b]))

    paid = Counter({int(k): v for k, v in caller["fixed_caller_paid"].items()})
    targets = Counter({int(k): v for k, v in caller["fixed_root_targets"].items()})
    for _, a, b in graph["moves"]:
        paid[grade(a, b)] += calls
    for i, (a, b) in enumerate(zip(graph["initial"], graph["final"], strict=True)):
        if i in main:
            continue
        if a == b:
            check(i == 1093 and a == 117, "lost return")
            continue
        targets[grade(a, b)] += calls
    lo, hi = margin(paid, targets, "0.00391884648")
    if repair is None:
        check(Fraction("6.2262414658") < lo / calls <= hi / calls < Fraction("6.2262414659"),
              "grade-four retirement repair allowance")
        check(graph["initial"].count(3) == 907 and graph["final"].count(3) == 0
              and sum(a == 3 for _, a, b in graph["moves"]) == 906, "source cut at chart 3")
        check(graph["initial"].count(127) == 0 and graph["final"].count(127) == 901
              and sum(b == 127 for _, a, b in graph["moves"]) == 900, "receiver cut at chart 127")
        return
    check(Fraction("2.2479132437") < lo / calls <= hi / calls < Fraction("2.2479132439"),
          "price available for an independent repair")
    _, a, b = repair["transfer"]
    check(grade(a, b) == 5, "repair grade")
    paid[grade(a, b)] += calls
    lo, hi = margin(paid, targets, "0.00391884648")
    check(Fraction("-2.7206502971") < lo / calls <= hi / calls < Fraction("-2.7206502969"),
          "complete repaired per-call margin")


def main():
    parent, delta = read("WORD.json.gz"), read("Q117-DELTA.json.gz")
    layout = read("LAYOUT.json.gz")
    check(canonical(layout) == read("EVENTS.json")["child_layout_sha256"], "layout")
    graph = reconstruct(parent, delta)
    main_roles = replay(graph, layout, delta, 5696, 305)
    repair_budget(graph, layout, main_roles, delta["repair"])
    print("PASS 1289 projected and repaired rows, one grade-five transfer and both per-call budgets")
    alternate = read("Q117-99044.json.gz")
    check(99044 in alternate["retire"], "retired transfer")
    other = reconstruct(parent, alternate)
    roles = replay(other, layout, alternate, 5695, 303)
    repair_budget(other, layout, roles)
    print("PASS 5695 transfers, 1289 projected rows, both cut deficits and the second repair budget")


if __name__ == "__main__":
    main()
