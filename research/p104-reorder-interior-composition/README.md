# P104 reordering and paid-interior composition

This package composes Joel Pulikkan's [PR5666](https://github.com/CrocSwap/integer-mult-bounds/pull/5666)
with [PR5678](https://github.com/CrocSwap/integer-mult-bounds/pull/5678).
Both use PR5660 as their source. PR5678 does not consume PR5666 or claim
that their gains add. This package emits the combined physical word and
checks its complete price.

The target conditional kappa is **0.00310819637**, up **0.00000003788**
over PR5678. This is a small gain. It is not a large improvement or a
global optimum. The source assumptions and helper endpoint contract remain.

Helper473 starts at chart3 instead of594. Its two moved contacts retain
the four paid interior bases. The complete change from PR5678 is
`Pnew-Pold=-36H1+18H3` and `Tnew-Told=-18H8+18H9`.
There are 66,384 paid calls, 12,958 productive targets and rank surplus760.
All 743 independent local inputs, nine owned fields and eleven owned
lifetimes remain. The peak owned occupancy is two.

The continuation check covers all 30,096 eligible common-frame cut
contexts and every helper boundary group. It permits a distinct physical
cut and boundary role. Every prefix/suffix helper pairing is zero.
Thus this exact word has no further independent split under that contract.
This does not exclude a changed word, paid native frame conversions,
a nonlocal cut interface or another caller.

## Reproduction

```sh
python3.12 -B scripts/research_package.py --package research/p104-reorder-interior-composition
```

The manifest registers normal and explicit optimized commands.
Expected: `PASS_P104_REORDER_INTERIOR_KAPPA_00310819637`.
Python3.12 and the standard library suffice. Native verification can take
several minutes. All source files are authenticated before retained code
is imported. The two nested source exports also authenticate their complete
inventories. The verifier reproduces the emitted word and result exactly.

The check includes every source edit, scalar dependency order, 565,504
independent transpose entries, all signed native maps and phases, the
complete physical caller, two rational moment engines, compiler costs and
all 58 outer conditions. Wrong source starts, lost endpoint changes,
missing contacts, unpaid requests and extra targets reject.

[PROOF.md](PROOF.md), [WORD.json.gz](WORD.json.gz), [RESULT.json](RESULT.json),
[INPUTS.json](INPUTS.json) and [SOURCES.json](SOURCES.json) contain the proof,
word, complete evidence and source pins. [PRICING-EDITS.json](PRICING-EDITS.json)
lists every change to the retained pricing code.

## Limits and credit

The five inherited framework assumptions and free helper endpoint contract
remain. The enormous uniform word and rounded multiplier are not executed.
These are author checks, not independent human acceptance or a proof-kernel
theorem. The selected frontier is unchanged. Full-history ancestry checking
is unavailable in this shallow checkout. Hosted CI is a separate status.

Credit Joel Pulikkan (GamingPuzzled), with Claude Opus 5.5 assistance, for
PR5666's reordering and endpoint change. Credit Ali Romman, with OpenAI
Codex assistance, for PR5678's paid interior bases and PR5660's complete
input. Credit zxv PR5643/5650, PR5644, PR5627, PR5599 and every P104,
Clifford, helper, controller, productive-recursion, interval and compiler
author in the retained notices. The new work is this exact composition,
its full cut check and complete conditional admission. Ali Romman with
substantial OpenAI Codex assistance; no subagents. Apache-2.0.

## Open next questions

- Change the scalar word to create a new independent helper cut. A larger
  census of this same common-frame domain cannot supply another split.
- Integrate a compatible native-frame or coded-caller construction, with
  every additional cost. Check active scopes in Discussions5645,5600 and5125.

Exact first steps and limitations are in [TASK.json](TASK.json).
