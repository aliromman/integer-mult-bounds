# Paid subdivisions for productive helpers

This is an author-checked checkpoint with **conditional kappa 0.00311224341**.
The seventeen cuts overlap the earlier paid-meeting claim in
[Discussion5600](https://github.com/CrocSwap/integer-mult-bounds/discussions/5600#discussioncomment-19197520).
That investigation has a stronger 23-split candidate. This checkpoint supplies
verification and completion evidence; no separate construction priority is claimed. It improves the complete
claim in [PR5718](https://github.com/CrocSwap/integer-mult-bounds/pull/5718),
0.00311024570, by **0.00000199771** (about 0.0642%). This is a modest gain.
The selected frontier is unchanged.

The new word adds seventeen paid transport subdivisions. One replaces H2 by
two H1 calls. Sixteen replace H6 by H5 and H1. Each change permits a productive
H9 helper to become two independent targets, H8 and H1. Strict concavity gives
a positive signed moment after the extra paid calls. These changes use actual
values at each cut. Input basis changes and output repairs have explicit costs.

The complete bill has 66,654 paid calls, 13,282 productive fields and rank
surplus 760. Its raw crossing is in (0.00311224343152, 0.00311224343153).
The local word has 761 public inputs, nine owned labels, eleven owned lifetimes,
peak two owned fields and 27,502 instructions. No new owned lifetime is needed.

[PROOF.md](PROOF.md) gives the construction, signed invoice and completion
parameters. [PLAN.json](PLAN.json) specifies the seventeen fixed steps;
[WITNESS.json](WITNESS.json) records their exact boundary data.
[RESULT.json](RESULT.json) records the complete replay and price.
[VALIDATION.md](VALIDATION.md) states the checks and their limits.

## Reproduce

Use Python 3.12 with its standard library. From the repository root, run:

```sh
python3.12 -I -B research/p104-paid-helper-subdivisions/verify.py
python3.12 -I -O -B research/p104-paid-helper-subdivisions/verify.py
```

Both commands must print:

```text
PASS_P104_PAID_HELPER_SUBDIVISIONS_KAPPA_00311224341
```

The registered command is:

```sh
python3.12 -B scripts/research_package.py --package research/p104-paid-helper-subdivisions
```

The verifier authenticates the retained source export before it imports source
code. It reconstructs the PR5718 word, applies all seventeen steps, and compares
all generated WORD, WITNESS and RESULT fields. It checks forward rational rows,
separate dense transpose columns, native maps and absolute phases, the complete
physical caller, two exact moment engines, completion, compiler and outer
inequalities. Actual invalid words and invoices must fail. The verifier does
not use Git, a network connection, external packages or a cached native receipt.

The finite local envelope changes from at most 750 to at most 768 public fields
and from fewer than 26,000 to fewer than 28,000 instructions. Every other local
predicate is retained. [LOCAL-EDITS.json](LOCAL-EDITS.json) lists these two edits.
The larger word fits the retained 1,000-slot route allowance and ordinary budget.
[PRICING-EDITS.json](PRICING-EDITS.json) lists all eight price edits, including
the new 16/17 completion parameters. The verifier checks both edit lists exactly.

## Scope and attribution

The five inherited native_all_columns, owned_erasure, compiler, precision_tape
and analytic_recovery assumptions remain. The finite replay does not constitute
formal acceptance of a complete multiplication theorem. The huge uniform word
and complete rounded multiplier are not materialized. Independent review of
this new package is pending. Source review receipts apply only to their stated
source words.

This branch starts from canonical main df0137efd35e11b4445734bf612bdd4c7ebc330c.
The source is flyingviper's PR5718 at a1534ac45401c9a9c80c07aa65c89e00e33594f8,
with substantial OpenAI Codex assistance. It retains Ali Romman's PR5713 and
the full P104 lineage. All source notices remain in the frozen export and
[NOTICE](NOTICE). The new construction, proof, scripts and text are by Ali Romman
with substantial OpenAI Codex assistance, under Apache-2.0.

Coordination and scope changes are in
[Discussion5645](https://github.com/CrocSwap/integer-mult-bounds/discussions/5645),
with the
[paid-subdivision claim](https://github.com/CrocSwap/integer-mult-bounds/discussions/5645#discussioncomment-19197695).
This package changes the paid query graph. It does not repeat fixed-query
boundary searches or the active-mode restoration audit.

## Open next questions

[The structured questions in TASK.json](TASK.json) give first steps and blockers:

- Independently check the seventeen cuts, finite envelope and completion choice.
- Test new intermediate frames or nonempty fresh decoders under a complete bill.
- Discharge a specific retained framework assumption.
