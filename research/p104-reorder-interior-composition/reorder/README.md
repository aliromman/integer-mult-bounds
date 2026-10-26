# One helper start lowered by an exact reorder on PR5660's P104 word: conditional kappa .00310680117

Complete conditional kappa **.00310680117** (claim line, tested at .003106801175). The exact raw crossing is in
(0.00310680119019, 0.00310680119020). That is **+3.784e-8** over the raw crossing 0.0031067633449 of
[PR5660](https://github.com/CrocSwap/integer-mult-bounds/pull/5660) (Ali Romman, claim .00310676333), the newest
P104 head when this was written.

**Conditional on** exactly the assumptions of PR5660 and its inherited P104 line, including free helper endpoints
(PR5561, audited in PR5567). This word changes one helper start frame. The selected record is not changed.
No independent human or proof-kernel acceptance is claimed.

## What changes

It is two re-framed ADDs plus one helper start, found as the optimum of R22's exact reorder+re-frame MIP
(earlier package: the R22 branch on PR5580, +4.09e-6). On PR5627, PR5650 and PR5660 the optimum is the same
move, both on the restricted domain and on the wide (dim gap <= 3) domain.

1. Helper 473 starts on the empty chart 3 (dim 0) instead of chart 594 (dim 1). Its credit rises from H8 to H9.
2. `213 += 1/2*473` now runs on chart 3 instead of 594. To do this it is ordered before the T/T-commuting
   `213 += 2*706` together with its `MOVE 213 3 -> 594`. Nothing else on 213 or 473 is crossed.
3. 473 waits on chart 3 and moves straight to 716. `473 += 100` runs on chart 716 (dim 3) instead of 718 (dim 4),
   so `MOVE 100 716 -> 718` disappears and 100 later moves 716 -> 783.

Per local call: paid MOVE jobs H1 -2, H3 +1 (66474 -> 66456 paid calls in all); helper credits H8 -1, H9 +1;
the rank-mass surplus is unchanged at 760 and there are still 12958 productive fields.

## Check

```
python3 -W error -B research/p104-reframe-r22/verify.py          # ~4 minutes; compares with RESULT.json
python3 -W error -B research/p104-reframe-r22/verify.py --write  # regenerates RESULT.json
```

`verify.py` does the following:

- It pins PR5660's `WORD.json.gz` (SHA256).
- It regenerates `WORD.json.gz` from it with `construct.py`. Every edit replaces one unique event run.
- It checks that only helper 473's start changes, and that the new scalar order is a linear extension of the
  dependence order: T/T and S/S pairs commute, every other pair on a shared register keeps its order, and owned
  registers are exclusive.
- It then runs PR5660's own unmodified stack from its exported sources (`sources.py`):
  - exact replay of all rows;
  - every dense transpose column;
  - every native map and phase;
  - the complete physical caller;
  - PR5660's `pricing.py`, with only the four literals in `LITERALS` changed: claim/tested kappa, previous
    claim, and the paid-call total.
- It runs six rejected controls:
  - old helper start;
  - read left on the old chart;
  - missing helper MOVE;
  - unpaid query;
  - extra target credit;
  - helper read swapped after its write.

The construction applies unchanged to PR5650's word, giving .00310648252 (verified on PR5650's stack
during development).

## Attribution

All retained checkers, words and pricing are by the PR5660 authors and their sources; see `../p104-contact-rewrite/NOTICE`.
The R22 search tools (exact commutation-poset + subspace MIP, `loop104.py`/`mip104.py`) are in the R22 package
on branch `r22-p104-reframe-optimum`.

Joel Pulikkan (GamingPuzzled) with Claude Opus 5.5 (Anthropic) via Claude Code. Apache-2.0.
