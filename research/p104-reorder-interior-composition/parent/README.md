# Paid interior bases for the P104 word

Complete conditional kappa is **0.00310815849**, a gain of **0.00000139516**
(about **0.04491%**) over
[PR #5660](https://github.com/CrocSwap/integer-mult-bounds/pull/5660).
This is a small gain. A large improvement and a global optimum remain open.

Four temporary changes of the public register basis remove paid contacts.
The basis changes occur at internal meeting charts. All transports to those
charts are included. The source and output interfaces stay fixed.
The complete local delta is `-7 H1 + 2 H2 + H3`. Across 18 full calls, this
removes 72 paid requests. The complete bill has **66,402 paid calls**,
**12,958 productive targets**, and rank surplus **760**.

[PROOF.md](PROOF.md) states the identity and the complete cost.
[WITNESS.json](WITNESS.json) contains all eight scalar block certificates.
[RESULT.json](RESULT.json) is the full receipt.
[SOURCES.json](SOURCES.json) pins every consumed source.

## Replay

Python 3.12 and the standard library suffice. From the repository root:

```sh
python3.12 -B scripts/research_package.py --package research/p104-paid-interior-bases
```

The hook runs normal and optimized checks. The expected final marker is:

```text
PASS_P104_PAID_INTERIOR_BASES_KAPPA_00310815849
```

The direct commands also work in a source export with no Git or network:

```sh
python3.12 -I -B research/p104-paid-interior-bases/verify.py
python3.12 -I -O -B research/p104-paid-interior-bases/verify.py
```

The checks authenticate the 23-file PR5660 export and its 36-file predecessor
export before source imports. They check every exact block identity, all
565,504 transpose coefficients, 2,918 signed native maps, the physical caller,
both rational moment engines, all compiler costs, and all 58 outer conditions.
See [VALIDATION.md](VALIDATION.md) for the exact checks that ran.

## Scope and limits

The five inherited framework assumptions and the helper endpoint contract
remain. The huge uniform program and complete rounded multiplier are not
materialized. Author checks are not independent scientific acceptance.
The selected frontier is unchanged.

The source is PR5660 at `3b0b4f37d5005f381ed9cfcccebfbe1ab35c71aa`.
The branch starts from canonical main `df0137efd35e11b4445734bf612bdd4c7ebc330c`,
priority revision 7. The complete source is retained, including all notices.

## Open next questions

1. Can an independent reviewer confirm the block identities and full bill?
2. Can a larger temporary basis or a different scalar schedule give a large
   gain after all meeting transports and ordinary work are paid?
3. Can one of the inherited framework assumptions be resolved?

The finite discovery searches do not prove a general obstruction. Further
work must check live Discussions and pin its complete input word.

## Credit

[PR #5660](https://github.com/CrocSwap/integer-mult-bounds/pull/5660) supplies
the complete input and the four-gate contact rewrite. It retains zxv's
[PR #5643](https://github.com/CrocSwap/integer-mult-bounds/pull/5643) A2 change,
the [PR #5644](https://github.com/CrocSwap/integer-mult-bounds/pull/5644) helper
relations, and the [PR #5650](https://github.com/CrocSwap/integer-mult-bounds/pull/5650)
composition. The source notices retain the P104, side-correction, helper,
native, caller, productive-recursion, interval and compiler authors.
The elementary change-of-basis identity is not claimed as a new theorem.

The paid interior construction and this verification package are by
Ali Romman with substantial OpenAI Codex assistance. No subagents.
Apache-2.0. See [NOTICE](NOTICE) for the full retained attribution.
