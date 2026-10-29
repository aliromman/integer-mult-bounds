# Exact forward cyclic tariffs with arbitrary workspace

For every field K and every p>=2, independent targets z^d_j in the full
formal algebra K[z]/(z^p-1) need total positive-z query rank at least
sum d_j, where 0<=d_j<p. Coordinate calls attain the bound. Arbitrary
constant linear mixing, feedback, clean fields and dirty workspace do not
reduce it. This includes characteristics which divide p.

[FORWARD-TARIFF.md](FORWARD-TARIFF.md) gives the proof. A finite nilpotent
state representation removes cyclic coefficient aliasing. A block determinant
then bounds the first nonzero determinant degree. This proves, for example,
that four inverse C3 targets still need eight forward rank-one calls with
any finite workspace and any scalar field. It closes the workspace question
in Nick Petalas's [C3 package](../c3-four-inverse-tariff/README.md), while
retaining that result and its credit.

[PROOF.md](PROOF.md) also gives a divided-difference rank budget. At a common
cyclic stage of period p, a target with power multiplicities m_j needs rank
sum m_j min(d_j,p-d_j). Each full call on q ports supplies at most 2q.
This bound permits both stage orientations and arbitrary constant matrices.
The stronger forward theorem charges the residues d_j themselves.

These statements exclude the six-field/eight-call and eight-field/ten-call
period-seven proposals with all main stage ranks equal to 24, at every
coupled width. Both proposals pass determinant congruence. Cyclic wrapping
does not supply their missing stage rank.

For the central-query, full-only contract on the selected PR6030 invoice,
the favorable complex saving is below **0.00511060056** at every length.
This is a ceiling, not an attained improvement. The selected kappa remains
**0.003125804**. The donor's requested large kappa gain has not been achieved.
No native route, caller, compiler, precision or outer admission is claimed.

## Evidence and replay

Use CPython 3.12 and its standard library, from the repository root:

```sh
python3.12 -B scripts/research_package.py --package research/cyclic-stage-rank-budget
python3.12 -I -B research/cyclic-stage-rank-budget/verify.py
python3.12 -I -B research/cyclic-stage-rank-budget/state_check.py
```

The registered commands run both checkers normally and with optimization.
They write only generated receipts to stdout. The universal claims depend
on the written proofs, supported by these exact finite checks:

- Every scalar cyclic power for periods 2 through 19, with up to four
  character partitions per period, over fields with distinct cyclic roots.
- Fifteen literal affine products, their telescoping expansions, and six
  complete coupled target matrices, including the tight open profile.
- 216 state, alias-removal and determinant-degree cases in characteristics
  2, 3, 5 and 7, including cyclic orders divisible by the characteristic.
- Twelve direct attaining mixed-power words with unchanged dirty targets.
- Source, missing-character, false polynomial-rank, determinant-congruence
  and nonfaithful-representation controls; exact rational price intervals.

The two scientific inputs are authenticated before use. The interval helper
comes from the preceding all-length package; its unrelated verification
routine is not run. Full provenance history checks need the relevant Git
objects. Local, exported and hosted replay states are in VALIDATION.md.

## Limits

The tariff assumes formal cyclic equality, or a faithful native action.
Free maps have coefficients constant in the selected stage. A second paid
generator, inverse calls, address maps which do not commute with the stage,
nonlinear instructions and promised input subspaces change the model. The
native application also keeps every helper autonomous. Using a transformed
helper as an extra stage-dependent channel requires a new rank budget. The
inverse-call extension of the divided-difference bound does not make the
stronger forward tariff an inverse-call theorem. A concave higher-rank price
also needs its own physical rank-to-cost relation.

## Open next questions

The balanced twelve-field/fourteen-call profile with eight powers 2 and four
powers 3 meets both rank budgets with equality. Its favorable price can exceed
the selected complex saving. It still needs an exact complete query endpoint.
The separate joint-query-stage-search archive records unsuccessful real and
complex searches and their conditioning; it is not part of this exact replay.

A useful changed route is a fully specified ordinary native map which does
not commute with the stage, or a changed paid generator. First test its full
endpoint and optimistic complete bill. The current BIT cap is also a limit
on any subsequent integer kappa claim. Structured questions are in TASK.json.

## Credit

Nick Petalas's C3 four-inverse and Q8 cyclic cubic packages supply the prior
special cases and register questions. Wentao Li's cyclic 2-power theorem
uses a different unit-cost power alphabet. Their predecessors, including
eumemic, keep their credit. The native stages, central query normalization,
P104 source, helper refinements and compiler retain all source attribution,
including Cinna Davis. The two-full border rediscovered during numerical
search is prior work by Dillon/PR4623 and Yiyang/PR4966.

The state representation, Vandermonde argument and block determinant identity
are standard tools. The new contribution is the general tariff proof and its
caller application. Ali Romman prepared this work with substantial OpenAI
Codex assistance. No subagents, independent human review, proof-kernel
acceptance or worldwide priority claim. New work is Apache-2.0.
