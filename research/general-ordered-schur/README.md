# Whole ordered Schur normal form for arbitrary H

[PROOF.md](PROOF.md) closes the general-H gluing question after the nilpotent
and H=1 primary results. Every sharp common-identity-pivot deletion with a
fixed compressed H and fixed original ordered sector dimensions has one
ordered residual orbit, over every field. The proof gives exact existence
capacities and an explicit invertible map for both residual operators.

The new step splits the original sector00 factor system by the Fitting
components of I-H. It handles arbitrary factor ranks and cross-component
data. The representative retains nonsplit components and characteristic two.
It does not infer ordered conjugacy from equality of central spectra.

[normal.py](normal.py) implements the map over prime fields. The independent
full-header checker covers1244 complete-domain H,90 transported/padded factor
systems,12 larger mixed-component fixtures and three boundary controls.
[RESULT.json](RESULT.json) retains its deterministic receipt.

```sh
python3.12 -B scripts/research_package.py --package research/general-ordered-schur
```

Python3.12 and its standard library suffice. Both registered modes print
`PASS_GENERAL_ORDERED_SCHUR`. Source bytes and the complete typed receipt
are checked. The finite runs corroborate the all-field written argument.

The complete kappa calculation is unchanged. The sharp-sector criterion is
retained as a cited premise. Nonsharp states, physical headers, transition
weights, later guards, source preparation, precision and supplier admission
remain open. Different H/profile data can give the same residual orbit;
the result does not claim that these are unique orbit labels.

## Open next questions

Compute the fixed-H complete-header weight for arbitrary H, including the
unrestricted sector00 factor count. Start with an exact rank-stratified count
and compare a complete small factor domain. Sum coincident residual orbits
before using those weights as an unconditional transition probability.

## Credit and provenance

Base7c472563474e60185e6ca101a17622133d047e56, priority7. All consumed sources
are pinned in [SOURCES.json](SOURCES.json). The source full-rank factor lemma
belongs to MatShep and codyw912; Nick Petalas supplies the cited sharp-sector
criterion. ssafwann's PR5267 supplies the nilpotent specialization and its
three-copy residual formula. Ali Romman's PR5348 supplies the H=1 primary
removal; PR5287 supplies the earlier central classification. This result
proves gluing at the actual original factors and removes the nilpotence
restriction from the ordered-pair theorem. Classical Fitting decomposition
and basis extension are used. Original source notices remain unchanged.

Live handoffs, all-state PR bodies, Discussion text and relevant main notes
were checked. Sparse notes and arbitrary historical diffs were not fully
searched. Ali Romman with substantial OpenAI Codex assistance in proof,
implementation, review and verification; Apache-2.0. No independent human
or proof-kernel acceptance is claimed.
