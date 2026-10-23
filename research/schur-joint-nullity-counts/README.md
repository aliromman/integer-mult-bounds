# Polynomial exact Schur header counts

The existing binary header formula enumerates2^(d^2) compressed matrices to
get their joint rank law. This package replaces that enumeration with an
exact triangular joint-nullity recurrence. It has O(d^4) summands and works
over every finite field. The unchanged one-block header ensemble then has
exact complete-header counts without materializing any compressed matrix.

[PROOF.md](PROOF.md) gives the incidence law, inversion, source-bound count
and full scope. [counts.py](counts.py) uses integer arithmetic only.
[RESULT.json](RESULT.json) records25,291 matrix comparisons,11,527 factor
pairs,382 product targets, five source header comparisons and20,160 literal
complete-header checks. Extension fields F4 and F8 are included. The executed
d32 binary case counts2^1024 possible matrices through561 strata.

```sh
python3.12 -B scripts/research_package.py --package research/schur-joint-nullity-counts
```

Both registered Python3.12 modes print `PASS_SCHUR_JOINT_NULLITY_COUNTS`.
The standard library suffices. Source bytes and the typed receipt are checked.
The theorem uses a supplied prime-power q; the API does not certify primality.

The complete kappa calculation is unchanged. The all-field sharp-sector
criterion is retained as a cited premise. These are one-block counts, not
later-guard or physical admission probabilities. The residual need not commute,
and different matrices with the same nullities need not give the same residual
orbit. No physical header generator or complete consumer cost is supplied.

## Open next questions

Use this exact count to select one source-bound profile with sufficient
one-block mass, then preserve its full residual states in a later-guard
recursion. Start with the two small binary profiles. A nullity-only repeated
transition is not justified by this theorem. Physical selection and emission
need a separate paid contract.

## Credit and coordination

Base7c472563474e60185e6ca101a17622133d047e56, priority7.
[SOURCES.json](SOURCES.json) pins ssafwann's PR4678 binary header theorem and
checker, and the PR5267 all-field sharp-sector statement with its retained
Nick Petalas attribution. The factor counts and complete-header completion
factor are predecessor results. Gaussian incidence inversion is classical.
The new result removes the matrix enumeration and binds the resulting joint
law to the source count, including extension fields.

A concurrent general-H gluing claim by Mihir Srivastava was found before
publication. An independently completed gluing map was shared on the public
fork as corroboration. This focused package makes no priority claim for that
gluing theorem and does not contain its implementation. The counting scope
is complementary. Live PR/Discussion bodies and relevant main notes were
checked; arbitrary historical diffs and unpublished work remain unsearched.

Ali Romman with substantial OpenAI Codex assistance in derivation,
implementation, review and verification. Apache-2.0. Original notices remain
unchanged. No independent human or proof-kernel acceptance is claimed.
