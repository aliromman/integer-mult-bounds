# Query frames beyond the shortest central tuple

No complete positive caller or new kappa was found. This directory saves
two exact query constructions, their literal checks and failed or unresolved
stage searches. It is exploratory research, not a registered caller package.

The construction in QUERY-FRAMES.md supplies n+2 rank-one queries with
product -I on n fields. The asymmetric construction has three image
directions and only two covector directions in its replacement triple.
Thus it leaves the common nondegenerate symplectic-transvection family.
Neither construction by itself supplies full rank-two stages.

The useful price test is n=12 with 14 full calls and target powers
(2,2,2,2,2,2,2,2,3,3,3,3). Its raw floating price crosses near .00343212
on the selected local profile. Native binding, ordinary fees, compiler and
outer closure are omitted. The retained admitted BIT cap is .0034. No
price in this directory is an achieved exponent.

The 12-field asymmetric instance passes the exact query product and its
first two scalar coefficient screens. Four complete complex numerical
endpoint runs have nonzero residuals; the best norm is about 2.57284.
The six-field instance also failed its eight numerical starts. These
outcomes are failed searches, not proofs of impossibility.

The four-field integer, dyadic and tested asymmetric frames fail the exact
second-coefficient equations. The six- and eight-field rational frames
remain unresolved after that coefficient screen. The six-field scalar
Groebner test timed out. It must not be described as an exclusion.

## Replay

Use Python 3.12 and the adjacent frozen-source proof helper:

    python3.12 -B research/corank-two-query-frame/frame_verify.py
    python3.12 -O -B research/corank-two-query-frame/frame_verify.py

The verifier uses a separate literal rational matrix product, elimination
rank and recursive Pfaffian. It checks ten integer/dyadic frames through
n=16 and all 80 stated asymmetric expansions at n=4,6,8. The changed-product
control destroys the required Pfaffian and rank. The two-triangle frame
and numerical endpoint are positive controls.

The 12-field coefficient command is

    python3.12 -B research/corank-two-query-frame/discovery.py --n 12 --mode expand-6-0 --degree 2 --seconds 180 --output /tmp/n12-chart.json

With NumPy 2.5.4 and SciPy 1.18.1, the full endpoint command is

    OPENBLAS_NUM_THREADS=1 python -B research/corank-two-query-frame/numerical.py --chart /tmp/n12-chart.json --profile '{"2":8,"3":4}' --trials 4 --steps 500 --output /tmp/n12-endpoint.json

The numerical model checks every Laurent coefficient, every projector,
the full stage differential identity and target multiplicities. It checks
24 analytic derivative columns against finite differences and includes
an exact triangle control and changed-target rejection. A numerical zero
would still require exact reconstruction and native/cost verification.

discovery-results retains the charts and receipts. Earlier coefficient
receipts made before chart export was added retain their original coverage.
The equations can be reconstructed from discovery.py; a free affine chart
does not assert that its remaining quadratic equations hold.

## Open next questions

Can a different member of either frame family have a full stage endpoint?
Can the new noncentral S3 contract in the adjacent directory provide one?
A positive answer must include every endpoint coefficient and the complete
native, helper, owned, ordinary, route, compiler and outer bill. Longer
central words are not excluded by the saved failed numerical runs.

Base: canonical main 3f0152a32eeebac0df11946dbd50cabe5a7f7839, priority
revision8, selected kappa .003125804. Canonical Discussion5593, comment19199752.
The query-Gram and projector machinery builds on PR6168 and the retained
Cinna Davis/P104/native/controller lineage. New work by Ali Romman with
substantial OpenAI Codex assistance. No subagents. Apache-2.0.
