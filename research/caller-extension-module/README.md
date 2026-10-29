# A two-to-four full-call extension test

No positive stage module or new kappa was found. This directory saves a
complete query construction and three exact first-coefficient failures.

The proposed module replaces two disjoint full calls on four old fields
with four full calls on those fields and two new fields. In output, work,
source order, its plus endpoint would be

    diag(z,z,z^2,z^2,z,z).

Its minus endpoint would also send the two old source fields to their
old outputs with coefficient -2 and negate the two new work fields. A
repeatable module of this form could add two Gaussian targets for two
additional full source calls. This is only the reason for the test;
repeatability, physical binding and a complete price have not been proved.

## Exact query construction

Start with four rank-one transvections on two work directions whose product
is -I. Write their nilpotents as u_i v_i, put W=[u_i], V=[v_i], G=VW and
D=I+2lower(G). Then W D^(-1)V=I. The complement Z=D-G has rank two.

Choose two independent columns A of Z. The four-by-four matrix

    E=[D^(-1)V, D^(-1)A]

is invertible. Let B be the last two rows of E^(-1). The six-field queries
use image columns from [B;W;0] and covector rows from [0,V,A]. Direct
chronological multiplication gives exactly the required old transfer and
new minus-identity block. There is no omitted query coefficient: S^2=I
reduces the query test to its two characters.

discovery.py constructs a rational two-parameter family of four-query
work words. It then uses every compatible rank-two projector at each of
the four sites. The full plus and minus first stage coefficients form a
rational linear system in36 chart variables. At (b,t)=(1,1),(1,0),(2,1),
that system is inconsistent. The (1,0) case includes adjacent coincident
work-query directions. No quadratic or higher-order condition can repair
these three fixed frames.

These finite failures do not exclude all work-query parameters, matrix
port coupling, additional stages or a different extension contract.

## Replay and open next questions

Use Python3.12 and the adjacent PR6168 polynomial helper:

    python3.12 -B research/caller-extension-module/discovery.py --b 1 --t 1 --output /tmp/module.json
    python3.12 -B research/caller-extension-module/controls.py
    python3.12 -O -B research/caller-extension-module/controls.py

The controls use literal rational products for the query construction and
a known five-call completion with new target power three. They reject an
omitted query. They also reconstruct all three stated linear failures.

Can a different rational work frame pass the first stage coefficient?
Can a useful module retain a more general old endpoint without paying for
another thin query? A positive construction still needs a proof that it
can be composed, native routes and phases, independent helpers, owned
state, every ordinary operation and the complete compiler/outer bill.

Base: main3f0152a32eeebac0df11946dbd50cabe5a7f7839, priority revision8,
selected kappa .003125804. Canonical Discussion5593, comment19199752.
The Gram and rank-two projector machinery uses PR6168 and the credited
Cinna Davis/P104/native/controller lineage. New module test by Ali Romman
with substantial OpenAI Codex assistance. No subagents. Apache-2.0.
