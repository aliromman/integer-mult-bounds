# A stage eigenvalue bound for shortest query words

This is a written proof under the stated contract. The local controls are
in parity_check.py. It is not a new caller or a kappa increase.

Let the base field have characteristic zero. There are nq main fields,
n+1 rank-q square-zero central queries, and no initialized or discarded
fields. The complete query word at stage value z=1 is -I. There are F full
calls, each with a compatible rank-2q stage projector. Other queries are
thin. Allow extra stage-only maps, with total projector rank R. Ordinary
constant invertible changes of frame are absorbed into the projectors and
query factors. The same scalar z is used in every stage map.

If F>=2, then

    dim ker(I + W_plus(-1)) <= Fq + R.

The bound does not require a formal query identity away from z=1. Thus it
also applies to an even faithful finite stage period, where z=-1 is an
actual character. For a constant monomial stage endpoint, the left side is
the number of odd stage powers, counted with multiplicity.

## Proof

The block normal form in PR6168 gives factors R_i=U_i V_i with

    V_i U_j = sign(i-j) I_q.

Factor each stage projector as P_i=A_i B_i, where B_i A_i is the identity
of its rank. This includes the extra stage-only maps. Let A concatenate
the A_i, and let B stack the B_i. Put G=BA, let L be its strict block-lower
part, and let D=I+2L. Chronological expansion at z=-1 gives

    W_plus(-1) = I - 2 A D^(-1) B.

Let k be dim ker(I+W_plus(-1)). The maps A and D^(-1)B give inverse
isomorphisms between the eigenvalue-one spaces of A D^(-1)B and
D^(-1)BA. Hence

    rank(D-G) = 2Fq + R - k.

Write Z=D-G. Its diagonal blocks are zero, and its off-diagonal blocks
are sign(i-j) B_i A_j. For a full call, compatibility gives P_i U_i=U_i
and V_i P_i=V_i. Compress Z at full calls using the block row factors
V_i A_i and block column factors B_i U_i. The (i,j) compressed block,
for distinct full calls, is

    sign(i-j) V_i P_i P_j U_j
      = sign(i-j) V_i U_j = I_q.

The diagonal compressed blocks are zero. Thus the compressed matrix is
(J_F-I_F) tensor I_q, with rank Fq for F>=2. In particular,

    Fq <= rank Z = 2Fq + R - k.

This is the stated inequality. The argument uses scalar coefficients and
matrix ranks. It makes no commutation assumption about the projectors.

The condition F>=2 is necessary. With two main fields, one full boundary
stage P=I and two thin queries, the stage endpoint is z I. Its two odd
powers exceed F=1. The usual three-query triangle supplies the query word.

## Consequence for the balanced alternating proposal

For n=2m fields per port, the alternating query word has F=m+2 full calls
and m-1 thin calls. The proposed profile with four degree-two stages and
n-4 degree-one stages has (2m-4)q odd powers. The bound therefore requires

    R >= max(0, (m-6)q).

For n>=14, that profile cannot have zero added stage rank. This holds even
with arbitrary coupling of all q ports. It does not exclude other power
profiles or identify an optimal invoice. At n=14 and q=104 the necessary
added rank is at least 104 stage coordinates. A physical rank-eight
advance on each such coordinate has its own cost. It cannot be omitted
from a kappa comparison.

The rank bound alone does not price the complete caller. It does not cover
extra query rank beyond the shortest-word contract, noncentral stage maps,
initialized work or an odd finite period without the -1 character. No
universal obstruction to kappa improvement is asserted.

## Open next questions

Can a different parity profile have a complete low-cost endpoint? First
check the native stage ranks and all extra maps. Can the rank argument be
extended to unequal query ranks or initialized work? The block normal form
used above does not establish either extension.

Credit the query/controller lineage of Cinna Davis and the retained native
and caller contributors. The block normal form is in PR6168. New stage-rank
argument: Ali Romman with substantial OpenAI Codex assistance. No subagents.
Apache-2.0. Independent scientific review is pending.
