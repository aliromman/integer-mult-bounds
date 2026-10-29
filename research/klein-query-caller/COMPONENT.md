# Three equal-rank commuting query involutions

This is a native component construction and a changed caller hypothesis.
It is not a complete P104 route, caller, or kappa result.

For three binary m-bit blocks a,b,c in Fourier coordinates, define

    q1(a,b,c)=a dot b,
    q2(a,b,c)=a dot c,
    q3(a,b,c)=a dot (b+c).

Let S_i be Fourier conjugation of the diagonal phase (-1)^q_i. Then

    S_i^2=I,       S1 S2=S3,       S1 S2 S3=I.

The phases commute and q1+q2+q3=0 over F2. Each polar matrix has rank 2m.
Thus each operator has the same Gaussian rank 2m. Swapping b and c, or
replacing b by b+c, relates the three operators by ordinary binary address
permutations. This gives three native query types with one common rank.

## Faithful characters

Put D=2^(3m) and T=2^(2m). The trace of each S_i is T. To see this, sum
(-1)^(a dot b) over all a,b,c. The sum over b is zero unless a=0, in which
case all 2^(2m) choices of b,c remain. The other two phases have the same sum.
The four simultaneous character multiplicities are therefore

    (D+3T)/4 for (S1,S2)=(1,1),
    (D-T)/4  for each of the other three characters.

They are positive for every m>=1. A complete endpoint test must include
all four characters. No character makes all three S_i negative.

## Relation to the retained P104 query

The preceding noncentral-query-stage component starts with

    C=((1+i)I+(1-i)X)/2,    Q=C^(tensor f),
    G=(I+X_all)/2-i(I-X_all)/2.

For even f, S=QG has Fourier phase

    q_f(y)=sum_(i<j) y_i y_j.

G is the retained two-address ordinary normalization. Its native proof and
source credit remain in that predecessor. This section gives an explicit
ordinary change from q_(2m) to the hyperbolic phase a dot b.

Index the input pairs by k=0,...,m-1. Put

    t_k=sum_(j>2k+1) y_j,
    a_k=y_(2k)+t_k,       b_k=y_(2k+1)+t_k,
    ell(a,b)=sum_k (k modulo 2)(a_k+b_k).

This triangular binary change is invertible. Directly separating the first
pair gives

    q_(2m)(y)=a_0 b_0+q_(2m-2)(y_tail)+sum y_tail.

Induction yields

    q_(2m)(y)=sum_k a_k b_k+ell(a,b).

The linear Fourier phase ell is an ordinary physical translation. The
binary change of coordinates is an ordinary physical address permutation
with the inverse-transpose rule. After those two maps, the normalized native
query has phase a dot b. Adding the independent c block gives S1; the
two address changes above give S2 and S3. No unsigned phase convention or
small-width extrapolation is used in these identities.

This uses source width f=2m and an additional m-bit spectator. Three
disjoint advancing blocks of 8f bits, together with the 3m query bits, fit
in a 51m-bit field. The source word has to be lifted to the spectator and
all routes have to use this changed coordinate assignment. Those physical
and compiler obligations have not been replayed. Equal component Gaussian
ranks do not make the address adapters or the complete source call free.

## Caller question and favorable price

Take three independent groups of six main fields. Their required query
targets are S1, S2 and S3 respectively. Each group has the retained exact
seven-transvection canonical query word. Interleave the three groups in
round order. This gives eighteen independent query targets in twenty-one
calls, before attaching stages.

The complete stage target would have twelve powers two and six powers three
of the common period-six stage. Every full call has a compatible rank-two
stage projector. The total target stage degree is 42, equal to twice the
full-call count. Query and stage ranks are both multiplied by the source
width scale. The favorable selected-source local moment is the same per
call as the twelve-field/fourteen-call profile. Its raw crossing is about
.00343212, before ordinary, native, precision, compiler and outer charges.
The admitted BIT cap is .0034. Neither number is an achieved saving here.

The source can retain every helper as an independent autonomous target in
the proposed algebraic boundary. A complete construction must prove that
boundary with the actual native coordinates and include every helper,
adapter, ordinary operation, owned lifetime and restoration charge.

## Credit and limits

The Q/G normalization, P104 query and source lineage retain their credit,
including Cinna Davis and the preceding native/caller packages. Hyperbolic
quadratic forms, Fourier conjugation and character decomposition are standard
tools. The changed three-query caller contract and tests here are by Ali
Romman with substantial OpenAI Codex assistance. No subagents, independent
scientific review or new kappa. Apache-2.0.
