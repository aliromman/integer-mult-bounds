# A stage-degree loss with a minimum number of queries

This result allows any finite number of initialized work fields and a
constant output projection. It applies to a formal forward stage and a
central involution. It does not exclude signed stages, noncentral frames,
or identities that use a finite-order stage relation.

Let the scalar field have characteristic other than two. There are n input
fields, any number of work fields, and exactly n rank-one query calls. All
other circuit matrices are independent of the query symbol. The work starts
at zero. A constant input embedding E and output matrix L specify the n
retained outputs. The work outputs can be discarded.
Constant linear maps can also occur between calls. They can include copies,
projections and clears. Their exterior maps have degree zero, so they do not
change the argument. The proof grants these maps without a cost; an actual
implementation must still obey the ownership and cost rules.

## A query-free transfer test

Replace the involution S by a formal scalar t. Every such linear query
circuit has the transfer form

    M(t) = D + t U (I-tK)^(-1) V,

where K is strictly lower triangular of size n. To see this, let h_i be the
scalar input to query i before that query acts. Then

    h = V x + t K h,       y = D x + t U h.

Strict chronology makes K strictly lower triangular. All matrices may
have entries in the rational function field of the stage variable z.
D is the complete transfer with every query disabled. The inverse is a
finite polynomial because K is nilpotent.

Suppose the required endpoint is M(1)=T and M(-1)=-T, with T invertible.
Put H=(I-K^2)^(-1). Addition and subtraction give

    T=U H V,       D=-U K H V.

There are n queries and n outputs, so U and V are square and invertible.
It follows that

    D T^(-1) = -U K U^(-1).

Thus D T^(-1) is nilpotent. In particular det D=0. More strongly,

    M(t) = U (tI-K)(I-tK)^(-1) H V,
    det M(t) = t^n det T.

Here det(tI-K)=t^n, det(I-tK)=1 and det H=1 because K is strictly
lower triangular. This is an identity in t, not only a comparison at the
two query characters. This argument permits
arbitrary initialized work, arbitrary constant output mixing, and all
query-to-query feedback allowed by chronology. It uses the exact two
query characters, not a linearized query test.

## Degree available from the full calls

Now each call is the native factor

    F_i(z,t) = I-P_i + z P_i(I+tN_i),

where P_i has rank two, N_i has rank one and is square-zero, and
P_i N_i=N_i P_i=N_i. The ambient field count can exceed n. Let

    M(z,t)=L F_n(z,t) ... F_1(z,t) E.

The induced map on the n-th exterior power of each F_i has total degree
at most two in the variables z and t. To prove this, choose a basis in
which P_i is diag(I_2,0) and N_i is a square-zero map on the first plane.
The induced map is the direct sum of terms that use zero, one or two
vectors from this plane. The zero-vector term is independent of z and t.
The one-vector term is z(I+tN_i), of total degree at most two. The
two-vector term is z^2 det(I+tN_i)=z^2. A constant change of basis
preserves this degree bound.

The Cauchy-Binet formula now shows that det M(z,t) has total degree at
most2n. The rectangular E and L cause no problem: their exterior maps
are constant. The exact query-free transfer identity gives

    det M(z,t) = t^n det T(z).

Consequently

    degree(det T(z)) <= n.

In particular n full calls cannot produce n independent S z^2 targets
by any such initialized-work caller. More generally, targets S z^(e_j)
with constant invertible input and output bases must satisfy

    sum_j e_j <= n.

The proof is independent of the number of work fields. It is stronger
than the square-word ordinary-residue bound in ORDINARY-RESIDUAL-BOUND.md.
The bound is sharp. Use n independent pairs (x_i,y_i), with y_i initially
zero. One full call on each pair gives (z x_i,z S x_i). Retain the second
output from each pair, then install these outputs on the original fields
by ordinary swaps and clear the owned fields. This gives T=z I_n and
sum_j e_j=n. It uses n work fields. This exact circuit is a bound control,
not a profitable native caller.

## Mixed calls and coupled ports

The same proof applies to a query of rank q and a full projector of rank2q.
Split its query response into q scalar components. The chronological K is
strictly block lower triangular. Its total size is the sum of the query
ranks. If this sum equals the number n of independent required outputs,
the square U,V determinant proof is unchanged.

On the full projector image, a basis puts the square-zero query into q
two-dimensional Jordan blocks. Each block contributes total exterior
degree at most two, by the preceding argument. Thus a full call contributes
at most2q to the total degree. A thin rank-q query contributes at most q.
If A is the sum of the query ranks in the full calls and B the corresponding
sum in the thin calls, then n=A+B and

    degree(det T) <= (2A+B)-n = A.

This permits arbitrary constant mixing between the ports and any finite
initialized work. It does not require n to be a multiple of104. It does
require the same central query and formal common advancing variable. No
other nonordinary stage map is granted.

## Complete P104 tariff consequence

Use the exact selected source word in SOURCE-PINS.json. One full call has
query rank104 and the paid and independent-helper profiles

    paid:    {1:2203,2:459,3:493,4:92,5:168,6:18,7:14},
    helpers: {1:2,2:5,4:89,5:11,7:5,8:16,9:412}.

Write c_r=paid_r-helpers_r. The source has

    sum r*c_r=1720,
    sum c_r*max(r-j,0)<=0 for each integer j>=1.

The verifier reconstructs the profiles from every MOVE and each actual
independent helper boundary in the authenticated source word. It checks
every nonzero hinge sum. The concave-tariff argument is retained from
shortest-controller-stage-capacity/PROOF.md: for any nondecreasing concave
g with g(0)=0 and g(1)=1, these hinge inequalities give sum c_r*g(r)>=1720.
Also g(r)<=r at each nonnegative integer r.

Take F full source calls and B rank-one thin queries, and require exactly
n=104F+B independent central-query outputs. Each full call has fresh
independent helpers and its complete source bill. No helper is mixed into
the main query code. There are no separately supplied stage or idle maps.
The degree bound gives sum e_j<=104F. Since one stage has Gaussian rank8,
each target S z^(e_j) has rank at most1+8e_j. Its actual rank can be lower.
The complete main-plus-helper target moment minus paid moment is at most

    n + 8*104F - 1720F - B = -784F.

Thus F>0 fails every stated concave tariff, including g(r)=r^(1-kappa).
At F=0 the bound is zero and gives no strict recursive saving. Additional
nonnegative ordinary, ownership, route or compiler costs cannot help.
This closes the minimum-total-query-rank formal caller route on this
fixed source, even with initialized work and cross-port scalar mixing.
It does not close a caller with extra query rank, signed or cyclic stages,
noncentral ordinary operators, changed helper use or a different source.

## A quantitative bound with extra query rank

The same transfer gives a bound beyond the minimum. Let the total query
rank be Q, the number of independent outputs be n=Q-h, and the number
of chronological query layers be m. Each layer can have an arbitrary
query rank. Then K is strictly block lower triangular and K^m=0.
As before, T=U H V is invertible and D=-U K H V, with H=(I-K^2)^(-1).
Put W=H V T^(-1), so UW=I_n. Select a complementary pair A,B with

    [W A]^(-1) = [U;B],       BA=I_h.

For R(t)=(tI-K)(I-tK)^(-1), the transfer is M(t)T^(-1)=U R(t) W.
The complementary-minor identity gives

    det(M(t)T^(-1)) = t^Q det(B R(t)^(-1) A).

The inverse R(t)^(-1)=(I-tK)(tI-K)^(-1) has a pole of order at most m
at zero, because K^m=0. Its h-by-h compressed determinant therefore has
pole order at most mh. Consequently

    t^max(0,Q-mh) divides det M(z,t).

For full-query rank sum A_full and thin-query rank sum B_thin, the
exterior total-degree bound is2A_full+B_thin=Q+A_full. Evaluation at t=1
now yields the necessary stage-degree bound

    degree(det T) <= A_full + min(Q,mh).

For h=0 this is the sharp minimum-query theorem above. The argument
allows rectangular inputs, discarded initialized work and arbitrary
constant maps at every gap. It does not assert that the bound is attained
by a native caller when h>0.

For F full P104 calls and no thin queries, Q=104F and m=F. If
0<=h<=104, the complete concave-margin upper bound becomes

    -784F + h(8F-1).

In particular every 0<=h<=98 gives a strictly negative bound, at most-98.
Thus a formal all-forward caller with n=104F-h independent central-query
outputs cannot improve the source at any stated concave tariff when it
loses at most98 output fields from its total query rank. This is an
unbounded exclusion in F, not a finite search. The first possibly positive
linear cases are h=99,F>=13; h=100,F>=7; h=101,F>=5; h=102,F>=4;
and h=103 or104,F>=3. These are necessary screens only.

There is also a bound at the current exponent scale. For an upper bound,
add h hypothetical degree-zero targets to reach104F targets. Keep the
stage-degree budget F(104+h). Discrete concavity then gives

    actual main moment <= F*((104-h)g(9)+h*g(17))-h,
    complete margin <= F*((104-h)g(9)+h*g(17)-sum c_r*g(r))-h.

The hypothetical fields are only a proof device; no caller receives their
productive credit. The verifier certifies the following power tariffs
with rational intervals whose endpoints pass exact integer-power tests:

- For h=101 the expression in parentheses is negative at kappa=.003.
- For h=102 it is negative at kappa=.0034.

Each expression has negative coefficients only at ranks1 through7 and
positive coefficients only at ranks8,9,17. Division by8^t makes it
nondecreasing in t, so the exclusions persist for every larger kappa<1.
Reducing h also decreases the expression. Therefore a caller that beats
the selected .00310045 needs h>=102 in this model, and one that reaches
.0034 needs h>=103. These bounds hold for every positive full-call count,
not only the searched sizes. They do not construct a caller at either
threshold.

The same total-degree argument permits different independent formal
commuting stage variables: set each one equal to z. A repeated finite-order
stage relation, a signed stage, and a noncentral coefficient remain outside
the model. Fresh helper independence and the fixed source invoice remain
required for the P104 price.

## Implication for the retained search

The four-input, two-work-field, four-call target with all degrees two is
impossible in this formal central model. Its floating search had residual
norm0.025872003514433046 and maximum parameter113.34257277821634 after the
retained trials. Those numbers alone were not an exclusion. This proof is
the reason to close that exact target.

The briefly proposed degrees(2,1) target with two calls is excluded too.
The final sharp bound supersedes the preliminary degree2n-1 bound used
to choose that test. A next profitable caller must use an additional
query, signed stages, a noncentral frame or an actual cyclic relation.
Each option has a
different full price. No target receives credit before its exact endpoint
and native source binding are proved.

## Checks and open next questions

- Check the nilpotent transfer identity and exterior coefficient identity
  in independent exact arithmetic. The retained verifier checks scalar and
  rank-two query blocks and the actual source tariff.
- Determine whether the extra-query degree bound can be attained by a
  profitable native caller. The two-field three-call triangle is a
  positive control with one extra query.
- Test whether a signed or cyclic initialized-work caller has a profitable
  exact endpoint and complete cost accounting.

Ali Romman with substantial OpenAI Codex assistance; no subagents. The
native full-call model and rank-one chronological method retain the
predecessor attribution in LAB.md and SOURCE-PINS.json. The returned dirty
scratch question in Discussion4902 is related context with a different
boundary contract. No worldwide priority, independent acceptance or new
kappa is claimed.
