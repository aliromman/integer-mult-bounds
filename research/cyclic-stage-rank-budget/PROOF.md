# A rank budget for a common cyclic stage

This is a scoped algebraic proof. It gives no caller or new kappa.
[FORWARD-TARIFF.md](FORWARD-TARIFF.md) proves the stronger exact forward
tariff over every formal cyclic order and every coefficient field. The
divided-difference proof here also applies to inverse-stage factors, since
their differences contribute the extra separable factor -1/(xy).

## 1. Affine factors

Let K be a field with p distinct p-th roots of unity. In characteristic zero
one can pass to a field extension to include the roots. Let

    W(z) = F_L(z) ... F_1(z),     F_i(z) = A_i + z B_i.

All matrices have size N. No commutation, idempotence or invertibility of
the coefficient matrices is needed for the next rank bound. Put b_i=rank B_i.
Choose disjoint lists lambda_1,...,lambda_a and mu_1,...,mu_b. Form the
block matrix

    L(W)_(j,k) = (W(lambda_j)-W(mu_k))/(lambda_j-mu_k).

Then rank L(W) <= sum_i b_i. To prove this, use the exact telescoping identity

    (W(x)-W(y))/(x-y)
      = sum_i F_L(x)...F_(i+1)(x) B_i F_(i-1)(y)...F_1(y).

Factor B_i=X_i Y_i with b_i inner coordinates. The i-th contribution to
L(W) is the product of a vertical block matrix and a horizontal block matrix
with inner dimension b_i. Rank subadditivity proves the claim. This also
covers degree-zero factors: they have b_i=0. Free constant input and output
changes of basis multiply L(W) by invertible block diagonal matrices.

This is the standard divided-difference rank argument. The contribution here
is its cyclic native-stage application and the resulting source-bound budget.
For the established relation between Loewner matrices and state realizations,
see B. D. O. Anderson and A. C. Antoulas,
[Rational interpolation and state-variable realizations](https://doi.org/10.1016/0024-3795(90)90140-8),
Linear Algebra and its Applications 137-138 (1990), 479-509. The direct proof
above states every fact used here; the reference is background and credit.

## 2. All cyclic powers at once

Partition all p roots into sets of sizes floor(p/2) and ceil(p/2). For an
integer 0<=d<p, the scalar divided-difference matrix of z^d has rank

    rho_p(d) = min(d,p-d).

For 1<=d<=floor(p/2), expand

    (x^d-y^d)/(x-y) = sum_(j=0)^(d-1) x^(d-1-j) y^j.

This factors through two rectangular Vandermonde matrices of full rank d.
Their inner square factor is invertible, so the product has rank d. For
d>p/2 put h=p-d. At the selected roots z^d=z^(-h), and

    (x^(-h)-y^(-h))/(x-y)
      = -x^(-h)y^(-h) (x^h-y^h)/(x-y).

Invertible row and column scalings reduce to the preceding rank h case.
For d=0 the matrix is zero. This proof works over every field in which the
p roots are distinct.

Suppose, at every p-th root z, the word has endpoint

    W(z) = C diag(z^d_1 I_(m_1),...,z^d_t I_(m_t)) D,

where C and D are fixed invertible matrices. Reorder the rows and columns of
L(W) by target coordinate. It is a direct sum of the scalar matrices just
computed, with the stated multiplicities. Therefore

    sum_j m_j min(d_j,p-d_j) <= sum_i rank B_i.       (1)

Only the residues d_j modulo p matter. The conclusion holds at every word
length and every positive multiplicity. It is not a finite search result.

## 3. The full-call application

A full call on q coupled ports has stage factor I-P+z P with rank P=2q.
This uses its reduced main interface. Every native helper remains autonomous;
its input is independent and its original output is retained. Accessing a
transformed helper as another stage-dependent channel changes this premise.
Its query, ordinary basis changes and any other allowed maps in this model
are constant in z. Multiplying the stage by these maps leaves the slope
rank at most 2q. Thus F full calls, with arbitrary constant maps between
them, satisfy

    sum_j m_j min(d_j,p-d_j) <= 2q F.                (2)

The constants can be noncommutative matrices. This includes query-coordinate
operators after a fixed finite-dimensional specialization, provided they
commute with the common stage and the stated slope-rank budget is retained.
No hypothesis fixes a shortest query frame or an alternating interior.

For the three-block advance of period six and the retained period-seven
native cell, the stage power d has Gaussian rank 8 min(d,p-d). The separate
query contributes rank one to each main target. The period-seven rank list is
8,16,24,24,16,8,0. The period-six list is 8,16,24,16,8,0.

The period-seven cell and its phase are specified in the predecessor
five-call-projector-rigidity/PERIOD.md, section 3. That proof establishes all
seven characters at every positive width. The three-block period-six stage
is U=R H_0; its powers act on one, two or three disjoint Hadamard blocks.
At every positive width, each block Hadamard has at least two independent
vectors of each sign. A length-three rotation orbit of positive vectors
gives all three cube-root characters of U. A rotation orbit of negative
vectors gives their negatives. Thus all six characters occur. The separate
query has both signs at widths at least two, which suffices for a contract
required at every positive width. These component facts do not replace a
full source-route replay.

Consequently cyclic wrapping cannot increase total stage rank beyond 16qF
in this model. Two proposed escapes are excluded for every q:

| Proposal | Required divided-difference rank | Allowed rank |
| --- | ---: | ---: |
| Six fields, eight full calls, powers 3,4,4,4,4,4 at period seven | 18q | 16q |
| Eight fields, ten full calls, five powers 3 and three powers 4 | 24q | 20q |

The twelve-field, fourteen-full-call profile with eight powers 2 and four
powers 3 has rank 28q on both sides of (2). It survives this necessary test.
The search archive has no complete endpoint for it. Equality in (2) is not
an existence proof.

Changing stage coordinates between calls, using maps which do not commute
with the stage, accessing transformed helpers, changing native rank profiles,
or changing the required
endpoint lies outside this application. A coefficient that is itself a
nonconstant function of z cannot be counted as a degree-zero ordinary map.

## 4. A favorable selected-source price ceiling

This part has an additional central-query, full-only hypothesis. All N main
inputs must receive the same independent query involution S, times a cyclic
stage power. At z=1 and S=-1, each full query changes identity by rank at
most q, and the target is -I_N. Thus N<=Fq by the rank bound for a product
minus identity. Ordinary constant basis changes do not alter this bound.
This step is not asserted for an arbitrary noncentral query boundary.

Put t=1-kappa with 0<=t<=1. Equation (2) and concavity give

    sum_main r^t <= N (1+16Fq/N)^t <= Fq 17^t.      (3)

For the last inequality, x(1+a/x)^t is increasing for x>0 when a>=0 and
0<=t<=1. Its derivative is (1+a/x)^(t-1)(1+(1-t)a/x)>0.

Use the selected PR6030 full-source invoice with q=104. Retain every helper
target of every call. Let P_t and H_t be its local paid and helper moments.
Even with all further costs omitted, the per-full margin is at most

    g(t)=104*17^t + H_t-P_t.

The authenticated profile gives

    g(t)= -2205 -449*2^t -492*3^t +7*4^t -149*5^t
          -9*6^t +3*7^t +18*8^t +385*9^t +104*17^t.

There is exactly one zero on 0<=t<=1. To check this, divide by 6^t and
differentiate. Every term except +7*(4/6)^t has nonnegative derivative.
The derivative of -2205*6^(-t) alone exceeds 2205/6, since log(6)>1.
The magnitude of the one adverse derivative is less than 7, since
log(3/2)<1. The whole derivative is positive. Also g(0)=-2787 and g(1)=48.

The exact rational replay locates the corresponding kappa zero between
0.00511060055 and 0.00511060056. All larger savings fail this favorable
moment. This is a ceiling, not an attained saving. The selected kappa stays
0.003125804, and its admitted BIT cap stays 0.0034. Ordinary maps, source
bindings, compiler completion, precision, stock and outer costs can only
reduce the favorable allowance used here.

## 5. A separate derivative screen

There is also a useful test that does not require the full rank calculation.
Suppose det W(z)=c z^D and W agrees with diag(z^d_j) at all p roots, where
0<d_j<p, counting every coordinate. If degree W<=L and L<p+min_j d_j, then

    D=sum_j d_j.

For L<p, polynomial interpolation already gives W=diag(z^d_j) when all
degrees fit below p. Otherwise write W=diag(z^d_j)+(z^p-1)K with
degree K<=L-p. At a p-th root lambda, Jacobi's determinant derivative gives

    D-sum_j d_j = p trace(diag(lambda^(-d_j)) K(lambda)).

Average over all roots. Every exponent on the right is strictly between
-p and zero, so its root sum is zero. This proves the statement. Constant
invertible boundary bases can first be removed. The Loewner bound already
rejects the two examples above; this derivative test explains the first
failed cyclic determinant relaxation.

## Credit and remaining work

The prior native stages, central involution, P104 source and current invoice
retain all predecessor credit, including Cinna Davis and the source authors.
The two-full border is prior work by Dillon/PR4623 and Yiyang/PR4966. The
present searches rediscover it; they do not supply a precision escape.

The new result is the all-length cyclic-stage rank budget, its exact source
ceiling, and its use to reject the two proposed targets. Ali Romman prepared
this work with substantial OpenAI Codex assistance. No subagents, independent
human review, proof-kernel acceptance or worldwide priority claim. Apache-2.0.
