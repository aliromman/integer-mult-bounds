# Exact forward tariff for every formal cyclic order

Let K be any field and let p>=2. In the full formal algebra K[z]/(z^p-1),
suppose a finite linear program must return independent outputs

    y_j = z^d_j x_j,       0<=d_j<p.

It can use arbitrary K-linear maps, any finite clean or dirty workspace,
and forward z calls of ranks r_i. All data paths, feedback, copies and
erasures are permitted in this algebraic relaxation. A rank-r forward call
has an affine matrix A+zB with rank B<=r; ordinary maps are constant in z.
No inverse-z call or nonconstant coefficient is free.

Then

    sum_i r_i >= sum_j d_j.                         (1)

Coordinate calls attain equality. Thus the exact minimum total query rank
is the sum of the least nonnegative residues. The theorem covers every
field, including characteristics which divide p. No root-of-unity or
semisimplicity assumption is needed. It concerns formal coefficient equality,
or a native representation faithful for that equality.

This extends the model boundary of Nick Petalas's four-inverse C3 result
and its restored-workspace question. It also contains the single positive
order-four inverse tariff as a special case. Those prior results and their
register-model lineage keep their credit. The allowance of every nonidentity
cyclic power at unit cost in Wentao Li's cyclic 2-power theorem is a different
alphabet; (1) does not assign unit cost to those powers.

## 1. A finite state representation of the literal word

First use a fixed finite ambient space for all live and private fields.
Initial embeddings and final projections are constant maps. Thus clean
initialization and observed-output selection do not add query rank. Ignore
zero-power targets and restrict to the N positive-power inputs and outputs.
Granting zero values to every other original input only weakens a lower bound.

For the chronological factors F_i(z)=A_i+zX_iY_i, use r_i inner coordinates.
There is a strictly lower block triangular matrix T of size
r=sum_i r_i, constant matrices B and C, and a constant D such that

    W(z)=D+z C(I-zT)^(-1)B.                         (2)

The blocks are explicit:

    T_(j,i) = Y_j A_(j-1)...A_(i+1) X_i   for j>i,
    B_i     = Y_i A_(i-1)...A_1,
    C_j     = A_L...A_(j+1) X_j,
    D       = A_L...A_1.

Put zero in the other blocks of T. Constant input/output maps can be
absorbed in B,C,D. Expanding (I-zT)^(-1) as a finite geometric series lists
exactly every ordered nonempty selection of zX_iY_i terms in the product.
This proves (2). T is nilpotent, including when intermediate constant maps
are singular. Larger rank factorizations give a weaker bound but remain valid.

## 2. Remove cyclic coefficient aliasing

Set

    B_star=(I-T^p)^(-1)B,
    V(z)=z C(I-zT)^(-1)B_star.

Both inverse expressions are finite sums because T is nilpotent. If W_k
is the coefficient of z^k in (2), the coefficient of residue k, for
1<=k<p, in its formal reduction modulo z^p-1 is

    sum_(a>=0) W_(k+ap)
      = C T^(k-1)(I-T^p)^(-1)B.

It is precisely the coefficient of z^k in V. The ordinary feedthrough D
contributes only to residue zero. The assumed endpoint therefore gives

    V(z)=diag(z^d_1,...,z^d_N)+O(z^p),              (3)

where every d_j is now in 1,...,p-1. This is an exact formal-power-series
statement, not a limit or an approximate endpoint.

## 3. Count the first determinant term

Divide column j of (3) by z^d_j. The result is I plus a matrix whose entries
have positive valuation: p-d_j>=1. Consequently det V is nonzero and its
first nonzero term is

    z^(d_1+...+d_N).

On the other hand, the block determinant identity gives

    det V(z) = det [[I-zT, z B_star], [-C, 0]].      (4)

Here det(I-zT)=1. Each of the first r rows in the block matrix has degree
at most one, and the final N rows are constant. Thus the polynomial in (4)
has degree at most r. Its first nonzero degree cannot exceed its degree.
This proves (1), in every characteristic.

No divisibility, algebraic closure, invertible scalar gate, bounded
coefficient or absence of workspace was used. The endpoint on independent
original inputs and the forward affine alphabet are essential hypotheses.

## 4. Consequences for the current research

- Four inverse targets for a positive order-three generator still need
  eight rank-one calls with any finite clean or dirty workspace and any
  coefficient field. More generally, N inverse targets need (p-1)N.
- Mixed cyclic powers have the exact additive tariff sum d_j. Determinant
  congruence modulo p is not an adequate budget screen.
- A common forward full stage on q coupled ports has slope rank at most 2q.
  Thus a target with stage residues d_j needs sum d_j<=2qF. This strengthens
  the absolute cyclic rank budget when a target uses a residue above p/2.
- The twelve-field, fourteen-call profile with eight powers 2 and four
  powers 3 is still tight and still open. The theorem does not construct it
  or settle its query endpoint.

An inverse stage or another paid generator changes the alphabet. A constant
in the proof must commute with the chosen native operator so that scalar
substitution represents the physical word. Native address maps which fail
that condition, nonlinear instructions, nonfaithful quotients and restricted
input promises are outside the theorem. A concave price for one higher-rank
invocation also needs its own conversion from total query rank to money or
the repository's recursive moment.

The state representation and block determinant identity are standard linear
algebra. The new contribution is this all-order, all-field tariff proof and
its application to the cyclic caller search. No worldwide priority claim,
new kappa or complete supplier admission is made. Ali Romman, with substantial
OpenAI Codex assistance. Apache-2.0.
