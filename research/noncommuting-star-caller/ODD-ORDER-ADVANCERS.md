# Exact advances of orders seven and nine

These are component constructions for the current caller search. They do not
supply a profitable caller or a new kappa. They change the stage relation from
the earlier whole-block permutation advances of orders four and six.

## Three-coordinate cells

Use binary coordinates x=(x_0,x_1,x_2). Let H_0 be the normalized real
Hadamard on coordinate zero. Define the binary address map

    A x = (x_0, x_0+x_1, x_1+x_2).

The additions in this address formula are modulo two. Define ordinary
monomial maps by

    R_7 |x> = i^(x_0+x_2) |A x>,
    R_9 |x> = i^(x_1+x_2) |A x>.

Put V_7=R_7 H_0 and V_9=R_9 H_0. Literal Gaussian-integer multiplication,
with a denominator sqrt(2) for each factor, gives

    V_7^7 = (1-i)/sqrt(2) I,
    V_9^9 = -(1+i)/sqrt(2) I.

The complete Gaussian-rank sequences of their positive powers are

    V_7: 1,2,3,3,2,1,0,
    V_9: 1,2,3,3,3,3,2,1,0.

Two checks give these ranks. One multiplies binary symplectic matrices and
computes the rank of their upper-right block. The other multiplies all exact
8-by-8 complex matrices and counts the nonzero entries in each row. Every
row has 2^r entries, with the same r as the symplectic check. The complete
power at the stated order is a scalar, not a residual Pauli map.

## A rank-eight stage at every positive width

Place eight disjoint copies of a cell on24 axes. Set

    U_7 = V_7^(tensor 8), U_9 = V_9^(tensor 8).

The eighth power of either scalar phase above is one. Consequently
U_7^7=I and U_9^9=I, with no residual phase. Both operators have Gaussian
rank eight. Their second powers have rank16 and their third powers rank24.
The eighth tensor power also removes the sqrt(2) denominator issue: the
entries are Gaussian dyadic numbers. The identities persist under another
tensor power of any positive digit-column width f.

In execution order, either U is the existing eight-axis H_8 followed by
eight ordinary monomial R cells. Thus the source wrapper in
ADVANCING-CONTRACT.md can use these maps in place of its final whole-block
permutation. The new monomial acts on24 transverse axes. It does not act
on the central coordinate. It therefore commutes with the central
involution S defined in that note.

The inverse stage has the same native Gaussian price. If U=R H_8, then

    U^-1 = R^-1 U R^-1,

because H_8^2=I. The two R^-1 wrappers are ordinary. This statement supplies
an advancing-stage inverse, not an inverse of a source COPY or CLEAR.

For e=2,3, the target U^e S is an ordinary equivalent of a Gaussian job of
rank17 or25. Here S=QG, and G is the ordinary whole-column map proved in
ADVANCING-CONTRACT.md. The standard native Gaussian normal form uses the
checked symplectic rank. The actual104-port routes, their phase anchors and
their ordinary allowances still need a source-bound caller replay.

## Search use and boundary

`clifford_advancers.py` enumerates all10,752 three-coordinate cells
R H_0 with an invertible binary address map and a symmetric binary phase
matrix. It retains22 power-rank profiles and one literal example for each.
Normal and optimized runs give identical receipts. The displayed two
examples have the exact phases above. No claim about all larger Clifford
advances follows from this finite census.

A constant-coefficient caller with an order-p stage can be tested in the
algebra with D^p=1 and S^2=1. For example, eight fields and nine full calls
could request six rank17 and two rank25 main jobs. With the pinned source
inventory, this unconstructed profile has linear target-minus-paid margin
328. It is only a price screen. The exact endpoint has not been found.
The admitted BIT cap remains .0034 even if a complex caller passes.

CYCLIC-REDUCTION.md removes redundant short searches. For seven calls at
period seven and target degrees two and three, cyclic reduction cannot
change the formal-stage identity. Nine calls at period seven lie outside
that lemma for the balanced target. A separate signed test uses four inverse
stages; its normalized target has degrees zero and one. Passing any finite
field screen still requires a characteristic-zero word and full pricing.

## Open next questions

- Find an exact caller in one of the domains that the cyclic lemma leaves open.
- Bind a surviving word to all104 actual source ports, with fresh helpers,
  target normalizers, COPY/CLEAR, routes and precision costs.
- Test ordinary frames that depend on R. They require a noncommutative operator
  algebra; the constant-frame cyclic polynomial checks do not cover them.

Credit the P104/native source authors, Cinna Davis's full-query model and
controller lineage, and the source pins in this package. Binary symplectic
normal forms and Clifford rank calculations are standard tools. The cell
census and this application were prepared by Ali Romman with substantial
OpenAI Codex assistance. No global priority, independent human acceptance,
complete native caller or selected-frontier change is claimed. Apache-2.0.
