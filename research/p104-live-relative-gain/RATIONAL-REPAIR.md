# One missing direction at chart 117

Five paid transfers in the rational circuit can be replaced by one transfer from chart 127 to chart 117. This leaves 5,696 paid transfers and changes the return for helper 1093 to a local identity at chart 117. The modified graph delivers 1,288 independent directions. To use its quoted tariffs in the 48-call circuit, a repair must deliver the remaining direction.

Let M be the required invertible map on the 1,289 scalar fields. In the modified graph the exact map is

```text
Lx = Mx − (Mh)(qᵀx),    qᵀh = 1.
```

Thus Lh = 0, and Lx = Mx for every input with qᵀx = 0. All 122 nonzero coordinates of h start at chart 3. All 305 nonzero coordinates of q are helper fields from chart 3 to chart 127. For example, an encoder at chart 3 can solve qᵀx = 0 for one of these fields. This gives a valid restriction with one fewer independent grade-nine return. At the trial that return is worth more than 8.92 per call.

The original circuit owes 36.15087585448 across 48 calls at κ = 0.00391884648. The modified graph has more than 2.2479132437 and less than 2.2479132439 per call available to restore the missing direction if it preserves all the listed endpoints. This calculation omits income for the local return at chart 117. It includes all other endpoint income on the condition that a repair makes those endpoints independent. A direct grade-nine delivery exceeds the available price. A different encoder can change the endpoint prices and must account for those changes.

## Reproduce the equation

```sh
python3.12 -I -B research/p104-live-relative-gain/verify_projection.py
```

The delta removes paid IDs 91075, 91373, 91580, 91582 and 103243. It adds ID 95 and specifies each changed scalar coefficient. The checker rebuilds the graph from the verified parent and checks its SHA-256. It uses integer arithmetic to check the map on every row and column, including the donor chart of each operand. The forward replay visits 13,380,291 terms. It then prices every transfer and listed endpoint with outward rounded arithmetic. [Compact delta and both vectors](Q117-DELTA.json.gz) · [Checker](verify_projection.py)

## A repair with one grade-five transfer

Chart 3 forms the scalar qᵀx from its public inputs. Send it to chart 117 on a fresh transfer, ID 96, before the other transfers. Add a multiple of this scalar to six existing rows at chart 3 and nine at chart 117. The delta gives the 15 coefficients. In the repaired circuit, let tₖ denote the coefficient of the extra signal at row k. Each tₖ is the sum of its local correction and the t values of its donors times their original coefficients. The checker evaluates this recurrence in dependency order and compares all 1,289 terminal values with Mh over the rationals. Thus the 5,697-transfer circuit computes M on independent inputs at the stated endpoints.

The new transfer has grade five. At the trial its price exceeds the available margin by between 2.7206502969 and 2.7206502971 per call. The calculation includes all 48 occurrences. This gives a constructive price for a full repair. A cheaper supply of the same signal or a separate transfer saving is still needed.

## Open next question

Can existing receipts deliver the correction (Mh)(qᵀx) within the available price? A repair needs to supply the correct coefficients to each affected output in a valid order. The input coform is available at chart 3. A low price for moving it once does not establish that the recipients can use it at their required charts.
