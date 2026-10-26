# Paid helper subdivisions

## Fixed source and claim

The source is PR5718 at a1534ac45401c9a9c80c07aa65c89e00e33594f8.
Its canonical word SHA256 is
`a1c8277189a811546bdc7b147d72e6a594504ff9585c9bf59a5681aa873abb01`.
The final word and all intermediate words have fixed canonical hashes in
PLAN.json. A canonical hash uses sorted JSON keys and compact separators.

All matrices below act on decoded scalar values. A MOVE changes the native
frame, with its complete rank and absolute phase checked separately. The source
public endpoint matrix is T = I + sum e_(x+112) e_x^T over the 104 live ports.
The public helper r has an identity row and column in T. Owned input columns
are erased. All fresh helper inputs are arbitrary and independent.

## One paid split

Choose an original public MOVE from frame A to frame B. Choose a cut c inside
a chronological window with no scalar read or write of r and no second MOVE of
r between the original MOVE and c. Choose an existing frame F with strict
containments A < F < B. This permits the original MOVE to occur at c.

Let R be the scalar prefix at c and S the suffix. The suffix column d = S e_r
records the effect of an injection at the cut. Set a = e_r^T R and
D = T^(-1) d. Choose two corrections v and u that read other fields in the
actual frame F at c. Define

```text
h = a + v R,     j = a + u R,     k = h - e_r^T,
g = (j - e_r^T) T^(-1).
```

The exact constructor checks:

- h_r = j_r = d_r = D_r = 1 and hD = jD = 1.
- All entries of h and D are public and have the same initial frame as r.
- All entries of d and g are public and have the same final frame as r.
- g_r = 0. Neither v nor u reads r.
- Every nonzero term in v and u uses its actual current value at c.

For the seventeen fixed witnesses, h = e_r^T and u = 0. The general identities
are retained in the constructor and checked exactly; no search occurs in replay.

Apply the initial basis change

```text
Q = (I + (D-e_r)e_r^T)(I-e_r k)
  = I + (D-e_r)e_r^T - D k.
```

First emit r -= sum k_s x_s. Then emit x_s += D_s r for s != r.
All operands share the declared initial frame. Since k_r = kD = 0,
these are legal invertible ordinary shears, Qe_r = D, hQ = e_r^T and jQ = j-k.

At c, replace the original MOVE by these operations:

1. Pay for MOVE(r, A, F).
2. Add v times the actual coframed fields to r. Its value becomes x_r; retire r.
3. Start a new independent public helper with value z in F.
4. Subtract u times the actual coframed fields from the new helper.
5. Pay for MOVE(new, F, B).
6. Use new in every subsequent source occurrence of r.

The original MOVE is removed. Owned identifiers increase by one. Their complete
INIT/COPY/use/CLEAR lifetimes remain unchanged. This operation creates no owned
zero and makes no assumption that a historical scalar value is still available.

## Endpoint repair

Replacing the continuation value adds d times (z-jQx) to the source suffix.
Use TQ = T + (d-e_r)e_r^T - d k. Every other public endpoint is then

```text
(Tx)_s + d_s tau,       tau = z - (j-e_r^T)x.
```

The new helper has value tau; the retired helper has value x_r. At the original
final frame, subtract d_s times new from every other affected output s. This
restores its exact endpoint (Tx)_s. Then add g times the restored outputs to
new. Since g_r = 0 and gT = j-e_r^T, the new value is z. The order is necessary:
the g correction reads the repaired endpoints. All dirty owned input columns
still map to zero. The verifier checks every coefficient after each step.

## Seventeen fixed changes

The helper roles, in order, are
731, 723, 722, 725, 724, 729, 727, 728, 730, 715, 712, 714, 713, 717, 718, 719, 720.
PLAN.json specifies each source hash, event index, cut, intermediate frame and
correction. WITNESS.json records h, j, d, D, g and all new helper identifiers.
The verifier computes these quantities over the rational numbers. It also checks
the final word with a separate reverse-chronology dense integer implementation.

Write H_r = r^t, where 0 < t < 1. One paid H2 becomes H1+H1. Sixteen paid H6
calls become H5+H1. Each change replaces a productive H9 by H8+H1. Thus

```text
Delta paid    = 18H1 - H2 + 16H5 - 16H6,
Delta targets = 17H1 + 17H8 - 17H9,
Delta(targets-paid)
  = (H2-H1) - (H9-H8) + 16[(H6-H5) - (H9-H8)] > 0.
```

The increments of x^t strictly decrease with x, which proves the last inequality
throughout 0 < t < 1. Rank mass is unchanged. The eighteen local caller
occurrences multiply both deltas by eighteen. All 15 query shears and 3,120
whole idle jobs remain in the complete bill.

The final local paid profile is
`{1:2202, 2:458, 3:495, 4:92, 5:175, 6:6, 7:15}`.
The productive helper profile is
`{1:19, 2:5, 4:91, 5:6, 7:6, 8:32, 9:394}`.
The complete bill has 66,654 paid calls, 13,282 targets and rank surplus 760.
The maximum paid rank is 72; the maximum target rank is 145.

## Complete price and finite bounds

Two exact interval engines verify the new moments. The raw crossing lies in
(0.00311224343152, 0.00311224343153); a negative derivative bound proves uniqueness
on [0, 0.0034]. The declared kappa is 0.00311224341 and the tested saving is
0.003112243415. Both lie strictly below the lower endpoint of the crossing.

The old completion parameters fail after the target count increases. The new
choice is cutoff c=16, floor b=17, tail offset u=7/50, power v=1024, reserve
10^(-12), decay 122/125 and tail block 128. At the tested saving, certified
upper bounds for the three charged ratios are

```text
alpha <= 0.999999999983223032330514874511 < 1 - 10^(-12),
gamma <= 0.975623625359102786173183551859 < 122/125,
delta <= 0.936724728289985271710889269085 < 122/125.
```

The verifier also checks decay^128 < 1/16, c^64 < b^63, exit and complement
capacity, the finite uniform threshold, rounding charges and all 54 strict plus
four functional outer inequalities. The normalized paid upper bound is
1048576000976562501/1048576001073741824 < 1. The inherited ordinary BIT admission
0.0034 remains above the tested saving.

The local envelope is explicit: at most 768 public fields, nine owned labels,
peak two live owned fields and fewer than 28,000 instructions. The final word
uses 761 public fields and 27,502 instructions. LOCAL-EDITS.json changes only
the old public and operation caps; all semantic and cost checks remain.
The 770 total field labels fit the retained 1,000-slot allowance. The complete
local ordinary price is 244,338 scalar units plus 5,217,000 route units, below
10^7. The complete nonrecursive bound is 184,647,552, below 10^9. The scalar
growth bound is updated to 2^(16*28000). The native paid rank mass remains 5,987.
The verifier checks these charges before it accepts the exponent.

## Verification boundary

The finite checks cover 592,900 endpoint coefficients in both directions, all
2,912 distinct native maps and absolute phases, the complete physical caller,
all paid and target occurrences, and the full conditional price. Fault controls
remove actual paid or repair operations, alter a real coefficient, omit a CLEAR,
misstate boundary frames and falsify invoice entries. They must fail under normal
and optimized Python.

The five inherited framework assumptions remain. No large uniform word or full
rounded multiplier is materialized. Source review receipts do not establish
independent review of this new word. A bounded discovery search over existing
intermediate frames and windows without pivot scalar use found these seventeen
changes; it is not an optimality proof for other frames or words.

Prepared by Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.
