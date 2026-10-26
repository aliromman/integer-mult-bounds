# Exact construction and complete price

## 1. One temporary basis

Let `a` and `b` be distinct public registers. Set `R = I + c E_ab` for a
nonzero dyadic scalar `c`. Then `R^-1 = I - c E_ab`.

Choose two internal cuts where both registers can reach one common chart.
Insert `R` at the first cut. Between the cuts, replace each scalar block
`M_j` by `R M_j R^-1`. Insert `R^-1` at the second cut. The complete interval
is unchanged:

```text
R^-1 (R M_k R^-1) ... (R M_1 R^-1) R = M_k ... M_1.
```

A block that does not use either register commutes with `R`. The same is
true for each owned initialization and erasure, since `a` and `b` are public.
The checker requires one certificate for every affected block. It checks
the exact old event list and every rational matrix coefficient. It also
rejects extra certificates. Thus this is an identity for arbitrary inputs,
not a test on selected vectors.

Each scalar block has one physical chart. The compiler moves each operand
to that chart before it is used. At every step, its physical value is the
declared native chart applied to its logical value. A scalar basis change
at a common chart therefore implements the same rational matrix `R`.
Native transport changes the chart and preserves the logical value.
The final basis removal restores the original logical interval.

## 2. Fixed witnesses

The four certificates are applied in the order in WITNESS.json. Block
indices refer to the complete word produced by the preceding step.

| Public pair | Coefficient | Local paid delta |
| --- | ---: | --- |
| 438, 733 | -4 | `-H1-H2+H3` |
| 119, 125 | -1 | `-2H1+H2` |
| 125, 128 | -1 | `-2H1+H2` |
| 128, 139 | -1 | `-2H1+H2` |

There are eight changed scalar blocks. Each certificate gives the full
before and after lists. The new contact graph permits four fewer local
native calls. The transport builder includes both meeting charts and all
remaining scalar contacts. It does not assign free movement to a basis
operation. All public paths remain nested.

The separate endpoint checks cover 743 arbitrary public inputs and nine
owned register slots. There are eleven owned lifetimes and peak occupancy
two. Every owned instruction and transport stays in the same order. All
535 productive helper targets and all public boundary charts stay fixed.

## 3. Native and caller binding

The retained compiler lowers each chart change to a signed Gaussian-dyadic
recipe. A separate interpreter checks every signed Pauli image and the
global phase. It checks 2,918 distinct maps and 208 main endpoint adapters.
The helper normalizers retain the exact productive target profile.

The p16/a7 caller is unchanged. Its physical binding checks all 7,488
contacts and all 53,248 central-square columns. It has 18 full forward
calls, 15 thin query shears and 3,120 paid whole-idle jobs. The 3,328 main
targets and 9,630 helper targets give 12,958 productive fields.

These exact finite identities use the retained factorwise tensor rule.
They do not add a new analytic or all-size admission assumption. The five
existing framework assumptions and helper endpoint contract remain.

## 4. Complete cost

For `H_r(t) = r^t` and `0 < t < 1`, the paid change per local call is

```text
delta P = -7 H1 + 2 H2 + H3,
delta T = 0.
```

Strict concavity gives `1+2^t>3^t` and `2>2^t`. Thus the complete signed
moment increases by `18 * (7 - 2*2^t - 3^t) > 0` throughout that interval.
The new bill has 66,402 paid calls and rank surplus 760. The local paid rank
mass remains 5,986; the helper rank mass remains 4,266. Their difference is
the 1,664 main progress plus 56 owned descent debt.

Two rational interval engines verify the full moment and put its unique
raw crossing in `(0.00310815850980, 0.00310815850981)`. The claim is
`0.00310815849`, and the tested saving is `0.003108158495`.

All scalar changes are priced. The local scalar allowance is 220,290 and
the route allowance is 5,167,000. The retained complete raw ordinary bound
is 184,647,552, below its 1,000,000,000 budget. The power-1024 compiler,
cutoff 17, floor 18, tail parameter 1/5, reserve 10^-12 and decay 39/40 are
unchanged. All 54 strict and four functional outer conditions pass.

## 5. Verification limits

The full source archive is authenticated before any retained module is
loaded. The nested predecessor archive is authenticated again before its
local or native checker is loaded. The new checker verifies the source
word and every intermediate complete word. Its separate dense transpose
engine checks 565,504 entries.

False basis signs, false block coefficients, a missing block certificate,
a missing native transport, a missing owned erasure, an unpaid query and
extra target credit are rejected. The inherited phase, caller and price
controls also run.

The huge uniform word and rounded multiplier are algorithm-defined and
are not materialized. These are author checks, not independent scientific
acceptance or a proof-kernel theorem. The numerical improvement is small.
No large gain or global optimum is established.
