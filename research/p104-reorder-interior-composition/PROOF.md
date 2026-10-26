# Exact composition and remaining split domain

## 1. Frozen sources

PR5678 at af814f6f3c477fdd13b8f04380a33cf32004966b supplies the four
paid interior bases. PR5666 at693a2203a4f84f8ff3d99daef07432004b717d06
supplies the helper473 reordering. Both start with PR5660 at
3b0b4f37d5005f381ed9cfcccebfbe1ab35c71aa. The verifier first rebuilds
all four interior changes from that source and compares the full word.
It also rebuilds the published PR5666 word from the same source.

The new constructor applies PR5666's seven exact, disjoint source runs
to PR5678's word. Each old run must occur exactly once. It changes only
helper473's initial chart, the specified events and the schema label.
The four interior bases remain in the executed word. No published decimal
gain is used to infer compatibility.

## 2. Reordering, frames and endpoints

The scalar and owned operations have the same multiset as PR5678.
For each register, two source reads commute, and two ADDs with that
register as destination commute. Every other shared-register relation
is exclusive. Owned operations are exclusive, and their lifetimes share
one ordering token. The verifier matches all occurrences and proves that
every dependent block retains its order.

The read `213 += (1/2)*473` moves to chart3, before the read from706 and
the associated chart change. Helper473 starts at3. Its update from100
moves to716. The adjacent native paths use the exact source edit list.
All scalar contacts have one actual common chart. Every public move is
nested in its required direction, and every owned move follows its
source lifetime and opposite nested direction.

Forward rational rows and a separate dense integer transpose check prove
the same complete endpoint on all743 independent inputs and all owned
columns. Native lowering checks every distinct chart change, each helper
target, the main adapters, every signed Pauli image and the global phase.
The frame induction gives the retained all-width statement. Finite
checks alone are not presented as a replacement for that statement.

## 3. Full independent-helper cut check

Let A be a scalar prefix and B its remaining suffix. Let H be one group
of helpers with equal initial and final charts. At a cut bank F in one
strict intermediate chart, let a be a functional supported on F and b
a direction supported on F. An independent split needs

    a A = c,    B b = d,

where c and d use only H, and cd=1 after normalization. The cut and
boundary pivots may be different physical roles. Their ordinary exchanges
do not remove this necessary independent-coordinate pairing.

The checker eliminates every input outside H from the bank's prefix
rows, retaining both the helper part c and the field coefficients a.
It separately eliminates every output outside H from the suffix columns,
retaining d and b. These two exact rational kernels span every possible
pure pair at the cut. If every basis pairing cd is zero, bilinearity
makes every pair zero. This excludes arbitrary rational bases, and thus
also dyadic bases. Empty owned coordinates do not supply an independent
helper input or output.

The forward pass records each affected bank after every operation.
A move affects both endpoint banks. An ADD, SCALE or COPY affects its
common bank. INIT and CLEAR affect the old and new banks and their exact
zero rows. The reverse pass applies the actual transpose operations;
COPY accumulates into its source and clears its destination column.
Every unchanged bank retains both matrices until one of those events
changes it. Thus the touched-bank schedule covers every operation cut.

All47 helper boundary groups are considered. Proper intermediate charts
produce30,096 checked contexts. Every helper pairing is zero. No restriction
requires the boundary helper itself to be present in the cut bank. The
result is about this emitted word and this independent common-frame
split contract. It does not exclude new paid transports or a changed word.

## 4. Complete bill

Per local invocation, the paid change is -2H1+H3. The independent target
change is -H8+H9. Across all18 calls, these become -36H1+18H3 and
-18H8+18H9. The two rank-mass changes are both18, so the complete surplus
remains760. The full bill has66,384 paid calls and12,958 productive targets.
Thin queries, all whole idles, main targets and physical boundary programs
are included without change to their source contracts.

The local scalar allowance remains220,290. All entry and exit native
maps, target normalizers, main adapters, owned lifetimes, internal routes,
parent parking, precision, finite stock and compiler work remain charged.
The complete ordinary-work allowance stays below10^9. It is not a claim
that any ordinary primitive has zero cost.

The price verifier uses two rational interval engines and checks every
compiler and outer condition. It tests.003108196375 and admits
conditional kappa.00310819637. The exact root interval, signed margins,
finite thresholds, controls and all58 outer conditions are in RESULT.json.
The margin above PR5678 comes from the executed new invoice. More digits
of the same crossing are not counted as a separate result.

## 5. Verification and limitations

Every retained file is hashed before import. The source export loaders
check paths, types, uniqueness, exact inventories and each payload digest.
The declared pricing edits are checked byte-for-byte against the retained
parent code. Normal and optimized replay must reproduce the same word
and complete JSON result. Explicit exceptions protect the new construction
checks under optimization; the original R22 assert-based constructor is
retained for provenance and its edit data, not used as an acceptance gate.

The source framework assumptions and free helper endpoint contract remain.
The enormous compiled word and rounded multiplier are not materialized.
There is no claim of independent scientific acceptance, a global optimum,
or a selected-frontier change.
