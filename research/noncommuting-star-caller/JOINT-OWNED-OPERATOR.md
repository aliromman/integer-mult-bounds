# The complete eight-input owned extension

This calculation uses the pinned PR5599 word. It does not use the newer
selected PR6030 word. All eight former owned inputs are now arbitrary.
Replace COPY by ADD and omit INIT and CLEAR. Keep every MOVE, ADD and SCALE.
The full scalar endpoint on 756 inputs is

    J = [T B; C D].

Here T is the original public shear on 748 fields. Exact rational forward
and reverse-transpose interpreters agree on all 571536 coefficients. The
matrices B and C have rank eight. The complete Schur complement is

    D - C T^-1 B = D
      = diag(1/2,1/2,-1/2,-1/2,1/4,1/4,1/4,-1/4).

Thus the joint operator is invertible. Every new input affects old main
outputs. None is an independent new productive output under the old ABI.
The retained JSON gives every coefficient and literal boundary chart.

## A concrete completion and its price

The source has two useful common-frame properties. Every nonzero row of
T^-1 B has initial frame zero. Every column used by C T^-1 has final frame
equal to the full nine-dimensional space. All eight new entrance frames
have dimension eight; all eight new exit lines lie in their entrance spaces.

The block factorization is

    J = [I 0; C T^-1 I] diag(T,D) [I T^-1 B; 0 I].

Apply the inverse right shear before J and the inverse left shear after J.
This gives exactly diag(T,D). In physical coordinates, each new input
needs one rank-eight inverse Gaussian to supply its zero-frame preparation.
Each new output needs one rank-eight inverse Gaussian from the common full
frame to its exit line. Scalar gathering, addition, scaling, temporary
copies, cleanup and routing remain positive ordinary costs. These are a
specified implementation, not an optimality claim for all completions.

After the dyadic diagonal normalization, the eight added outputs have
Gaussian rank seven. Relative to the unchanged old public interface, the
recursive change is therefore

    added paid = 16 H8,       added targets = 8 H7.

The linear margin changes by -72. For every positive nondecreasing tariff
g, the change 8g(7)-16g(8) is strictly negative. Thus this direct completion
cannot improve kappa, even before its ordinary costs are included.

The chart induction extends the scalar endpoint to all tensor widths under
the retained native all-column premise. This package does not emit the
complete physical completion word or validate its compiler. The losing
recursive price is the stop test for this implementation.

## Open next questions

- Can a complete caller use the joint leakage instead of canceling both
  block shears? Such a caller must prove all old and new output jobs.
- Can another source produce common-frame leakage corrections already paid
  inside the caller? Reuse of a paid call must preserve its original input.

The earlier eight-entrance probe in this package and the released
Discussion5600 tests of PR5786 are predecessor evidence. Those tests keep
all old public outputs fixed. This calculation adds the complete joint
matrix, its exact block factorization and the price of one explicit repair.
Credit the pinned P104 sources in SOURCE-PINS.json and NOTICE.md.
Ali Romman with substantial OpenAI Codex assistance. No new kappa.
