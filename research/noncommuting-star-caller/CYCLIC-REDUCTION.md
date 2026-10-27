# When a cyclic stage gives no new caller identity

This lemma applies over every field. It does not reduce complex coefficients
modulo two. It gives a necessary test, not a new kappa.

Let z be a central indeterminate. A full factor has the form

    F_i(z) = C_i + z L_i,

where the images of C_i and L_i are disjoint, and their kernels have zero
intersection. The native factors meet these conditions. For a forward stage,

    C_i = I-P_i, L_i = P_i(I+c_i N_i),

where P_i is a projector, N_i is square-zero, and P_i N_i=N_i P_i=N_i.
For an inverse stage, multiply its entire matrix by z; this interchanges
C_i and L_i. The same image and kernel conditions hold. Fixed invertible
boundary matrices do not change the conclusion below.

## Lemma

Take a product M(z) of l such factors. Reduce it modulo z^p-1, with p>=2.
Suppose its reduced coefficients of degrees 0 through k-1 are zero, where
1<=k<p. If

    l <= p+k-1,

then M has degree less than p before reduction. In particular the reduced
identity is already an identity with a formal stage variable. The cyclic
relation supplies no extra cancellation in this range.

Proof. Induct on k. If k=1 and l<p, there is nothing to reduce. If l=p,
the constant reduced coefficient is M_0+M_p. The first term is in the image
of C_l; the second is in the image of L_l. These images are disjoint, so
both terms are zero. Thus the degree is less than p.

For k>=2, write M=(C+zL)B, and let B_j be the reduced coefficients of B.
For j=0,...,k-1, with indices read modulo p,

    C B_j + L B_(j-1) = 0.

Disjoint images give C B_j=0 and L B_(j-1)=0 separately. For
j=0,...,k-2, both C B_j and L B_j are zero. Their kernels have zero
intersection, so B_j=0. The product B has l-1<=p+(k-1)-1 factors.
The induction hypothesis gives degree B<p and zero formal coefficients
0 through k-2. Hence M has degree at most p and formal constant coefficient
zero. Its reduced constant coefficient is also zero, so its degree-p
coefficient is zero. This completes the induction.

The bound on l cannot be increased in this general class. A product of
p+k factors zI is z^(p+k) I. It reduces to z^k I but is not that formal
polynomial. This is a sharpness control, not a profitable caller.

## Determinant test for full native calls

There are n fields and l full calls. Let a calls use the inverse stage.
Each native query shear has determinant one. Multiplication by z on all
n fields at each inverse call makes the word polynomial. Its determinant
is

    z^(2l+(n-4)a).

Let e_j be the requested signed target powers. Put d_j=(e_j+a) mod p,
with 0<=d_j<p, and k=min d_j. If k>=1 and l<=p+k-1, the lemma gives
the necessary integer equality

    sum_j d_j = 2l+(n-4)a.

A congruence modulo p is not sufficient in that range. Independent constant
input and output bases can change a nonzero determinant scalar, but cannot
change this degree. The test applies separately at either central-query
character S=1 or S=-1, when those characters are distinct. It also applies
to the query-free channel of the characteristic-two tangent model.

For eight fields and nine forward calls with p=7, the target powers
(3,3,3,3,3,3,3,4) have k=3. The lemma applies, but their sum is25 rather
than18. This target is impossible in the stated model over every field.
The timed-out binary search of that target is therefore unnecessary for
this conclusion. The six-field/seven-call balanced test at p=7 is also
equivalent to its formal-stage test; this observation does not exclude it.

The eight-field/nine-call balanced degrees (2,2,2,2,2,2,3,3) lie just outside
the lemma's length bound. A signed test with four inverse calls and target
powers (3,3,-3,-3,-3,-3,-3,-3) has normalized degrees (0,0,1,1,1,1,1,1),
so k=0. The lemma does not decide either test.

## Checks and open next questions

`cyclic_reduction_check.py` checks all two-dimensional binary projector words
through length five. It tests the hypothesis for periods two through four
and retains a sharpness control. These finite checks support the proof;
they do not prove a native implementation or a complete caller.

- Search the first domains beyond the bound, rather than repeat an equivalent
  formal-stage search.
- Extend the test when the reduced constant coefficient is nonzero, including
  the caller with a removable ordinary input residual.
- Determine whether an endpoint beyond this bound has a characteristic-zero
  lift, actual native routes and a favorable complete invoice.

Prepared by Ali Romman with substantial OpenAI Codex assistance. The native
full-factor model and source caller lineage retain their attribution in
SOURCE-PINS.json and LAB.md. The proof uses elementary polynomial and linear
algebra. No worldwide priority or independent acceptance is claimed.
