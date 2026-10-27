# An exact two-call identity with a small residual

This is an approximate-output construction. It is not an exact Gaussian
caller or a kappa improvement. Its recursive precision cost remains open.
The two-full border is predecessor work: credit Dillon/PR4623 and Yiyang's
two-full Laurent border, paid-precision audit and bounded-lane continuation.
The retained local identity is a control in the present normalization. It
is not claimed as a new border method. Main contains the relevant audits in
research/two-full-border-paid-precision and research/bounded-lane-border-core.

Let U commute with S, with S^2=I. The advancing contract supplies both

    F_1 = U [[1,-S],[0,1]], F_2 = U [[1,0],[S,1]].

They use two unchanged full native words, with fresh independent helpers.
Let B=[[0,-1],[1,0]] and T_e=diag(1,e). For every nonzero scalar e,

    T_e B^-1 F_2 F_1 T_e^-1 = U^2 [[S,0],[-e,S]].

Thus the first output is U^2 S x. The second output is U^2 S y-e U^2 x.
The residual has rank one in the bank index. In Euclidean norm its operator
norm is |e|, because U is unitary. The input y is scaled by e^-1.
The exact Laurent check in border_caller.py retains this term. It rejects
both its omission and its wrong sign. It does not replace a limit by an
exact endpoint.

With the three-block advancer, U^2 S is ordinary equivalent to rank17.
If the residual could be removed without a recursive fee or a precision
loss, the selected source would have104 rank17 main targets per full call.
Its linear target-minus-paid margin would be48 per call. This is an
optimistic counterfactual price, not a supplied exact recursion.

## Precision boundary

Take e=2^-k. If ||x||_2<=2^b, the structural error is at most2^(b-k).
A requested error at most2^-p therefore permits the sufficient choice
k>=b+p. The scaled input can then require k more integer bits. Native
approximation errors also need a full propagation proof through both calls.

One cannot conclude that this cost is a fixed constant over a recursive
tree. At a descendant, the larger input magnitude can affect the next k.
Repeated scaling can increase the required precision. An exact scalar
identity at each finite e does not bound that increase.

There is also an exact repair: apply S to the first output and add e times
that value to the second. It costs one extra central query per port. With
the selected source, the linear margin then changes from96 to-8 for the
two-call block. That repair does not improve kappa.

## Open next questions

- Can a coupled return contract keep the residual correlated through recursion?
- Can a source-bound precision analysis control the scaling without consuming
  the complete moment margin?
- Can an exact larger word cancel the residual with already-paid work?

No border-rank theorem is assumed to answer these questions. A theorem about
an exponent of n does not automatically preserve an exponent of log n.
