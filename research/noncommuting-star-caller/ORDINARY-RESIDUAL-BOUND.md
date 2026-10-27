# An ordinary residue does not save a query in the formal forward model

Let the scalar field have any characteristic. Work in its algebra with a
central symbol S and relation S^2=1. Keep the two coefficients in this
algebra separate, also in characteristic two.

There are n independent fields and l full calls. Call i has a rank-one
square-zero query N_i and a projector P_i, with P_i N_i=N_i P_i=N_i.
The stage variable z is central and formal. The call is

    F_i(z,S) = I-P_i + z P_i(I+S N_i).

Suppose the complete word has the exact endpoint

    F_l ... F_1 = C + S T(z),

where C is an ordinary constant matrix, T(0)=0, and T(1) is invertible.
Then l>=n+1. The conclusion allows any rank for C. In particular, a
rank-one ordinary input residue cannot reduce the number of full calls to n.

Proof. At z=1 write the prefix after j queries as C_j+S B_j. Then

    C_j=C_(j-1)+N_j B_(j-1),
    B_j=B_(j-1)+N_j C_(j-1).

The image of B_j lies in the sum of the images of N_1,...,N_j, so
rank B_j<=j. The constant coefficient in z of the complete word is C.
Its image lies in ker P_l, hence N_l C=0. Square-zero N_l and the first
recursion give N_l C_(l-1)=N_l C_l=N_l C=0. The second recursion now
gives B_l=B_(l-1), of rank at most l-1. Since B_l=T(1) is invertible,
n<=l-1. This proves the claim.

The result also holds for independent constant input and output basis
changes: multiply the endpoint by their inverses before the same argument.
It does not apply to signed stage factors, a cyclic identity with actual
wraparound, noncentral stage-dependent frames, or an endpoint that projects
away extra initialized work fields. Those are separate contracts.

## Fixed chain used in the search

For even n, set N_i=e_i(e_(i-1)^T-e_(i+1)^T), omitting out-of-range terms.
The product of the n raw queries is C+S B, where

    C=(e_0+e_2+...+e_(n-2)) e_0^T.

The even row i of B is -e_(i+1)^T. Each odd row i<n-1 is
e_0^T-e_(i+2)^T, and the last row is e_0^T. Thus det B is a unit. The
chain attains full query rank with only one ordinary residue before the
stages are inserted. The theorem explains why this raw identity alone
does not give a formal forward caller with all target stages positive.

The retained searches in chain_caller_search.py test cyclic and signed
variants. None of the current receipts is an exact endpoint. This proof
does not convert a numerical failure or a timeout into a wider exclusion.

## Open next questions

- Can a finite cyclic stage remove the residue without losing target rank?
- Can signed stages give an exact endpoint with a complete favorable bill?
- Can initialized work fields and a final ordinary output projection change
  the bound, while every original input has its required independent output?

Ali Romman with substantial OpenAI Codex assistance. The rank-one query
method and source interface retain the predecessor credit in LAB.md and
SOURCE-PINS.json. No independent acceptance or new kappa is claimed.
