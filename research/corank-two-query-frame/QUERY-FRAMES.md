# Two exact central-query constructions

These are query identities over characteristic zero. They do not assert a
complete stage endpoint, native caller or new exponent.

## 1. The rank-complement test

Write T_i=I+u_i v_i with v_i u_i=0. Let U contain the image columns, V
the covector rows, G=VU, and D=I minus the strict lower part of G. Direct
chronological expansion gives

    T_(L-1) ... T_0 = I + U D^(-1) V.

Suppose rank(G)=n and G factors through n fields. Then this product is -I
if and only if Z=D+G/2 has rank L-n.

For necessity put M=D^(-1)Z=I+D^(-1)G/2. The endpoint equality implies
M^2=M and UM=0. Its image is ker(U), of dimension L-n. For sufficiency,
M fixes ker(G), which has dimension L-n. The assumed rank of M makes
its image exactly ker(G)=ker(U). Thus UM=0 and
(UD^(-1)V+2I)U=0. Since U is onto, UD^(-1)V=-2I.

This is the same chronological Gram mechanism used in the shortest-query
normal form. The present construction uses complement rank two.

## 2. A rational rank-two complement

Let n be even and L=n+2. Choose nonzero t_i with

    product_(i even) t_i = - product_(i odd) t_i.

Set G_ii=0 and

    G_ij=sign(j-i)*(t_i+t_j)/t_j.

Then Z_ij=(1+t_i/t_j)/2, including the diagonal. Its rank is two unless
all t_i are equal. The skew matrix H=G diag(t) has Pfaffian

    pf(H)=2^(L/2-1)*(product_even t + product_odd t).

This formula follows by expansion along the first row and induction; the
intermediate products cancel in pairs. If the first n parameters are
positive rational numbers, the leading n-by-n Pfaffian is nonzero. The
full Pfaffian vanishes, and skew rank is even, so rank(G)=n.

Let K be that leading block, V=G[:,0:n], and U=K^(-1)G[0:n,:]. Then VU=G.
Its columns and rows give the required queries by the rank-complement test.
The code separately checks the complete literal product.

The integer family uses t_i=i+1 up to the last index; the dyadic family uses
t_i=2^i. The last parameter is fixed by the displayed product condition.
This is a family, not a classification of all rank-two complements.

## 3. An asymmetric two-to-three replacement

Let A=I+u_a v_a and B=I+u_b v_b be adjacent queries. Choose a nonzero vector
w in ker(v_a) intersect ker(v_b). Put

    A'=I+(u_a+w)v_a,       B'=I+(u_b+w)v_b,
    r=(1+v_b u_a)v_a+v_b,
    c=-r A^(-1) B^(-1),  C=I+w c.

Every new update is square zero, and the exact chronological identity is

    C B' A' = B A.

Indeed B'A'=BA+w r. Also BA fixes w and r w=0. Multiplication by C
therefore cancels the added term, with no residual higher product.

For an adjacent pair in the canonical shortest tuple, a nonzero w in the
common covector kernel is independent of u_a and u_b. The three new image
vectors have rank three. All three covectors lie in the original rank-two
covector span. They therefore cannot all be symplectic transvections for
one nondegenerate alternating form: that form would give an invertible
map from the image directions to their proportional covector directions,
which would preserve their span rank.

This replaces n+1 queries by n+2 while retaining their complete -I product.
The stage projectors are separate unknowns. The saved 12-field instance
passes only necessary low-order stage equations, and fails the completed
numerical searches. No full caller follows from the replacement identity.

The earlier shortest-query Gram and complete projector charts retain their
PR6168 attribution. The new rational family, asymmetric replacement and
checks are by Ali Romman with substantial OpenAI Codex assistance. No
independent review or worldwide priority claim. Apache-2.0.
