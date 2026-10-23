# Whole ordered Schur residual for arbitrary H

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.

## Contract and conclusion

Work over any field F. Let P0,Q0 be commuting idempotents on V, with ordered
joint sectors V00,V10,V01,V11 of dimensions r00,r10,r01,r11. Let U=F^d,
d>=1, and let i:U->V, pi:V->U satisfy

    pi i = pi P0 i = pi Q0 i = I,       H = pi P0 Q0 i.

S,T are the common-identity-pivot Schur complements on ker pi. Require
rank(ST-TS)=2d. The all-field sharp-sector criterion from Nick Petalas,
as stated in the pinned nilpotent source, is a cited premise: the sector
factors X_j=pi|V_j, Y_j=E_j i have X10,X01,X11 onto and Y10,Y01,Y11
injective. Their products are -W,W,W,H, where W=I-H. Sector00 has no extra
rank restriction. The factor-extension principle of MatShep and codyw912
is retained and is implemented explicitly below.

Set a0=nullity H and a1=nullity W. A sharp split exists exactly when

    r00 >= d-a1,
    r10 >= d+a1,  r01 >= d+a1,  r11 >= d+a0.                 (1)

For every H and every such profile, all admitted ordered residual pairs
are simultaneously similar. The representative below determines that entire
pair. No splitting field, nilpotent-H assumption, or odd-characteristic
assumption is needed. Different H/profile data can still give the same
residual orbit; we do not claim that these parameters are a unique orbit
label. Counts of headers and physical admission are separate questions.

## Representative

Take the Fitting decomposition for W:

    U = U_star direct-sum U_one,
    U_star = im W^d,   U_one = ker W^d,
    d_star=dim U_star,  d_one=dim U_one.

Write H_star,W_star and H_one,W_one for the restrictions. W_star is
invertible; W_one is nilpotent; H_one=I-W_one is invertible. Choose columns
N0 that complement im H_star in U_star and N1 that complement im W_one
in U_one. They have a0 and a1 columns. Empty components are omitted.

On U_star^3 direct-sum F^a0, use coordinates (x,y,z,t) and the ordered pair

    S_star(x,y,z,t) = (H_star(x-z)-N0 t, 0,
                       W_star(z-x)-N0 t, t),
    T_star(x,y,z,t) = (0, H_star(y-z)-N0 t,
                       W_star(z-y)-N0 t, t).               (2)

This is the predecessor three-copy formula. Its derivation needs only
invertibility of W_star, not nilpotence of H_star.

On E direct-sum E, E=U_one direct-sum F^a1, put

    A = [ W_one N1 ],
        [   0    0 ]
    S_one = [ I A ],       T_one = [ 0 0 ].                 (3)
            [ 0 0 ]                [ A I ]

The complete representative is (2) direct-sum (3), with scalar pairs

    (0,0)^(r00-d_star),
    (1,0)^(r10-d-a1),
    (0,1)^(r01-d-a1),
    (1,1)^(r11-d-a0).                                     (4)

All exponents are nonnegative by (1). Its dimension is sum r_j-d. The
choices of component bases and N0,N1 change the representative only by
simultaneous similarity. This follows from the constructive reduction of
any original split to the representative, not from matching central spectra.

## Existence and full-rank normalization

A factorization XY=B with X onto and Y injective requires ambient size
at least d+nullity B. The intersection im Y intersect ker X is Y ker B.
The construction X=[B,N,0], Y=[I;0;0], with N a complement to im B,
attains every size above that bound. An unrestricted factorization of -W
exists in dimension r00 exactly when r00>=rank W=d-a1. These independent
sector factors satisfy all three split/pivot equations, since their products
are -W,W,W,H. The cited sharp criterion proves necessity and sufficiency
of (1).

Here is a concrete basis for a full-rank factor pair. Start with the d
columns of Y. For each column v of a chosen complement N to im B, append
one solution l of Xl=v. These columns are independent: a dependence maps
under X to a relation between im B and the complement. Complete them with
vectors from ker X. This works because im Y+the lifts+ker X spans the
whole sector. In that basis X=[B,N,0] and Y=[I;0;0]. Thus normalizing the
three full-rank sectors introduces no extra extension data.

Use the Fitting basis of U first. W and H are then block diagonal. Choose
the W complement wholly in U_one and the H complement wholly in U_star.
The construction above therefore splits sectors10,01,11 into their two
component factor systems and unused scalar directions.

## The missing sector00 gluing step

Only sector00 lacks full ranks. Put X=X00,Y=Y00. Let e_star be projection
onto U_star along U_one, and let B=-W. The operator

    J = Y (B|U_star)^(-1) e_star X                         (5)

on V00 is an idempotent, where the inverse has domain U_star. Indeed XY=B,
so J^2=J. Moreover

    im J = Y U_star,  XJ=e_star X,  JY=Y e_star.

Consequently V00=Y U_star direct-sum ker(e_star X). The first summand
has dimension d_star and carries the invertible product -W_star. In the
second summand the factors map only to and from U_one, with product -W_one.
This is a decomposition of the actual original factor system, including
arbitrary cross-component and rank-deficient data. It does not assume a
preexisting direct sum of headers.

An explicit sector basis is the columns Y|U_star followed by any basis of
ker(e_star X). In those coordinates the first component has Y=I and
X=-W_star; all cross terms vanish by the displayed identities. The second
component has dimension r00-d_star and retains arbitrary factors X1,Y1
with X1Y1=-W_one.

Combine this sector basis with the three preceding full-rank sector bases.
The ambient space, i, pi and both original projectors split into a star
part, a one part and unused scalar sectors. Therefore ker pi and both
Schur residuals split too. This is the justified gluing implication.

## The explicit residual maps on both components

For the star part use coordinates (x,y,z,t) in the kernel through

    (v00,v10,v01,v11,extra11)
      = (x+y+W_star^-1(H_star z+N0 t), x,y,z,t).            (6)

Substitution into P0-P0 i pi P0 and Q0-Q0 i pi Q0 gives (2), exactly as
in the nilpotent predecessor. Formula (6) uses no inverse of H_star or
2I-H_star. In particular H_star=2I in odd characteristic is legal.

For the one part, write t in its full remaining V00, and let b,c be the
two a1 extra coordinates. Put C=(I+W_one)^(-1); it exists because W_one
is nilpotent, in every characteristic. The explicit kernel map is

    x = x_bar - C X1 t,    y = y_bar - C X1 t,
    z = -H_one^-1(X1 t+W_one x+N1 b+W_one y+N1 c).          (7)

Its original sector coordinates are (t,x,b,y,c,z). As in the one-primary
predecessor, substitution gives (3) on (x_bar,b,y_bar,c), while every t
coordinate is killed by both operators. Y1 disappears. This proves the
first scalar exponent in (4), including inequivalent original sector00
factor systems. The other three scalar terms are the unused directions
from full-rank normalization.

Equations (5)--(7), preceded by the sector basis changes, are an explicit
injective map from the declared normal coordinate space onto ker pi.
Their parameter count is sum r_j-d, so the map is invertible onto that
space. They intertwine both S and T. This proves the whole ordered-pair
claim and the independence from all remaining factor choices.

## Evidence and limits

normal.py constructs that map using exact prime-field elimination. verify.py
independently builds the original commuting projectors and full headers,
forms literal Schur complements, and checks both intertwining equations and
the invertibility of the map. It does not infer conjugacy from central
invariants. Complete domains cover every H in dimensions1,2,3 over F2,
dimensions1,2 over F3 and F5:1244 matrices. There are90 further transported,
padded factor systems and12 larger mixed-component fixtures, including an
irreducible quadratic block over F2. Empty00, H=0, H=I, higher primary blocks,
characteristic two and H=2 in odd characteristic are included.

The theorem holds over every field under the cited sharp-sector premise;
the code executes prime fields only. Original sector basis changes, kernel
choices and the residual matrices are explicit. Runtime and storage for
physical header emission are not priced. Nonsharp states, header/orbit
probabilities, later guard survival, native operations, precision and the
complete supplier/outer calculation remain outside the result. No complete
kappa improvement, independent human review or proof-kernel acceptance follows.
