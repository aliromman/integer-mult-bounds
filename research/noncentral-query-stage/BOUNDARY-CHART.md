# A rational query chart with an ordinary boundary shear

This note gives an exact chart and a complete finite screen inside it.
It does not claim a classification of every native caller or ordinary
coefficient algebra. All variables here commute and have characteristic
zero. There are no stage projectors in this calculation.

## 1. Standard-block endpoint

Use J=diag(1,-1) and S=[[-1/2,3/2],[1/2,1/2]]. At each site choose vectors
u_i^+,u_i^- and rows v_i^+,v_i^- with v_i^+u_i^+=v_i^-u_i^-=0. The standard
query block is

    T_i(alpha,beta)=delta_(alpha,beta) I
                     +(S_(alpha,beta)-delta_(alpha,beta)) u_i^alpha v_i^beta.

It is a rank-one unipotent. The trivial representation is I, and the sign
representation is I-2u_i^-v_i^-. These are the three representations of
one query conjugated by an ordinary bank frame in the J algebra.

Let L=n+1, with n even. The tested standard endpoint is

    diag(I,B) (S tensor I),     B=I+N, N^2=0, rank(N)=1.

The sign endpoint must independently be -B. For this normalized endpoint,
rank(M_standard-I)=n+rank(B-I). Thus n+1 rank-one queries permit at most a
rank-one boundary difference. The determinant then makes that difference
square zero. This is a rank test for this endpoint, not a native price.

## 2. A constructive rational chart

Choose nonzero t_0,...,t_n with product_even t=product_odd t. Put

    G_ii=0, G_ij=-1 for i<j, G_ij=t_j/t_i for i>j.

The matrix diag(t)G is skew. Its leading n-by-n Pfaffian is, up to sign,
t_0 t_2 ... t_(n-2), so that block K is invertible. Put

    V=G[:,0:n], U=K^(-1)G[0:n,:],
    h_j=2(t_j-1), g_i=2(1-1/t_i),
    r=K^(-1)g[0:n], w=h[0:n],
    N=r w/2,
    U_minus=U+r(1+h/2), V_minus=V-1 w.

Here 1 denotes the all-ones column or row of the needed size. The product
condition gives Vr=g and wr=0. Every diagonal entry of VU and
V_minus U_minus is zero. The plus query vectors are the columns of U and
rows of V; the minus vectors are from U_minus and V_minus.

Direct chronological multiplication gives the entire standard endpoint
diag(I,I+N)(S tensor I). boundary_chart.py checks this with rational matrix
products. The sign block is not automatic. In particular, t=(2,3,4,5,15/8)
gives a nonzero sign error in every entry of the four-by-four sign block.

One route to this chart starts from the standard endpoint equations. If
D=I+(3/2)lower(VU)+(1/2)lower(V_minus U_minus), their top-left equation is
UD^(-1)V=I. The rank-one complement D-VU has diagonal one and can be scaled
to the all-ones matrix. The other blocks then give
V_minus=V-1w and U_minus=U+r(1+wU/2). The two sector nilpotence equations
give the displayed t chart. The saved finite screen uses the chart as
constructed; it does not rely on an unreviewed global completeness claim.

## 3. A rank-two sign obstruction

Let G_minus=V_minus U_minus and D_minus=I+2lower(G_minus). If the sign word
is -I-N, the matrix Z_minus=D_minus-G_minus has rank at most two. This
follows from the chronological product and

    rank(I_L-AB)=L-n+rank(I_n-BA).

In this chart, H=diag(t) Z_minus is the symmetric matrix

    H_ii=t_i,
    H_ij=2t_j-t_i for i<j,
    H_ij=2t_i-t_j for i>j.

For ordered indices with values a,b,c, its principal three-by-three
minor is -3(a-b)(a-4c)(b-c). For four ordered indices with values a,b,c,d,
the minor with rows1,2,3 and columns2,3,4 is
-6(a-b)(b-c)(c-d). These formulas are direct determinant identities.

Remove adjacent repeated values from t. Four remaining values would make
the second minor nonzero. Thus at most three constant runs remain. With
three runs the first minor requires first=4*last. With two runs, the
product condition sets the value on the odd-length run to one. With
three runs, either all lengths are odd, giving middle=first*last, or only
one length is odd, giving that run's value one.

At n=4 this leaves four two-run patterns and six three-run patterns.
run_screen.py substitutes all ten into the literal sign word. It takes
the monic gcd of the rational polynomial numerators of every residual.
The two-run gcds have only the roots0 and1, or only1. Zero is outside
the chart; one makes N=0. Every three-run gcd is1. Therefore this chart
has no nonzero boundary shear at n=4. The one-run case also has t=1 and
N=0. No numerical approximation or generic Groebner timeout is used in
this finite conclusion.

The matrices and minors give useful exact constraints for larger cases.
An all-length sign exclusion, a global classification and a full native
caller are not claimed. Prepared by Ali Romman with substantial OpenAI
Codex assistance; predecessor attribution is in README.md. Apache-2.0.
