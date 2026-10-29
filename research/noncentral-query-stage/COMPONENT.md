# An ordinary noncentral phase and a finite endpoint test

This is a component identity. It does not give a complete caller or a new
kappa. The native source, ordinary-work budget and all 104 routes have not
been replayed with this extension.

## 1. The query at every positive width

Use the query C=((1+i)I+(1-i)X)/2 from the retained P104 interface. At width
f put Q=C^(tensor f), X_all=X^(tensor f), and

    G=(I+X_all)/2 - i(I-X_all)/2.

The predecessor ADVANCING-CONTRACT.md proves that G is ordinary. It has a
two-address implementation and G^2=X_all=Q^2. This G is a whole-width map;
it is not a tensor power of the width-one normalization.

Put d=2 floor(f/2), k=d/2, and define the ordinary translation

    L_f = X_last if f is 3 modulo 4, and I otherwise.
    A_f = G L_f,       S_f = Q A_f.

All these maps commute. Thus A_f^2=Q^2, S_f^2=I, and
Q^(-1) A_f^(-1)=S_f. The forward-source wrapper in the predecessor can use
A_f in place of G. It adds at most one address translation in each adapter.
The eight advancing axes are separate from this query coordinate, so they
commute with all the maps in this note.

In the Fourier basis, Q has eigenvalue i^w on a character of weight w.
The eigenvalue of QG is (-1)^(w(w-1)/2). Define

    q_f(y)=sum_(i<j) y_i y_j over F2.

The Fourier phase of S_f is q_f(y), with y_last added when f is 3 modulo 4.
This proves the claimed normalization, including its absolute phase.

## 2. An ordinary quadratic phase

On the first d bits of an address x, put

    j_f(x)=q_d(x) + ((k+1) modulo 2) * sum_i x_i,
    J_f |x> = (-1)^j_f(x) |x>.

The last bit is a spectator when f is odd. J_f is one diagonal sign map.
Its sign can be computed from the Hamming weight modulo four. It uses no
recursive Gaussian call. Its address and scalar work still need a complete
ordinary-work allowance in any proposed caller.

Then, for every positive f,

    S_f^2 = J_f^2 = (S_f J_f)^3 = I.

Usually S_f and J_f do not commute. At f=2 the literal matrices are

    J=diag(1,1,1,-1),
    S=(1/2) [[1,1,1,-1], [1,1,-1,1],
             [1,-1,1,1], [-1,1,1,1]].

## 3. Proof of the order-three relation

First let f=d=2k>0. The polar matrix of q_d is B=I+ee^T, where e is the
all-ones vector. B is alternating and B^2=I over F2. In particular it is
invertible. Set q_star(x)=q_d(Bx). Direct expansion gives

    q_star(x)=q_d(x)+((k+1) modulo 2) * sum_i x_i.

Thus J=diag((-1)^q_star). Let H_B be the normalized Fourier matrix with
kernel 2^(-k)(-1)^(x^T B y). The quadratic Gauss sum has the form
sum_x (-1)^q_d(x)=epsilon 2^k, with epsilon in {1,-1}. This follows by
squaring the sum and using the nondegenerate polar form. Translation of
the summation variable gives the full Fourier identity

    S = epsilon J H_B J.

The same translated sum gives

    (J H_B)^2 = epsilon H_B J,
    (J H_B)^3 = epsilon I.

Consequently (SJ)^3=(epsilon JH_B)^3=I. This argument retains the sign
epsilon; dropping it is not a valid native phase check.

For odd f=d+1, the radical of the polar form of q_f is the all-ones line.
Write y=(x_0+r,...,x_(d-1)+r,r). Then

    q_f(y)=q_d(x)+(k modulo 2)*r.

The translation L_f removes that last term. The associated ordinary
linear-address permutation conjugates S_f to S_d tensor I. It leaves
J_f unchanged because its first d physical coordinates are unchanged.
The even-width proof therefore applies. At f=1 both S_f and J_f are I.

## 4. A complete constant-size algebra test

These relations give an action of S3, the six-element symmetric group.
Over a characteristic-zero field its three irreducible representations
can be written over the rationals:

| Representation | J | S |
| --- | --- | --- |
| trivial | 1 | 1 |
| sign | -1 | -1 |
| standard | diag(1,-1) | [[-1/2,3/2],[1/2,1/2]] |

Checking an identity in all three representations checks every coefficient
in the group algebra. Thus it is sufficient at every positive native
width. Coefficients can also be rational functions of a separate commuting
stage variable.

For even f=2k, the native characters are

    dim=4^k,  trace(S)=trace(J)=epsilon 2^k,  trace(SJ)=1.

The standard multiplicity is (4^k-1)/3. The sum of the trivial and sign
multiplicities is (4^k+2)/3, and their difference is epsilon 2^k. For k>=2
both one-dimensional multiplicities are positive because

    (4^k+2)/3 - 2^k = (2^k-1)(2^k-2)/3 > 0.

At odd width these characters have an extra spectator factor of two.
Therefore all three representations occur at every f>=4. The test is
also necessary for an identity required at every native width in this
fixed coefficient algebra. A one- or two-bit matrix is not a complete
all-width acceptance test: the sign representation is absent there.

For example, boundary_chart.py constructs an exact standard-block query
word at t=(2,3,4,5,15/8). Its sign block fails. This is a literal example
of the error that the three-representation test prevents.

## 5. What the new contract permits

Matrices with entries a+bJ are ordinary scalar bank maps selected by the
two signs of the address phase J. An invertible such matrix has an inverse
of the same form. Thus a full call can use different invertible bank bases
on the two J eigenspaces. Its compatible rank-two projectors can differ
between those sectors. These bases need not commute with S.

This changes the central-coefficient premise in PR6168. It does not change
the paid native source word or give a free extra Gaussian. Every conditional
bank operation, phase evaluation, adapter, route and owned lifetime must
be included in a completed caller. The source wrapper above is algebraic;
no complete 104-port execution has been emitted here.

## Checks and open next questions

component_verify.py uses literal integer kernels and full matrix products
at widths 1 through 8. It checks squares, cubes, absolute Fourier phases
and all character multiplicities. It rejects omission of the odd-width
translation at f=3 and omission of J at f=4. The written proof, rather
than extrapolation from those finite checks, gives the all-width claim.

Find a complete endpoint in this algebra, including its ordinary boundary
maps. Then bind the actual source routes and complete cost. Larger ordinary
coefficient algebras, extra queries and changed stage schedules remain open.
No positive caller or kappa increase has been established.

The Q-to-involution normalization and forward wrapper are from Ali
Romman's earlier noncommuting-star-caller checkpoint. The P104 query,
controller and source lineage retain their credit, including Cinna Davis.
The noncentral phase, all-width S3 test and present checks are by Ali
Romman with substantial OpenAI Codex assistance. No subagents, independent
scientific acceptance or worldwide priority claim. Apache-2.0.
