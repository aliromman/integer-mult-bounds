# Exact joint-nullity law and polynomial Schur header counting

Ali Romman with substantial OpenAI Codex assistance. Apache-2.0.

## Result and inherited contract

Let q be a prime power, d>=0, and let H range over all matrices in M_d(F_q).
Let N_d(a,b) count matrices with nullity H=a and nullity(I-H)=b. We give an
exact triangular recurrence with O(d^4) summands. It replaces the q^(d^2)
matrix enumeration in the retained binary Schur header count. No matrices
or field elements are enumerated by the new counting algorithm. All arithmetic
is integer arithmetic. The field order is a supplied prime power; the code
does not prove that an arbitrary integer q is a field order.

The application counts complete headers for a fixed commuting idempotent
pair and one common identity pivot, including sharp commutator rank2d.
It retains the source's actual joint dependence through H. It does not
multiply marginal probabilities, assume later guards pass, or give an
unconditional residual-orbit distribution. The complete kappa is unchanged.

The binary factor counts and header law belong to ssafwann's PR4678.
Their field-independent counting arguments extend to F_q as shown below.
The all-field sharp-sector criterion is retained from the pinned nilpotent
source's credited statement of Nick Petalas's result. Gaussian binomial
counting and incidence inversion are classical. The contribution is the
explicit joint law, its source-bound use instead of matrix enumeration,
and exact evidence across prime and extension fields.

## Incidence moments

Write [n choose k]_q for the Gaussian binomial coefficient, zero outside
0<=k<=n. Put

    T_d(i,j) = [d choose i]_q [d-i choose j]_q
               q^(ij+d(d-i-j)),              i+j<=d.       (1)

For a fixed H, the number of pairs (A,B) with A an i-subspace of ker H
and B a j-subspace of ker(I-H) is

    [nullity H choose i]_q [nullity(I-H) choose j]_q.

The two eigenspaces have zero intersection: Hv=0 and Hv=v imply v=0.
Thus A and B must be disjoint. Choose A in [d choose i]_q ways. A j-subspace
B disjoint from A is the graph of a map from a j-subspace of F_q^d/A to A.
There are [d-i choose j]_q choices of the quotient subspace and q^(ij)
choices of its graph. For such A,B, require H|A=0 and H|B=I. A complement
of A+B has dimension d-i-j and can map to F_q^d arbitrarily. This gives
q^(d(d-i-j)) extensions. These choices prove (1) and the moment identity

    T_d(i,j) = sum_{a>=i,b>=j,a+b<=d}
                  [a choose i]_q [b choose j]_q N_d(a,b).  (2)

## Exact finite inversion

Process a=d,d-1,...,0 and, for each a, b=d-a,d-a-1,...,0. Then

    N_d(a,b) = T_d(a,b)
       - sum_{i>=a,j>=b,i+j<=d,(i,j)!=(a,b)}
             [i choose a]_q [j choose b]_q N_d(i,j).       (3)

Every term on the right has already been computed. Its diagonal coefficient
in (2) is one, so this is the unique solution. The actual matrix counts satisfy
(2); hence (3) equals those counts. This proves integrality and nonnegativity
without an approximation or a finite-dimension extrapolation. In particular
sum N_d(a,b)=q^(d^2), N_d(d,0)=N_d(0,d)=1, and N_d(a,b)=N_d(b,a) by H->I-H.
Terms with a+b>d are zero, including in characteristic two.

There are (d+1)(d+2)/2 strata and O(d^4) terms in all sums. Gaussian
coefficients can be prepared in polynomial work. The integers in the joint
law have O(d^2 log q) bits: they count matrices with chosen subspaces, and the
number of subspace pairs is at most q^O(d^2). Ordinary exact integer powers,
products and division therefore give polynomial bit work too. This is not
an O(d^4) bound on individual bit operations. The implementation uses Python
integer powers and exact divisions; their work is not declared free.

## The unchanged header ensemble, over any finite field

Fix the ordered original joint-sector dimensions r=(r00,r10,r01,r11), with
n=sum r_j, and 1<=d<=n. A complete ordered header is in GL_n(F_q). The event
requires both leading d-by-d blocks to be I_d. A split i,pi with pi i=I
has sector factors with products

    X00Y00=-(I-H), X10Y10=I-H, X01Y01=I-H, X11Y11=H.       (4)

These four equations already imply the split and both pivot equations.
The cited sharp-sector criterion requires X10,X01,X11 onto and Y10,Y01,Y11
injective. Sector00 remains unrestricted.

Define I_q(r,k)=product_{l=0}^{k-1}(q^r-q^l), zero for k>r, and G_q(r)=I_q(r,r).
For a fixed d-by-d matrix of rank rho, its factor counts through F_q^r are

    U_q(r,d,rho) = sum_{k=rho}^{min(d,r)}
       [d-rho choose k-rho]_q I_q(r,k) q^(d(r-k)),          (5)

    A_q(r,d,rho) = I_q(r,d) q^(rho(r-d)) I_q(r-d,d-rho)     (6)

when r>=d, and A_q=0 otherwise. Formula (5) allows all X,Y. To prove it,
choose the image of X, a k-space containing the rho-dimensional image of
the product, in [d-rho choose k-rho]_q ways. There are I_q(r,k) onto X
with that image. Each of the d columns of Y has q^(r-k) preimages under X.
For (6), choose injective Y first. X is then fixed on im Y. On a complement
of dimension r-d, its quotient map to F_q^d/im(product) must be onto,
giving I_q(r-d,d-rho) choices; lifts contribute q^(rho(r-d)). These are
the credited binary source arguments with q in place of two, not a new
factor-count principle.

For a=nullity H, b=nullity(I-H), put rho0=rho1=rho2=d-b, rho3=d-a. The exact
numbers of common-pivot and sharp complete headers are respectively

    Z = G_q(n-d) sum_{a+b<=d} N_d(a,b) product_{j=0}^3 U_q(r_j,d,rho_j),
    C = G_q(n-d) sum_{a+b<=d} N_d(a,b) U_q(r00,d,d-b)
                        A_q(r10,d,d-b) A_q(r01,d,d-b) A_q(r11,d,d-a). (7)

Each split has exactly G_q(n-d) ordered bases of ker pi. Therefore C/G_q(n)
is the unconditional one-block sharp probability. C/Z is the conditional
probability given both identity pivots only when Z>0. No quotient by GL_d,
a sector centralizer or an H similarity-class size belongs in (7). The
matrices H are labeled leading-coordinate matrices, and N counts each once.

Inserting (3) into (7) gives a polynomial exact algorithm. Reuse each factor
count for its sector and rank. Header integers and intermediates have
O(n^2 log q) bits, so their arithmetic is also polynomial in n and log q.
The result computes one-block events. Applying it again to a noncommuting
residual without a new theorem is invalid. Nullities alone need not label
the whole residual orbit.

## Evidence and discriminators

The independent field implementation enumerates25,291 matrices across F2,
F3,F5,F4 and F8. F4 and F8 use polynomial-basis arithmetic with irreducible
polynomials x^2+x+1 and x^3+x+1. Full elimination gives both nullities for
every matrix. All tables match (3), including zero strata.

Separate complete factor domains contain11,527 (X,Y) pairs and382 prescribed
product targets. They check both (5) and (6). Five source binary header
profiles match the unchanged predecessor predict function, including the
sharp densities3/496,9/496 and27/2763520. All20,160 literal GL4(F2) headers
are compared through the source's independent full-matrix Schur oracle.
At q3,d1, profile(1,1,1,1), the unique eligible H=2 contributes16 splits;
this is an odd-field calibration, not the binary scalar law.

The d32,q2 case replaces a space of2^1024 matrices with561 strata. The
d16,q4 case uses153 strata for4^256 matrices. Exact complete-header counts
are retained by bit lengths and SHA256 of canonical unsigned big-endian
integers, rather than printing thousands of decimal digits. These are
executed calculations, not estimates from small cases.

No production header sampler, later-guard chain, native operation, source
placement, precision, compiler, supplier, or outer multiplication proof is
supplied. Finite checks corroborate the quantified incidence proof. No
independent human review or proof-kernel acceptance is claimed.
