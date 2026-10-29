# Arbitrary query planes with owned work

These are discovery searches. No small floating residual is an exact endpoint.
No new kappa is claimed.

The input has r arbitrary fields and w zero fields. Each of r+1 full calls
has the form I+U A(z,S)V, with VU=I_2 and
A=(z-1)I_2+z(S-1)E_01. The retained r fields must equal
S times diag(z^2,...,z^2,z^3,z^3). Discarded fields have no endpoint condition.
Both query planes and stage projectors can change. All complex parameters
are free within this factorization. The numerical equation is evaluated at
r+2 distinct roots, so a zero residual would specify the full formal
polynomial. The optional period flag instead tests all p-th roots. It
changes the contract to the quotient by z^p-1.

The normalized option replaces V by (VU)^(-1)V at each evaluation. This
imposes the projector constraint exactly within the invertible chart. A
singular chart is outside this parameterization. It does not prove that all
solutions lie in a bounded numerical region.

The polynomial-chart option uses U=[I;A] and V=[I-BA,B], in a fixed
canonical frame. The query core is [[-xy,x^2],[-y^2,xy]], which squares
to zero. Thus VU=I exactly, with no variable matrix inverse. This chart
covers planes whose chosen two-row minor is invertible. The zero query
at x=y=0 is a degenerate boundary and is not an admitted full call. Any
candidate must have nonzero query cores and pass an exact reconstruction.

The first normalized period-six run completed one trial and then failed
because SVD did not converge. Its failure receipt is retained. This is
neither a caller nor an impossibility result. The polynomial chart avoids
the variable inverse that can make that parameterization singular.

The analytic Jacobian is checked against finite differences. The exact
three-call triangle is a positive control. Each receipt lists the seed,
iteration limit, numerical residual and elapsed time. Python 3.12 with
NumPy 2.5.4 and SciPy 1.18.1 was used for discovery. Set
OPENBLAS_NUM_THREADS=1 when running search.py.

The four-field one-work and two-work formal searches, and the six-field
two-work formal search, did not reach an exact candidate. Their nonzero
residuals do not prove impossibility. Large fitted coefficients can indicate
a border limit; they do not establish one.

## Open next questions

Can an order-six stage produce a finite exact caller on eight fields? Test
all characters, then reconstruct exact coefficients and check the complete
native operator, independent helper contract, ordinary maps, ownership,
precision, compiler and outer bill. Every discarded field has a real
initialization and erasure cost. No cost is inferred from these searches.

This continues Discussion5593 from canonical main
3f0152a32eeebac0df11946dbd50cabe5a7f7839. The shortest-query obstruction in
PR6168 is a separate proof. The source and query lineage credits Cinna Davis
and all retained P104/native/caller contributors. New discovery: Ali Romman
with substantial OpenAI Codex assistance. No subagents. Apache-2.0.
