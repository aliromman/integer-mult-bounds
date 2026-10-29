# Alternating caller continuation

No kappa improvement or complete positive caller has been found.

[BOUNDARY.md](BOUNDARY.md) gives an all-length obstruction for the fixed
alternating interior with positive boundary stages. It uses cyclic trace
certificates at n=4,6,8,10 and a uniform rank proof for n>=12. Matrix port
coupling is allowed. [PARITY.md](PARITY.md) gives a separate all-length
stage eigenvalue bound for shortest query words, including the extra-stage
rank charge and the necessary one-full-call exception.

The fixed family is closed at the selected exponent, including zero, one
or two full boundaries with either orientation. Equal-sign small cases use
trace certificates, mixed-sign cases use a pole-rank price bound, and the
remaining short cases fail a favorable balanced-power price bound. The
source profile is authenticated against the retained selected PR6030
certificate. No native, compiler or outer replay is claimed here.

Normal and optimized trace-witness replays pass, including literal 2-by-2
and 3-by-3 matrix traces and changed-coefficient rejection. The curvature
block and rational price checks also pass. The proof is author-side work.
Independent scientific review is pending.
The stored trace witnesses can be replayed without a search:

    python3 -B research/alternating-history-caller/trace_verify.py

The checker imports the adjacent proof and initialized-query chart files.
It uses the standard library. SymPy 1.14 is used only for the separate
history, literal parity and scalar Groebner controls. Numerical boundary
searches use NumPy 2.5.4 and SciPy 1.18.1; their finite differences and inverse
cancellation controls are recorded. These searches are not no-caller proofs.

The complete saved prefix span is tested in powers of z. The earlier z-1
coordinate error in PR6168 is not repeated. For six fields with both end
stages removed, the history contains one degree-zero query direction and
two degree-one query directions. It does not give six independent targets.
Longer order-six histories in the retained tests do not repair this defect.

The initial ten-field inverse-boundary profile put two targets at power -2.
It fails the exact leading-coefficient rank cut in BOUNDARY.md. The revised
profile {-2:1,-1:2,0:1,1:3,2:2,3:1} passes that cut, but all eight numerical
starts have a nonzero residual. No full native or cost claim follows.

## Open next questions

Can arbitrary initialized-work query planes or the faithful order-six stage
escape these fixed-interior obstructions? The adjacent free-query-work
archive holds those tests. A positive construction still needs a complete
Laurent endpoint, actual native forward and inverse interfaces, independent
helpers, every owned lifetime, scalar work, routing, precision, compiler and
outer accounting. The selected frontier is unchanged.

Canonical investigation: Discussion5593, comment19199752. Base: main
3f0152a32eeebac0df11946dbd50cabe5a7f7839, priority revision8. PR6168 contains
the separately registered shortest-query proof. The source-edit lane in
Discussion5600 remains separate. Credit Cinna Davis and all retained
P104/native/caller/compiler/BIT authors. New work by Ali Romman with
substantial OpenAI Codex assistance. No subagents. Apache-2.0.
