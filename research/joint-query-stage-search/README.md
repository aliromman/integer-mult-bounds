# Joint query and stage searches

No exact caller or new kappa was found. This archive lets future work replay
the searches, inspect their conditioning and avoid treating an approximate
border as a valid endpoint. The separate
[exact tariff package](../cyclic-stage-rank-budget/README.md) proves the rank
cuts; numerical failure in this archive is not a proof of impossibility.

For n independent fields, the program varies every n-by-2 matrix B_i and
2-by-n matrix C_i in n+2 full factors

    I + B_i [(z-1)I_2 + (s-1)z E_01] C_i,       C_i B_i=I_2.

It compares the complete product with s diag(z^d_j), for both s=1 and s=-1.
The default target has n-4 powers two and four powers three. Sampling at
length+1 distinct complex points determines the formal polynomial in exact
arithmetic. The period option instead checks each character of the specified
finite cyclic stage. Floating checks alone never certify either endpoint.

The soft constraint mode includes C_i B_i-I_2 in the residual. The exact-dual
mode uses (C_raw B)^(-1) C_raw. This makes C_i B_i=I_2 algebraically, subject
to floating roundoff. Complex parameters and an analytic real Jacobian are
supported. The controls check two exact triangles in floating arithmetic,
reject changed target powers and compare 24 Jacobian columns with finite
differences. The exact triangles come from the attributed PR6168 source.

## Retained evidence

All seven completed calls are in discovery-results. The conditioning file
reconstructs their effective projectors and query nilpotents. The best point
in each call was chosen by Euclidean residual norm, not maximum residual.

| Saved call | Fields / full calls | Stage | Best maximum residual | Largest query norm |
| --- | --- | --- | ---: | ---: |
| six | 6 / 8 | formal | .00274974 | 150.49 |
| six-long | 6 / 8 | formal | .00046365 | 891.53 |
| twelve | 12 / 14 | formal | .00541728 | 48.94 |
| six-complex | 6 / 8 | formal | .00294358 | 164.67 |
| twelve-period6 | 12 / 14 | period six | .00313617 | 244.62 |
| six-period7 | 6 / 8 | period seven | .29951899 | 122.84 |
| eight-period7 | 8 / 10 | period seven | .38195776 | 30.65 |

Residuals include the program's sample normalization. A bound on raw
parameters is not a bound on the effective projectors: the inverse of C_raw B
can be large. The small residuals, increasing query norms and small singular
values do not supply exact rational points or a uniform precision bound.
The known two-full border explains one observed limit. It is prior work by
Dillon in PR4623 and Yiyang in PR4966. Its paid precision and native binding
obligations remain; this archive makes no new border claim.

The period-seven targets are 3,4,4,4,4,4 and 3,3,3,3,3,4,4,4. Both pass the
determinant congruence and fail the exact cyclic-stage rank bound. The search
now rejects them before optimization. Use `--calibrate-excluded` only to
repeat a historical search; that flag does not weaken the theorem.

## Replay

The discovery environment was CPython 3.12.13, NumPy 2.5.4 and SciPy 1.18.1
on macOS arm64, with `OPENBLAS_NUM_THREADS=1`. Commands below run from the
repository root. Put the output outside the committed evidence directory.
All calls used seed 290610. Their saved parameters and metadata retain the
complete tested endpoints. Execution times are not portable.

```sh
export OPENBLAS_NUM_THREADS=1
python -B research/joint-query-stage-search/search.py --n 6 --trials 4 --steps 300 --output /tmp/six.json
python -B research/joint-query-stage-search/search.py --n 6 --trials 1 --steps 3000 --resume research/joint-query-stage-search/discovery-results/six.json --output /tmp/six-long.json
python -B research/joint-query-stage-search/search.py --n 12 --trials 2 --steps 400 --random-frame --bound 8 --output /tmp/twelve.json
python -B research/joint-query-stage-search/search.py --n 6 --trials 2 --steps 500 --random-frame --exact-dual --complex --bound 8 --output /tmp/six-complex.json
python -B research/joint-query-stage-search/search.py --n 12 --trials 2 --steps 600 --random-frame --exact-dual --complex --period 6 --bound 8 --output /tmp/twelve-period6.json
python -B research/joint-query-stage-search/search.py --n 6 --trials 3 --steps 1000 --random-frame --exact-dual --complex --period 7 --degrees 3,4,4,4,4,4 --bound 8 --calibrate-excluded --output /tmp/six-period7.json
python -B research/joint-query-stage-search/search.py --n 8 --trials 2 --steps 800 --random-frame --exact-dual --complex --period 7 --degrees 3,3,3,3,3,4,4,4 --bound 8 --calibrate-excluded --output /tmp/eight-period7.json
```

No formal replay is registered for these numerical searches. The adjacent
proof package has separate standard-library exact checks. The optimized
Python rank-gate control also rejected the eight-field period-seven target.
Neither a numerical endpoint nor the complete source/consumer bill passed.

## Open next questions

Can the twelve-field, fourteen-full-call profile attain an exact endpoint?
It meets the forward and absolute stage-rank bounds with equality and has
a favorable local price. First reconstruct an exact word; then charge all
native maps, helpers, ordinary work, restoration and outer obligations.
More seeds in an unstable border regime do not establish this construction.

Can a changed native address map, which does not commute with the stage,
escape the common-stage model? Specify the actual map and test its full
endpoint before assigning any free rank or cost. This archive supplies no
such map. The admitted BIT saving .0034 is another limit on an eventual
integer kappa improvement.

## Source and credit

Base: main `3f0152a32eeebac0df11946dbd50cabe5a7f7839`, priority revision 8.
The selected kappa is .003125804. Canonical Discussion5593, comment19199752.
Retain the P104, native query, stage, helper and compiler lineage, including
Cinna Davis and the credited predecessors. The border belongs to Dillon
and Yiyang. New searches by Ali Romman with substantial OpenAI Codex
assistance. No subagents or independent scientific review. Apache-2.0.
