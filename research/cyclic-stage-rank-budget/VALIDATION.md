# Validation record

Date: 2026-10-29. Host: macOS arm64. Runtime: CPython 3.12.13, standard
library only. Main was refreshed to
`3f0152a32eeebac0df11946dbd50cabe5a7f7839`. Priority revision 8 and the
selected record were unchanged. The live overlap read covered the 30 most
recent Discussions, their last five comments, the last 30 handoff comments
and the 50 most recently updated open PRs. It found no new overlapping
claim. This bounded read is not a proof of novelty.

The registered package replay passed all four commands in VERIFY.json:

```sh
python3.12 -B scripts/research_package.py --package research/cyclic-stage-rank-budget
```

The receipt ended with
`PASS_REGISTERED_RESEARCH_PACKAGE research/cyclic-stage-rank-budget outcome=negative-result`.
The universal results rely on the written proofs. The finite checks support
the formulas and failure controls; they are not formal proof acceptance.

A separate temporary export contained this package and its two input files.
It had no Git directory. Both programs ran from an unrelated working
directory, normally and with `-O`, with input files set read-only. The file
list and every file hash stayed unchanged. The normal and optimized stdout
hashes agreed:

| Program | SHA-256 of stdout |
| --- | --- |
| verify.py | ca56a3018c08a9ffa2d0d11dbf88d6ea925c5fa050e8c16b471c0bd7a1e5a93b |
| state_check.py | 524028978042c2bda6bf4976f47278632203161060b32d12605e90d05f74cc18 |

Two mutations were tested only in the temporary export. Appending one byte
to the source certificate was rejected with `A scientific input hash differs`.
Removing every nonzero cyclic alias from the state check was rejected with
`Cyclic coefficient removal differs`. Both checks ran with `-O`.

The input hashes and explicit source heads are in SOURCES.json. The full
historical provenance preflight was not run in this shallow checkout. No
selected native caller, compiler, BIT supplier, precision or outer replay
was run for this proof package. Hosted replay of this new package and
independent scientific review were pending at preparation. No new kappa
is claimed. The selected value remains 0.003125804.
