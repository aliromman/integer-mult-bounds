# Validation record

The scientific package at `b90438ec9d51e7babc0b6b6cb9998a9df81c123d` passes
the following author checks with Python 3.12 and the standard library.

```sh
python3.12 -B scripts/research_package.py --package research/p104-paid-interior-bases
python3.12 -B scripts/check_research_commit.py --package research/p104-paid-interior-bases
python3.12 -B scripts/export_research_package.py export --package research/p104-paid-interior-bases --ref HEAD --output /tmp/p104-interior-science.tar
python3.12 -B scripts/export_research_package.py check /tmp/p104-interior-science.tar
```

The registered hook runs both `-I -B` and `-I -O -B`. Both runs give
`PASS_P104_PAID_INTERIOR_BASES_KAPPA_00310815849`, followed by
`PASS_REGISTERED_RESEARCH_PACKAGE`. They check the source exports, the four
construction steps, all 565,504 transpose entries, 2,918 signed native maps,
the complete physical caller, both rational price engines, all compiler
terms, all 58 outer conditions, and the invalid-word controls in RESULT.json.

All 16 package files match the scientific commit. The official export passes
its integrity check. Its 17 files include the export manifest. An isolated
optimized replay of this export passes in 89.70 seconds with an empty PATH,
no Git checkout, no network request, and read-only source files. Every file
hash and mode stays fixed.

Separate temporary copies test archive-byte corruption, payload corruption
with a new outer digest, and a missing payload. Each case is tested on both
the outer export and the nested predecessor export, with and without `-O`.
All 12 cases reject before any retained research module is imported. These
are source-integrity tests; they do not replace the scientific replay.

The complete conditional claim is `0.00310815849`. The tested saving is
`0.003108158495`, and the exact raw crossing is
`(0.00310815850980, 0.00310815850981)`. The source and final bills are checked
separately. Their difference is `18*(-7 H1 + 2 H2 + H3)` in paid requests,
with no change in productive targets.

The branch starts from canonical main
`df0137efd35e11b4445734bf612bdd4c7ebc330c`. The source archive is retained in
its direct child `89c7e069d6b93b7802f946d01886d92970303f86`. Full-history
ancestry mode was not run because this checkout is shallow.

Hosted CI and independent scientific review are separate statuses. The five
inherited assumptions and helper endpoint contract remain. The huge uniform
word and complete rounded multiplier are not materialized. The selected
frontier is unchanged. The numerical gain is small; no large improvement or
global optimum is claimed.
