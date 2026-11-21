# Ten changes to the rational P104 circuit

Seven changes of bank basis, one new helper and three removed transfers raise the rational raw crossing from the interval [0.00388013644474, 0.00388013644475] to [0.00388332156403, 0.00388332156404]. These numbers price all 48 calls and their fixed caller. The 36-call binary circuit limits the parent conditional multiplication result to 0.00387992803556. PR 7863 uses the same native profiles and a 451-factor schedule to obtain a conditional 0.00387992804346. [Parent result and both raw intervals](https://github.com/CrocSwap/integer-mult-bounds/blob/87fbf4f7b35b0f18ab12b5767775a9a424cf0b15/research/p104-mixed-rank-supplier/RESULT.json) · [451-factor proposal](https://github.com/CrocSwap/integer-mult-bounds/pull/7863)

The first three edits change the forward graph. The next four change its transpose and then transpose it back. The eighth edit sends transfer 96070 through chart 1048. An independent input there splits the return for helper 449. The ninth edit removes transfer 92453 between charts 480 and 498. Its five receivers share the same input chart. Four independent differences cancel the transfer. The fifth helper, role 1255, uses the earlier receipt 91155 and returns at chart 480. The tenth edit gives up the old helper 324. One relation on its 404 input coefficients removes both transfers 92917 and 94384. One paid transfer has grade three and the other has grade one. This costs one grade-four helper. The new independent helper takes the vacant role number 324. Each column still has its required scalar coefficients at its stated output chart. The ordered journal records the original and replacement charts, the paid transfers, each basis change, and the coefficients of the new return. [Edit journal](EVENTS.json)

## Check the result

Run this command from the repository root with CPython 3.12. It uses the standard library.

```sh
python3.12 -B research/p104-live-relative-gain/verify_component.py
```

The check reads the supplied graph and verifies all 1,289 independent input columns by reverse accumulation over the rationals. It visits 28,670,764 adjoint terms. It checks every source chart and each nested paid or productive transition. It rebuilds the two complete histograms and checks the signs at both reported endpoints. Successive integer square roots give outward bounds for the powers. [Independent checker](verify_component.py) · [Histograms and result](RATIONAL-RECEIPT.json)

Clear1nsight independently replayed all 1,289 rows in both directions and checked the first reader of each of the 5,700 transfers. The 48 calls have 282,816 paid and 52,512 productive occurrences. An ordered check of the helper owners rejected a swap that preserved the sorted digests. [Independent audit in PR 7906](https://github.com/CrocSwap/integer-mult-bounds/pull/7906)

The repository hook runs the same check with its source pins. Run `python3.12 -B scripts/research_package.py --package research/p104-live-relative-gain`. [Source manifest](SOURCES.json) · [Registered check](VERIFY.json)

The forward replay uses the pinned parent word, layout and bank operations. To run it, use the three source files or directories listed below. Supply local copies at the indicated arguments.

```sh
python3.12 -B research/p104-live-relative-gain/replay.py \
    --source "$PARENT_GRAPH" --layout "$PARENT_LAYOUT" --engines "$PARENT_ENGINES"
```

The replay reconstructs the graph from the journal and checks its canonical digest. It checks 13,109,869 scalar products on an integer grid. Of those products, 7,556,199 use a denominator greater than one. The minimum common product grid is 35237207388978003467690311680000. The largest scalar literal uses 37 bits. [Forward receipt](AUDIT.json)

| Input | Source |
| --- | --- |
| Parent graph | [Clear1nsight, PR 7830, DOUBLE-CANDIDATE.json.gz](https://github.com/CrocSwap/integer-mult-bounds/blob/87fbf4f7b35b0f18ab12b5767775a9a424cf0b15/research/p104-chart306-chain/DOUBLE-CANDIDATE.json.gz) |
| Parent layout | [Clear1nsight, PR 7786, continuation/LAYOUT.json.gz](https://github.com/CrocSwap/integer-mult-bounds/blob/89a67c404226d51028e135de2b9169e9020d57af/research/p104-iterated-returns/continuation/LAYOUT.json.gz) |
| Parent engine directory | [Ali Romman, PR 7776, complex](https://github.com/CrocSwap/integer-mult-bounds/tree/4f7c67ae387cc91580f469aeda4d508b372b9ec1/research/p104-dual-supplier-lift/complex) |
| Fixed caller counts | [PR 7830, caller_census in DOUBLE-APPLICATION.json](https://github.com/CrocSwap/integer-mult-bounds/blob/87fbf4f7b35b0f18ab12b5767775a9a424cf0b15/research/p104-chart306-chain/DOUBLE-APPLICATION.json) |

## Binary caller test

A two-contact binary word leaves a mixed term. A literal replay of its two-level chronological product gives 50 wrong addresses on an input and output of the same Hamming weight. The calculation uses the actual paid moves and exact rational address coordinates. [Counterexample and reproduction](BINARY-PRODUCT.md)

Every triple of the nine frozen binary decoder returns needs total movement rank 21 on the retained stock. A separate lattice check gives the same lower bound for all 788,035 three-dimensional subspaces of their span, including changes of basis. A new network for one triple can save less than 0.162 per call at the trial. The binary circuit needs more than 0.416 per call. [Proof and exact replay](BINARY-DECODER.md)

## Open next questions

At the trial 0.00391884648 the rational signed margin is in [−36.1508758544792045, −36.1508758544792044]. The binary margin is in [−14.9873626241398933, −14.9873626241398932]. Both must become positive before qualification at this trial. The five framework assumptions apply to the parent complete result. A change that clears both tariffs also needs a complete compiler, precision, ordinary BIT and outer application replay. [Coordination](https://github.com/CrocSwap/integer-mult-bounds/discussions/7111#discussioncomment-19777325) · [Prior contribution and comparison](https://github.com/CrocSwap/integer-mult-bounds/pull/7776)

Ali Romman prepared this local research with substantial OpenAI Codex assistance. The circuit uses Clear1nsight's paid-chain contractions and the PR 7776 bank and helper operations. The retained construction credits DreamingOfClouds, CubeJerry, Rohan Arun, SYLforge, Cinna Davis and the other authors in the source notices. [Predecessor credits](https://github.com/CrocSwap/integer-mult-bounds/blob/87fbf4f7b35b0f18ab12b5767775a9a424cf0b15/research/p104-mixed-rank-supplier/NOTICE)
