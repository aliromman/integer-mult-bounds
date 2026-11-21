# A complete conditional application with 630 binary fields

The two changed circuits give **κ ≥ 0.00388743431214** under the five assumptions in the source framework. PR 7863 gives 0.00387992804346. The relative gain is 0.19346412%. The selected record on main has its own review. [Finite receipt and comparison](APPLICATION.json) · [PR 7863](https://github.com/CrocSwap/integer-mult-bounds/pull/7863) · [Selected record](../../certificates/selected-result.json)

| Checked quantity | Binary | Rational |
| --- | ---: | ---: |
| Native calls in the full invoice | 36 | 48 |
| Independent public fields per call | 630 | 1,291 |
| Paid transfers per call | 2,803 | 5,676 |
| Lower endpoint of the raw crossing | 0.003887434522233 | 0.00389749436514 |
| Upper endpoint of the raw crossing | 0.003887434522234 | 0.00389749436515 |

## Binary construction

The binary circuit uses the same 111 physical main port quartets as the iterated-return supplier. Each input pair has the scalar map (x, y) ↦ (x, y+x). The other 408 inputs have separate identity outputs. The construction first routes the 111 data inputs to their receivers. At each chart, independent input columns can fill unused directions of the outgoing coefficient matrix. The rest of the incoming independent space drains to the terminal bank. Choose a basis for these spaces, then carry each new input and its eventual return in the word. This completes the data routes to a map on arbitrary independent public inputs. [Literal word](BINARY-CANDIDATE.json.gz)

The choice of basis changes where an individual helper can return. On eleven terminal changes, a combination of deliveries from one chart affects an unprotected helper direction. Rotate the remaining deliveries to remove their component in that direction. A linear form at an earlier bank distinguishes the affected helper. Return the form there and rotate its independent source coordinates. If the form needs an additional delivery, charge the delivery between the actual charts. The resulting public returns consist of thirteen grade-four and 395 grade-eight helpers. The complete binary word contains 29,427 XOR additions, 3,433 initializations, 3,433 clears, 2,803 paid moves and 630 swaps. Its forward and reverse maps each check all 396,900 scalar coefficients. The omitted-clear control must fail.

The paired caller retains 36 native calls. Its literal readback accounts for the root targets, 32 idles and four queries across 111 ports. The 768-factor finite compiler starts at the admitted ordinary seed of 0.00385248. It charges discarded targets, checks the uniform and delayed-cover losses, then runs four finite bootstrap stages. It admits 0.003887434312223 at the ordinary BIT interface. [Finite counts, histograms and cap](APPLICATION.json) · [Predecessor finite compiler](https://github.com/CrocSwap/integer-mult-bounds/blob/89a67c404226d51028e135de2b9169e9020d57af/research/p104-iterated-returns/binary/finite.py)

## Rational consumer

The rational circuit uses the four-term quotient and two independent helper splits. Its 48-call invoice includes all caller work. The smallest common product grid is 154113755636630124932500335624192000. The basis replay includes 25,038,131 products. The bounded scalar compiler lowers 29,692 occurrences with literals wider than 24 bits. Its independent reader executes the lowered programs in four signed components and checks dirty register retirement, ordinary work and the depth bound. [Rational construction](BIDIRECTIONAL.md) · [Grid, fixture and ordinary work](APPLICATION.json)

The compiler passes the measured ordinary BIT cap to the rational consumer. It checks 54 strict and four functional inequalities. The smallest outer margin is 10⁻¹⁵. The five inherited hypotheses concern the native operators on all columns, owned erasure, the all-size compiler, tape precision and analytic recovery. The finite replays check the particular circuits and their interfaces. Independent scientific review of this changed application remains open. [Inherited source and scope](https://github.com/CrocSwap/integer-mult-bounds/tree/89a67c404226d51028e135de2b9169e9020d57af/research/p104-iterated-returns/continuation)

## Reproduce

Run from a checkout with CPython 3.12 and the standard library.

```sh
python3.12 -I -B research/p104-live-relative-gain/verify_application.py
```

The program extracts 134 source files into a temporary directory. Seven files come from the pinned PR 7786 package. Its original source lock checks the SHA-256 and Git blob of each of the other 127 files. The binary compiler rereads the proposed literal word. The rational compiler rebuilds the proposed program from its journals. The independent reader then rebuilds the grids, signed costs, ordinary costs and outer conditions. Use `--out FILE` to retain the full physical receipt. [Replay](verify_application.py) · [Offline sources with retained notices](APPLICATION-SOURCES.tar.gz) · [Source lock](https://github.com/CrocSwap/integer-mult-bounds/blob/89a67c404226d51028e135de2b9169e9020d57af/research/p104-iterated-returns/SOURCE-LOCK.json)

## Open next questions

The comparison with PR 7863 and the previous round's gain requires κ greater than 0.0039188464750082730562. Both native tariffs have deficits at 0.00391884648. A cheaper paid delivery or a new independent split must preserve each primary and auxiliary output. Can a change to the internal data routes make another low-grade helper available before a paid terminal delivery?

Ali Romman prepared the new circuits, composition and replay with substantial OpenAI Codex assistance. Clear1nsight supplied the pinned iterated-return programs and independent scalar audit. The project also credits Jerry Chen, DreamingOfClouds, zxv, SYLforge, Rohan Arun, Cinna Davis, Douglas Colkitt and the predecessors named in the retained notices. [Attribution](PARENT-NOTICE)
