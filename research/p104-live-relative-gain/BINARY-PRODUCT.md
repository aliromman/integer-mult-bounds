# A mixed output in the two-contact binary caller

The following short caller has two proposed rank-15 main outputs at each of 111 ports. One term from the second input also appears in the first output. The natural two-level chronological product sends that term between two roles of equal Hamming weight. A restriction to one weight therefore fails for these matchers.

## The single call

Let Q act on the shared head and let D₁ and D₂ each act on seven new body axes. The operators commute and Q² = I. At one port a native contact on the ordered pair (a,b) returns (Dₛa, Dₛb + QDₛa). The sum is over the binary field. Apply these instructions in order.

1. Swap the two fields.
2. Add the second field to the first.
3. Run contact 1 on the first and second fields.
4. Add the first field to the second.
5. Run contact 2 on the second and first fields.
6. Swap the fields.

For inputs (x,y) the result is

$$
D_1D_2\big(Qx+(I+Q)y,\;Qy\big).
$$

For each of the 222 main fields the diagonal job uses 15 axes. The unchanged helper jobs belong to the two separate native calls. The extra term is on the same body as the desired job. The full word must remove it to claim both target credits.

## Two-level check

Take the chronological product with the same two-contact caller at each paid position. Use the natural inverse coordinate matchers. Number both outer and inner roles 0 and 1 in the order (x,y). Put one input bit at role (0,1) of port 0. Inspect output (1,0) at the same port. Both role labels have weight one. The checker finds **50 distinct nonzero output addresses** in this component. Replacing each inner contact by its diagonal part gives zero in the same component. [Literal replay](verify_binary_product.py)

The address calculation is exact over the rationals. For the two exchanged address chunks, use their difference as the coordinate and fix their sum to zero. A swap on a native chart U acts by the reflection Rᵤ, with sign −1 on U and +1 on its complement. Only the shared inner head slice carries a nonzero difference. All other inner slices have zero difference. A paid outer move from chart A to chart B acts on this slice by RᵦRₐ. With the data in chart coordinates, the direct signal passes through. The mixed component gains the direct signal and its image under RᵦRₐ. The charts of a paid edge are nested, so these two reflections commute. Each binary ADD uses symmetric difference of the two finite address sets.

The proof uses both complete native executions, including the actual order of all 2,874 paid moves in each. At the boundary between the two contacts, a mixed value at a different port or helper has zero transfer back to port 0 through the next diagonal native endpoint. This permits the checker to retain only this port at that boundary. The input, one nonzero mixed output and the diagonal control are fixed in the script. Reduction modulo any sufficiently large odd prime preserves all distinct addresses in this finite computation.

Run the check with CPython 3.12 and its standard library.

```sh
python3.12 -I -B research/p104-live-relative-gain/verify_binary_product.py
```

The binary word comes from Clear1nsight's iterated-return supplier. The source digest is checked before each replay. Ali Romman prepared this test with substantial OpenAI Codex assistance. [Pinned predecessor word](https://github.com/CrocSwap/integer-mult-bounds/blob/89a67c404226d51028e135de2b9169e9020d57af/research/p104-iterated-returns/binary/WORD.json.gz) · [Chronological product contract](https://github.com/CrocSwap/integer-mult-bounds/blob/e54523517d0f03e58352e647c403a5af3e525a57/research/p104-shared-native-helper-splits/retained/compiler/PROOF.md)
