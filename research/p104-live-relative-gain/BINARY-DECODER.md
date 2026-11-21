# Integer movement cost of three frozen BIT returns

For each choice of three of the nine fixed decoder forms in the #7786 ordinary BIT word, the minimum sum of movement ranks needed to bring the three forms to the empty chart is **21** in the model below. The original three paths have rank seven each. The lower bound allows any linear binary network and any new rational relay charts in the same eight-dimensional space.

## Scope

Use the compressed `research/p104-iterated-returns/binary/WORD.json.gz` at commit `89a67c404226d51028e135de2b9169e9020d57af`. The verifier requires its compressed SHA-256 to be `d28f0d28d9d87d5fe57b8ef19e12faf98014f0bbef781987ce22daaa82f21cb9`. This proof uses the 111 independent live inputs `X_i`. The verifier runs the literal word and its transpose on all 680 public variables. It identifies the nine transfers that arrive at the empty chart. Their transfer forms, after removal of public inputs already available there, depend only on the `X_i`.

Use the input registers and paid receipts of the original 636 charts as local stock. Remove the 27 paid receipts whose backward row equals that of one of the nine final returns. Also omit every paid receipt affected by a live `Y_i`. At each chart grant every linear combination of the remaining stocks after projection onto the `X_i`. The projection relaxes the construction by discarding coefficients on other independent inputs. For the three requested forms, the verifier also checks the original nested paths on this projection. Each live input line that occurs in a stock lies in that stock's rational chart.

A fresh relay starts empty. Each transfer sends one binary linear combination of currently available data. Local linear operations and recursive branching are free. A movement between nested charts costs their dimension difference. More generally, the lower-bound proof applies if the cost of a transfer between `A` and `B` is at least `dim(A) + dim(B) − 2 dim(A ∩ B)`. The original paths use nested charts. The result concerns the integer sum of movement ranks. The tariff `r^(1−κ)` also depends on how these ranks are grouped. The scope fixes the three requested forms and the available stock. Joint work with other returns, changed stock, and larger blocks require separate checks.

## Lattice bound

Let `ℓ_i` denote the rational line of `X_i`. Take the primitive integral representative of each line and scale its Euclidean length to two. All 111 resulting root vectors have integral coordinates. Define `Λ` to be their integral span together with `2 Z^8`. For any rational chart `C`, let `C̄` be the image of `Λ ∩ C` inside `Λ/2Λ`. The sublattice `Λ ∩ C` is saturated in `Λ`, so `dim C̄ = dim C`. Inclusion of rational charts preserves inclusion of their images. In particular the image of `A ∩ B` lies in both `Ā` and `B̄`. Each `ℓ_i` reduces to a binary line with generator `r_i`.

Transpose a proposed network. At the empty chart, start three independent bits. At each chart `C`, let `V_C` be the span of the three-bit coefficient vectors ever available there. Consider the binary space

```text
S = sum over C of C̄ ⊗ V_C.
```

It starts with dimension zero. A local operation preserves `S`. Suppose `A` sends coefficient vector `u` to `B`. The previous `S` already contains `Ā ⊗ u`. The increase is at most `dim(B̄) − dim(Ā ∩ B̄)`, which is at most `dim(B) − dim(A ∩ B)`. The total paid rank bounds `dim S`.

Write `s_i` for the three output coefficients required on the original input `X_i`. In the transposed computation, `s_i` is the sum of the coefficient vectors of local stocks that contain `X_i`, counted with their coefficients. Their charts all contain `ℓ_i`. Consequently `S` contains `r_i ⊗ s_i` for each `i`. Let `D` be the span of these 111 tensors.

The verifier computes `dim D = 21` for one triple and `20` for the other 83. In a rank-20 network for one of those 83, equality forces `S = D`. Therefore a nonzero coefficient vector `u` can reach a chart `C` only if `C̄ ⊗ u` is contained in `D`.

Test that condition for all seven nonzero `u` and each stocked chart. For each of the 83 triples, every eligible stocked chart above dimension two has dimension four. Grant all such charts and every root in their stock for free. Between 58 and 103 required input positions still lie outside their combined support. Each stock at a chart of dimension one or two contains at most two distinct live input positions, even after granting all its local rows. A small chart must receive at least one positive-rank transfer before it can have a nonzero coefficient vector. Twenty units can thus reach small banks that cover at most 40 of the remaining positions. This contradicts the computed minimum of 58. The proof permits other empty relays because they cannot add a live input to a stock.

## Consequence for the binary price

At the trial κ = 0.00391884648, the binary caller needs a saving greater than 0.41631562 per native call. The original three paths cost at most 21 under the tariff `r^(1−κ)`. Each new paid rank is an integer from one through eight. Moving a unit from a smaller positive rank to a larger rank decreases the total tariff until the larger rank reaches eight. Thus a total rank of at least 21 costs at least `2·8^(1−κ) + 5^(1−κ)`. Even with the original charge rounded up to 21, the available saving is less than 0.162. Replacing any one triple in the stated stock falls short of the trial. The proof allows branching and recombination of the three forms at intermediate charts.

## Replay

Python 3.9 or later suffices. The script requires only the standard library.

```sh
python3 -I -B verify_decoder_triples.py
```

The output gives 83 triples of tensor rank 20 and one of tensor rank 21. For the lattice computation, the verifier begins with the basis of `2 Z^8`. When a scaled root has half-integral coordinates, its parity `p` gives the missing vector `sum(p_i b_i)/2`. Replacing a basis element with this vector enlarges the lattice by index two. The change of coordinates is integral on doubled coordinates. For each rational chart, the verifier takes primitive integral representatives in this basis. If a dependence exists modulo two, it replaces one row in the dependence by the primitive representative of the sum. The rows retain their rational span, and their index in the saturated lattice decreases. Once the residues have full chart rank, they span the reduction of `Λ ∩ C`.

The original binary word and its return paths are due to the contributors credited in [PR #7786](https://github.com/CrocSwap/integer-mult-bounds/pull/7786). The later mixed caller consumes this same native supplier. [PR #7830](https://github.com/CrocSwap/integer-mult-bounds/pull/7830)

## All three-dimensional return spaces

The same tensor bound proves a larger claim with 256 choices of lattice. Take any three independent forms in the span of the nine returns. Every linear network with root-line incidence that computes these three forms has total movement rank at least 21. This covers all 788,035 three-dimensional subspaces. Only the 84 original triples have an explicit cost-21 path construction here.

Let `L` be the lattice above. For each nonzero binary character `α` of `L/2L`, also use the index-two sublattice `L_α = {z in L : α · z = 0 modulo 2}`. A rational line has a primitive integral generator in each lattice. Choose an index `k` with `α_k = 1`. A basis for `L_α` consists of `2e_k` and `e_i − α_i e_k` for `i ≠ k`, in the chosen basis of `L`. For a root with coordinates `z`, doubled coordinates in this new basis are `sum_i α_i z_i` at index `k` and `2z_i` elsewhere. Divide by the gcd to get a primitive representative. Reducing modulo two gives the root label for the tensor test. The lattice argument above applies separately to each of these 256 lattices.

The finite replay finds the affine coordinate function of every return on the 111 root labels in the base lattice. For one three-space, it tries the base lattice first. It then takes each return combination in the space whose affine function has constant term one and nonzero linear part `α`. These characters only select lattices for the test. The independently computed tensor rank in a single tested lattice certifies the lower bound. Across all subspaces the first successful attempt has index one through five with respective counts 107,950, 605,441, 71,850, 2,711 and 83.

Each three-space occurs once in the enumeration. Choose three nonzero nine-bit rows with strictly increasing highest set bits. Set both earlier pivot columns to zero in later rows. These are the unique reduced-echelon basis rows of the subspace. The verifier expands all 111 binary tensors and reduces their 24-bit rows over `GF(2)`.

Run the complete check with Python and a C++20 compiler.

```sh
python3 -I -B verify_decoder_spaces.py
```

The output gives the number of subspaces, the certified lower bound and the five first-success counts. This stronger argument needs the input-line incidence of the stock at each chart. It does not need the finite bound on the number of small banks used in the first proof. A different decoder span requires a new calculation.
