# PolarAED — AED of polar codes with UTL / block-LTA automorphisms

Extension plan. Target: reproduce the UTL and BLTA automorphism-ensemble (AE)
decoding experiments for *plain* polar codes inside this repo, reusing the POD
decoder, the BSGS/Schreier–Sims tooling and the MClassifier stack, and then run
MClassifier on this new code family.

Reference papers (the only two needed):

* **[P21]** C. Pillet, V. Bioglio, I. Land, *Polar Codes for Automorphism Ensemble
  Decoding*, ITW 2021 — arXiv:2102.08250. Defines UTL, the admissibility
  criterion (Thm. 1/2, Cor. 1), and the UTL-friendly code design (§III-C).
  Already in `POD_itw2026/ref.bib` as `pillet2021AEDpolar`.
* **[P22]** C. Pillet, V. Bioglio, I. Land, *Classification of Automorphisms for
  the Decoding of Polar Codes*, arXiv:2110.14438 (+ journal version
  *Group Properties of Polar Codes for AE Decoding*, arXiv:2206.03342, T-IT 2023).
  Defines BLTA(S), the SC-absorption group [1], equivalence classes (EC), the
  PUL decomposition, EC-representative sampling and the automorphism-friendly
  design Algorithm 1.

---

## Status (2026-09-20)

| phase | state |
|---|---|
| P1 `polar_design.py`, `monomial_aut.py`, gates | **done** |
| P2 `blta_sample.py`, `verify_aut.py`, `n128_k100_highsnr` = [P21] Fig. 3 | **done** — `codes/n128_k100_highsnr/{README,RESULTS}.md`, `fig3_repro.png`. LTA byte-identical to SC; UTL ~0.9 dB better; AE4-SCL8 == SCL32 |
| P3 `design_utl.py` → (1024,512) `U1`/`U2`, [P21] Fig. 4 | not started |
| P4 `design_blta.py` → [P22] Figs. 3–5 | not started |
| P5 `blta_bsgs.py` cross-check | not started |
| P6 MClassifier on the AE code | not started |

---

## 0. Executive summary — what actually has to be built

**The C++ decoder needs no changes.** This was verified, not assumed (§1.4).
`AdjustPolarDecoder` (`src/ErrorCorrectionCode/AED_attemp_20260102.cpp`) already:

* derives the frozen / dynamic-frozen `relation_ship` from `polar_matrix · P · H`;
  feeding it a **plain polar code's** H and an **identity** `permutation_src`
  makes every constraint a *static* frozen bit — which is exactly the
  UTL/BLTA setting;
* loads a list of permutations from `automorphism_src` and runs one SC/SCL
  branch per permutation, selecting by best path metric (`AED_relation_check.cpp`).

So the work is **entirely on the Python side**: build the polar code, compute its
affine automorphism group, sample non-redundant automorphisms, and emit the two
files the decoder already knows how to read (`*_H.matrix`, `aut_*.txt`).

MClassifier is `-ini`-driven on top of the same decoder class, so it transfers
with new `.ini` files plus retraining.

---

## 1. What already exists, and the key reuse facts

### 1.1 Decoder (no change needed)

`AdjustPolarDecoder` constructor, in order:

1. builds `polar_matrix` from `operationArray` (all-`1` ⇒ full `F^{⊗n}`);
2. loads `permutation_src` → `received_order` (empty ⇒ identity);
3. `constraint = polar_matrix · P · Hmatrix`, transpose, Gauss–Jordan →
   `relation_ship[i]`: empty ⇒ information bit, `{i}` ⇒ frozen,
   longer ⇒ dynamic frozen;
4. `load_automorphism_set()` composes `received_order_set[ℓ][j] = received_order[h_ℓ[j]]`.

`AdjustPolarDecoderRelation::doDecode` runs all branches and keeps the minimum
path metric. [P21]/[P22] instead select the branch minimising ‖y − x̂‖². Under
min-sum SC the accumulated path metric is a monotone *approximation* of that,
not the same quantity — the internal LLRs it sums are themselves approximate.
Good enough to reproduce the qualitative results (it did, see
`codes/n128_k100_highsnr/RESULTS.md`), but a true LS combiner is a ~20-line
opt-in addition to `AED_relation_check.cpp` and should be in place before
quoting absolute gaps against the papers.

### 1.2 Index convention (verified, §1.4)

`polar_matrix[r][c] = 1 ⟺ c ⊆ r` (bitwise submask) — i.e. `F^{⊗n}` in natural
order, **no bit reversal**. This coincides exactly with the monomial-code
convention of [P21, Table I]:

* coordinate index `j` ↔ variable vector `(x_0, …, x_{n−1})`, `x_i = bit i of j`
  (**LSB-first**);
* row index `r` ↔ negative monomial `m_r` containing `x̄_i` **iff bit i of r is 0**;
  `r = N−1` ↔ the constant monomial (weight-N row), `r = 0` ↔ degree-n monomial.

Affine transformation → coordinate permutation:

```
π(j) = idx( A · bits(j) + b ),    bits LSB-first
```

This direction (not `Aᵀ`, not `A⁻¹`) is the one under which **lower-triangular A
gives automorphisms of every decreasing monomial code** — confirmed numerically
in §1.4.

### 1.3 Python tooling to reuse

* `project/POD/schreier_sims.py` — Schreier–Sims → BSGS JSON with transversals;
  also `sift_remainder_left` / `right_coset_key` for coset canonicalisation.
* `project/POD/convert_bsgs_to_perm.py` — uniform sampling from a BSGS chain,
  writes the `"<rows> <n>"`-header permutation file the decoder reads.
* `project/POD/generator_set_eBCH.py` — the template for a
  "build code + build automorphism generator set" script (`nullspace_mod2`,
  `rank_mod2`, `is_automorphism`, `save_matrix_txt`).
* `project/POD/*/plot_logs.py` — log plotting; reuse as-is.

### 1.4 Feasibility checks already run (evidence, not plan)

Run against the prebuilt `project/POD/m7t10_itw2026/POD` binary, unmodified.

1. **Plain polar code ⇒ static frozen set.** Built the [P21] §III-A example
   (16, 7), `I = {7,10,11,12,13,14,15}`, `H = nullspace(G)` saved transposed.
   With `permutation_src` empty, the constructor printed the relation string
   `1111111011000000` — no `2`, i.e. **zero dynamic-frozen rows**, and the
   information positions are exactly `I`. The existing decoder therefore
   handles plain polar codes with no code change.

2. **Admissibility test reproduces the paper.** For that same code the test
   "`A[i][j]=1` admissible ⟺ for every `r ∈ I` with bit `i` of `r` equal to 0,
   `(r | 1<<i) & ~(1<<j) ∈ I`" yields upper-triangular admissible set
   `{(1,2)}` and *all six* lower-triangular positions — matching [P21]'s worked
   example (`A_{1,2}` free; `A_{0,1}`, `A_{0,2}`, `A_{·,3}` not) and confirming
   LTA ⊆ Aut.

3. **Permutation direction settled.** For the (16,10) code
   `I = {5,6,7,9,10,11,12,13,14,15}` (UT-admissible = `{(0,1),(2,3)}`),
   `perm_from(A)` is an automorphism for the admissible lower entries while
   `perm_from(Aᵀ)` is not. Convention confirmed as §1.2.

4. **SC absorption of LTA, end-to-end through the C++ binary.** Same (16,10)
   code, AWGN @ 1.0 dB, SC (L = 1), identical seeds, single-branch AE:

   | run | BLER |
   |---|---|
   | plain SC | `100000/463428 = 0.215783` |
   | AE with one **LTA** element (`A[3,0]=A[2,1]=1`) | `100000/463428 = 0.215783` — **bit-identical at every checkpoint** |
   | AE with one **UTL** element (`A[0,1]=A[2,3]=1`) | `100000/463276 = 0.215854` — different |

   LTA is absorbed, UTL is not. The whole chain (index convention → permutation
   file → decoder) is validated.

All four are reproducible in one shot:

```
cd project/PolarAED/tools
python3 _feasibility_check.py --workdir /tmp/polaraed_feas \
    --binary ../../POD/m7t10_itw2026/POD        # --no-sim for algebra only
```

`tools/_feasibility_check.py` is the seed of the real tooling: it already
contains the conventions (`polar_matrix`, `admissible`, `perm_from_affine`,
matrix/permutation file writers, `.ini` template) that §4 factors out into
proper modules. Current output: `check 1..4 : PASS`.

---

## 2. The math to implement

### 2.1 Admissibility matrix

For information set `I ⊆ [N]`, `N = 2^n`:

```
adm(i, j) = ∀ r ∈ I with bit_i(r) = 0 :  ((r | 1<<i) & ~(1<<j)) ∈ I
```

(the variable change `x_i ← x_j` of [P21, Thm. 2], applied to binary expansions).
`adm(i,j)` for all `i > j` ⟺ the code is decreasing / UPO-compliant
(⟹ LTA ⊆ Aut, [P22, Thm. 2]). Assert this and warn loudly if it fails.

### 2.2 Block structure S

The admissible **upper** positions of a UPO-compliant code form the strict upper
triangles of diagonal blocks ([P22, Thm. 3]): the affine automorphism group is
`A = BLTA(S)`, `S = (s_1, …, s_t)`, `Σ s_i = n`, blocks laid along the diagonal
starting at variable index 0 (LSB). Recover `S` by scanning: block boundaries are
the indices where no admissible upper entry crosses. Cross-check with
`|BLTA(S)| = 2^{n(n+1)/2} · Π_i Π_{j=2}^{s_i} (2^j − 1)` ([P22, Lem. 14]) against
a brute-force automorphism count for `n ≤ 4`.

Sanity anchors: RM ⇒ `S = (n)`; the [P21] (16,7) example ⇒ `S = (1,2,1)`.

### 2.3 SC-absorption group and equivalence classes

* `[1] ⊇ LTA(n)`; if `s_1 > 1` then `[1] ⊇ BLTA(2,1,…,1)` ([P22, Lem. 13]).
* `[1] = BLTA(S₁)` for some refinement `S₁` of `S` ([P22, Thm. 4]). Working
  assumption, as in the papers: **`[1] = BLTA(2,1,…,1)`** when `s_1 > 1`,
  else `[1] = LTA(n)`. Verify empirically (§5, gate G4) — a larger `[1]` shows up
  as branches that produce bit-identical output.
* Number of ECs, `E = |BLTA(S)| / |[1]| = (1/3) Π_i Π_{j=2}^{s_i} (2^j − 1)`
  ([P22, Eq. 48]) — upper bound on a useful `M`.
  Anchors from the papers: `S = (4,1,1,1,3)`, n=10 ⇒ `315 · 21 / 3 = 2205`;
  `S = (3,5)`, n=8 ⇒ `21 · 9765 / 3 = 68355`. Both reproduce exactly.
* **PUL decomposition** ([P22, Thm. 5 / Lem. 17]): every `A = P·U·L` with
  `L ∈ LTA`; since `L` is absorbed, every EC has a representative `P·U` with
  `P ∈ A_P` (block-wise permutation matrices, `|A_P| = Π s_i!`) and
  `U ∈ A_U` (block-wise UTL, `|A_U| = Π 2^{s_i(s_i−1)/2}`).
  **The translation `b` is always absorbed** (`(I,b) ∈ LTA`), so all
  representatives can be taken purely linear, `b = 0`.
* **Exact EC test** ([P22, Lem. 6]): `π₁ ∈ [π₂] ⟺ A₁ · A₂⁻¹ ∈ [1]`.
  This is an `O(n³)` structural check on a ≤10×10 binary matrix — cheap and exact.

---

## 3. Directory layout

```
project/PolarAED/
  PLAN.md                    <- this file
  CMakeLists.txt             <- copy of project/POD/CMakeLists.txt (project name PolarAED)
  main.cpp                   <- copy of project/POD/main.cpp (unchanged)
  tools/
    polar_design.py          # reliability -> I -> G/H .matrix
    monomial_aut.py          # adm matrix, S, |BLTA(S)|, |EC|, A_G auxiliary matrix
    design_utl.py            # [P21] §III-C design
    design_blta.py           # [P22] Algorithm 1 design
    blta_sample.py           # EC-representative sampler -> aut_*.txt   (primary)
    blta_bsgs.py             # BLTA generator set -> schreier_sims.py -> BSGS  (cross-check)
    verify_aut.py            # automorphism + absorption verification
    data/5g_reliability_1024.txt
  codes/
    n128_k100_highsnr/       # [P21] Fig. 3
    n1024_k512_U1/  U2/      # [P21] Fig. 4
    n1024_k512_S41113/       # [P22] Fig. 3
    n256_k128_S35/  S53/     # [P22] Fig. 4 / Fig. 5
      <code>.matrix, <code>_H.matrix, aut_*.txt, *.ini, <run>/log.txt, plot_logs.py
```

Each `codes/<id>/` mirrors a `project/POD/<code>_itw2026/` directory so
`plot_logs.py` and the existing run scripts work verbatim.

---

## 4. Tools to write

### 4.1 `tools/polar_design.py`

* `ga_reliability(n, design_snr_db) -> list[int]` — Gaussian-approximation
  density evolution, returns bit-channel indices in decreasing reliability.
  (Do **not** reuse the C++ `bhattacharyya_value`; the papers design with DE/GA,
  and the C++ value is only used inside the decoder for reporting.)
* `load_5g_sequence(N)` — 3GPP 38.212 Table 5.3.1.2-1, truncated to `N`.
  **Needs a data file; not currently in the repo** (flagged in §8).
* `rm_information_set(n, r)`.
* `info_set_to_matrices(n, I) -> (G, H)` with `G = F^{⊗n}[I, :]`,
  `H = nullspace(G)`; write `G` as `<id>.matrix` and **`Hᵀ`** as `<id>_H.matrix`
  (`N × (N−K)`) — same convention as `generator_set_eBCH.py`.
* CLI: `--n --K --design {ga,5g,rm} --snr --out_dir --id`.

### 4.2 `tools/monomial_aut.py`

* `admissible_matrix(n, I) -> bool[n][n]` (§2.1).
* `is_decreasing(n, I)` — all lower entries admissible.
* `block_structure(adm) -> S` (§2.2).
* `blta_size(S)`, `num_ec(S)`, `absorption_profile(S) -> S₁`.
* `AG_matrix(n, I) -> list[list[set[int]]]` — for each `(i,j)`, the set of bit
  indices that must be **added** to `I` to free position `(i,j)`; this is the
  `A_G` of [P22, Alg. 1] and the `A_M` of [P21, Eq. 4]. Computed by collecting
  `((r|1<<i) & ~(1<<j))` over `r ∈ I, bit_i(r)=0` and keeping those not in `I`
  (closed under the decreasing-monomial requirement).
* CLI prints a report: `S`, `|Aut_aff|`, `|EC|`, `|A_U|`, `|A_P|`, admissible map.

### 4.3 `tools/design_utl.py` — [P21] §III-C

Input `n, K, s, reliability sequence R, target UT positions`. Start from
`I_s = R[0 : K−s]`; repeatedly pick a target admissible entry `A[i][j]`
(prefer the **bottom-right**, `n/2 < i < j`, per [P21]'s conclusion), add the
`p ≤ s` indices listed in `A_G[i][j]`, recompute `A_G`, until `|I| = K`.
Output `I`. Reproduces the paper's `U1` (3 free UT positions ⇒ 8 UTL elements,
upper-left) and `U2` (bottom-right concentrated) for `(1024, 512)`.

### 4.4 `tools/design_blta.py` — [P22] Algorithm 1

Two nested loops: outer over design SNR (`SNR_min : ΔSNR : SNR_max`), inner
decrementing `K_s` from `K−1`; free the blocks of the target `S` column by
column using `A_G`, recomputing `A_G` after each freed column; succeed when
`|G| = K`, restart with `K_s−1` on overshoot. Early stop via the
`Σ_{k≤d_max} C(n,k) < K` check (Line 8). Guarantees the decreasing-monomial
property, unlike §4.3. Targets: `S = (4,1,1,1,3)` for `(1024,512)`,
`S = (3,5)` and `S = (5,3)` for `(256,128)`.

### 4.5 `tools/blta_sample.py` — the ensemble (primary path)

```
sample_ec_representatives(S, M, mode, D=(dU,dP), seed) -> list[perm]
```

* `mode = "PU"` — [P22, §III-C]: draw `P ∈ A_P` (independent permutation of the
  variables inside each block) and `U ∈ A_U` (independent strict-upper-triangular
  fill inside each block), form `A = P·U`, `b = 0`.
* `mode = "UTL"` — [P21]: `U` only, support restricted to admissible UT positions.
* `mode = "random"` — uniform from `BLTA(S)` (random `GL(s_i)` blocks by
  rejection, free entries below the block diagonal, `b = 0`); the papers'
  "random A" baseline.
* **Dedup:** reject a candidate `A` if `A · A_prev⁻¹ ∈ [1]` for any accepted
  `A_prev` (exact, [P22, Lem. 6]). Optionally also enforce the Hamming-distance
  heuristic `D = (d_U, d_P)` on the `v` (UTL fill) and `p` (block permutation)
  descriptor vectors — this is what [P22] shows helps at small `M`.
* Emit `aut_<id>_<mode>_M<M>.txt` with the `"<M> <N>"` header the decoder expects
  (`load_automorphism_set` does **not** prepend the identity, so include it
  explicitly as row 0 if you want the baseline branch in the ensemble).

### 4.6 `tools/blta_bsgs.py` — BSGS cross-check (the requested reuse)

Builds a **generator set of `BLTA(S)` as permutations of `[N]`** — block `GL(s_i)`
generators, the elementary transvections below the block diagonal, and the `n`
basis translations — writes it in `generator_set_*.json` format, then pipes it
through the existing `schreier_sims.py` → `convert_bsgs_to_perm.py`.

> **Concern, stated once and then worked around.** For `N = 128` this is fine
> (chain storage ≈ `Σ|orbit_i| · N` ints). For `N = 1024` with `|BLTA(S)| ≈ 2^55`
> the stored transversals blow up to tens of millions of length-1024 permutations
> — the JSON alone would be multiple GB. And uniform BSGS sampling would still
> need the EC filter afterwards, because `[1]` is itself huge. Since `BLTA(S)` is
> known in *closed form*, §4.5 samples it exactly and in closed form, and the
> EC test of [P22, Lem. 6] is an `O(n³)` matrix check rather than a group sift.
> So: **§4.5 is the production path; §4.6 is the `N ≤ 256` validation
> cross-check** (agreement of sampled group order, membership of §4.5's output in
> the BSGS group) and stays available for the generic "automorphism group only
> known by generators" case — which is exactly the POD/eBCH situation it was
> written for.

### 4.7 `tools/verify_aut.py`

* algebraic: `rank([G ; G[:,π]]) == rank(G)` for every emitted permutation
  (reuse `rank_mod2`/`is_automorphism` from `generator_set_eBCH.py`);
* structural: every pair in distinct ECs;
* behavioural: run the decoder on a few hundred frames and assert that branches
  in distinct ECs produce differing codewords on at least one frame, and that an
  LTA element reproduces the SC output bit-for-bit (this is check 4 of §1.4,
  promoted to a regression test).

---

## 5. Validation gates

| gate | check | anchor |
|---|---|---|
| G1 | `admissible_matrix` on [P21]'s (16,7) | UT = `{(1,2)}`, all LT admissible |
| G2 | `blta_size`, `num_ec` | `S=(4,1,1,1,3)` ⇒ `E = 2205`; `S=(3,5)` ⇒ `E = 68355` |
| G3 | brute-force `Aut_aff` for `n = 3, 4` vs `|BLTA(S)|` | [P21, Table IV] for `N = 8` |
| G4 | LTA branch ≡ SC branch, bit-identical | done, §1.4 item 4 |
| G5 | `designed I` is decreasing | `is_decreasing` |
| G6 | BSGS group order (`N ≤ 256`) == `blta_size(S)` | §4.6 |
| G7 | AE-M BLER monotone non-increasing in `M` at fixed SNR | sanity |

---

## 6. Experiment matrix

`.ini` template = `project/POD/m7t10_itw2026/SC.ini` with
`permutation_src` **empty**, `Hmatrix_path` → the new polar `_H.matrix`,
`operationArray` = `'1' × n·(N/2)`, `bha_value_setting` = `'?' × N`.

| dir | code | design | runs |
|---|---|---|---|
| `n128_k100_highsnr` | (128,100) | DE/GA @ 10.5 dB | SC; AE32-SC with LTA / UTL / random; AE4-SCL8-UTL; SCL32 (5G) |
| `n1024_k512_U1` | (1024,512) | §4.3, upper-left | SC; AE8-SC |
| `n1024_k512_U2` | (1024,512) | §4.3, bottom-right | SC; AE32-SC |
| `n1024_k512_S41113` | (1024,512) | §4.4, `S=(4,1,1,1,3)` | AE8-SC `D=(4,3)`; AE32-SC `D=(2,2)`; AE8/32-SC random A; ML bound = SCL512 |
| `n256_k128_S35` | (256,128) | §4.4, `S=(3,5)` | AE8-SC; AE32-SC; SCL512 (ML bound) |
| `n256_k128_S53` | (256,128) | §4.4, `S=(5,3)` | AE32-SC (expected to be poor — the negative control) |

Expected reproductions: LTA gives *zero* gain under AE-SC; UTL/BLTA AE-32-SC
approaches the ML bound on the codes designed for it; the bottom-right UTL design
(`U2`) beats the upper-left one (`U1`); `S=(5,3)` stays poor.

---

## 7. MClassifier integration

`MClassDecoder` (`project/MClassifier/common/mclass_decoder_wrapper.h`) derives
from `AdjustPolarDecoderRelation` and only ever touches `received_order_set`,
`relation_ship`, `path_metric` and `Hmatrix` — all of which are populated
identically for a plain polar code. Every `mclass_*` binary takes `-ini`.
So integration is:

1. New `.ini` per code (`AE64SC` / `AE64SCL4` style, `M = 64`).
2. New target dirs under
   `project/MClassifier/adaptive_path_selector/pm/<polar_code_id>/` and
   `.../trace_learnt_path_selector/<polar_code_id>/`, add-only entries in
   `project/MClassifier/CMakeLists.txt` pointing at the *existing*
   `mclass_adaptive_m_dataset.cpp` / `mclass_adaptive_m_bler_sim.cpp` sources
   (no new C++ unless the feature set changes).
3. Regenerate datasets, retrain `train_adaptive_m.py`, rerun the sweeps.

**What is scientifically new here, and worth calling out.** In POD the `M`
branches are cosets `π·Aut(C_b)` of a *dynamic-frozen* transform; here they are
equivalence classes of a *static-frozen* polar code, and the ensemble has an
exact, finite, enumerable size `E` (§2.3). That makes two MClassifier questions
crisper than on eBCH:

* the adaptive-`m` target is defined against a *known-complete* ensemble
  (`M = E` is reachable for `S=(3,5)`-style codes only in principle, but the
  branch population is algebraically labelled — every branch carries its `(P,U)`
  descriptor);
* branches now have **structured features for free** — the `v`/`p` descriptor
  vectors and the Hamming distances of [P22, §III-D] — which is exactly the kind
  of per-branch side information the trace-learnt selector currently has to infer
  from path metrics. Feeding `(v, p)` in as branch features is the obvious first
  new experiment.

---

## 8. Known gaps / risks

1. **5G reliability sequence is not in the repo.** Needed for the `SCL32 - 5G`
   and `CA-SCL - 5G` baselines. Either add the 1024-entry 38.212 table as a data
   file or drop those baselines and compare against DE/GA-designed SCL.
2. **No CRC support in the decoder.** [P22] uses CA-SCL and a 6-bit CRC on some
   AE runs. Reproducing those needs a CRC concatenation in the decoder
   (candidate-side check in the combiner would be the least invasive place).
   Recommend deferring; the non-CRC comparisons are the ones that carry the
   paper's message.
3. **`operationArray` length.** `libParser.cpp` uses `value_name[10000]`.
   `n·(N/2)` = 5120 for `N = 1024` (fits) but 11264 for `N = 2048` (overflows).
   If the study ever goes to `N = 2048`, bump that buffer.
4. **Startup cost at `N = 1024`.** The constructor does two dense `N³`-ish
   `matrixMultiplication` calls over `char` (≈1.5·10⁹ ops). Tolerable once per
   run; if it becomes annoying, add an identity fast-path for
   `permutation_src == ""`.
5. **`[1]` may be larger than `BLTA(2,1,…,1)`** for some codes ([P22, Thm. 4]
   leaves the refinement `S₁` code-dependent). Gate G4 / `verify_aut.py`'s
   behavioural check catches this; the symptom is duplicate branch outputs and a
   flat BLER-vs-`M` curve.
6. **`A_G` for the design algorithms must preserve decreasingness.** [P22]'s
   Algorithm 1 is the version that guarantees it; [P21]'s §III-C design does not.
   Implement both, but use §4.4 for anything that feeds the BLTA claims.
7. **`end` in the `[AWGN]` section is effectively exclusive.**
   `AWGN::nextChannel()` stops when `start + step == end`
   (`(start+step-end)*(start-end) <= 0`), so `start=3.0 step=0.5 end=5.5`
   yields points 3.0 … 5.0 and *not* 5.5. Set `end` half a step past the last
   point you want (here `5.6`). Cost me one pass of the Fig. 3 sweep.
8. **GA stage order.** `ga_means()` must transform bit `n-1` first and bit `0`
   last — `(W^-)^+ != (W^+)^-`, and reversing the loop silently swaps the
   reliability of index pairs like 1 and 2. Gate `G-STAGE` in
   `tools/test_polar_design.py` pins this against the decoder's own
   Bhattacharyya recursion (which runs `stage-1 → 0`).

---

## 9. Phasing

| phase | deliverable | gates |
|---|---|---|
| P1 | `polar_design.py`, `monomial_aut.py`; `(16,7)` and `(16,10)` regression | G1, G2, G3, G5 |
| P2 | `blta_sample.py` + `verify_aut.py`; `n128_k100_highsnr` LTA-vs-UTL-vs-random AE32-SC | G4, G7 |
| P3 | `design_utl.py` → `U1`/`U2` at `(1024,512)`; reproduce [P21] Fig. 4 | — |
| P4 | `design_blta.py` → `S=(4,1,1,1,3)`, `(3,5)`, `(5,3)`; reproduce [P22] Figs. 3–5 | G2 |
| P5 | `blta_bsgs.py` cross-check at `N ≤ 256` | G6 |
| P6 | MClassifier on the best AE code (`M = 64`), incl. `(v,p)` branch features | — |

P1+P2 are the load-bearing ones: once `n128_k100_highsnr` shows the LTA/UTL split
in BLER, everything after it is code design plus bookkeeping.
