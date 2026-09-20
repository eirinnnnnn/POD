# Results — (128, 100) polar, high-SNR DE/GA design

Reproduction of C. Pillet, V. Bioglio, I. Land, *Polar Codes for Automorphism
Ensemble Decoding*, ITW 2021 (arXiv:2102.08250), **Fig. 3**.

Figure: `fig3_repro.png`. AWGN, BPSK, DE/GA design at 10.5 dB Eb/N0.
Every point is run to **100 block errors** (cap 600 000 frames, never reached
except where noted). All runs share the same message and channel seeds.

## BLER

| Eb/N0 [dB] | 3.00 | 3.50 | 4.00 | 4.50 | 5.00 |
|---|---|---|---|---|---|
| SC                    | 3.06e-01 | 1.56e-01 | 5.98e-02 | 2.12e-02 | 5.23e-03 |
| **AE32-SC  LTA**      | **3.06e-01** | **1.56e-01** | **5.98e-02** | **2.12e-02** | **5.23e-03** |
| AE32-SC  UTL          | 7.74e-02 | 2.41e-02 | 7.18e-03 | 1.33e-03 | 3.10e-04 |
| AE32-SC  BLTA (P·U)   | 7.67e-02 | 2.30e-02 | 6.46e-03 | 1.16e-03 | 2.72e-04 |
| AE32-SC  random A     | 7.80e-02 | 2.20e-02 | 6.00e-03 | 1.31e-03 | 3.00e-04 |
| AE4-SCL8  UTL         | 5.81e-02 | 1.53e-02 | 3.81e-03 | 7.64e-04 | 1.75e-04 |
| SCL32                 | 4.56e-02 | 1.35e-02 | 3.47e-03 | 7.57e-04 | 1.74e-04 |

Frame counts (all 100 errors):

| | 3.00 | 3.50 | 4.00 | 4.50 | 5.00 |
|---|---|---|---|---|---|
| SC / AE32-SC LTA | 327 | 642 | 1673 | 4723 | 19107 |
| AE32-SC UTL      | 1292 | 4156 | 13927 | 75497 | 323072 |
| SCL32            | 2193 | 7400 | 28837 | 132109 | 573084 |

## What reproduces

**1. LTA automorphisms give exactly zero AE gain.** The SC and AE32-SC-LTA rows
are not merely close — they are *byte-identical*, same error count and same
frame count at every SNR point (100/327, 100/642, 100/1673, 100/4723,
100/19107). All 32 LTA branches are absorbed by SC, so the ensemble decodes
32 copies of the same thing. This is [P21]'s central negative result, and it is
visible before simulating: `verify_aut.py` reports the LTA ensemble covers
**1 distinct equivalence class out of 32 branches**, while the UTL, P·U and
random-BLTA ensembles each cover 32/32.

**2. UTL automorphisms are not absorbed and close most of the SC→SCL gap.**
AE32-SC-UTL improves on SC by 4.0× in BLER at 3.0 dB growing to 16.9× at
5.0 dB, i.e. roughly **0.9 dB at BLER 5e-3**, at the latency of a single SC
decode. Against SCL32 it is left with a gap of about 0.25 dB.

**3. A small AE of list decoders matches a large list.** AE4-SCL8-UTL and SCL32
are statistically indistinguishable across the whole range
(1.75e-4 vs 1.74e-4 at 5.0 dB) — the same effective list size 4×8 = 32, but the
four branches are independent, so the AE version needs no inter-branch path
exchange.

## What is *not* resolved here

**UTL vs BLTA(P·U) vs random-A is inside the noise.** At 5.0 dB the three are
3.10e-4 / 2.72e-4 / 3.00e-4. With 100 errors per point the relative standard
error is ~10%, so a 14% spread is not a measurement. The ordering suggested by
[P22] — that BLTA equivalence-class representatives beat pure UTL because they
also use the block-permutation part `P` — would need ~1000 errors per point to
test. That is a cheap follow-up (≈10× the frames on the three AE32 curves only)
but it was not run.

## Deviations from the paper

1. **SNR range 3.0–5.0 dB**, not 3.0–6.0. SCL32 runs at ~95 frames/s in this
   implementation (`clonePage` copies the full `(stage+1)×N` page per path
   split, 32 paths × 100 information bits per frame), so 5.5 and 6.0 dB would
   have cost hours for the SCL reference alone while the AE-SC curves finished
   in minutes. Curves stop where all seven can be compared on equal footing.
   Extending: AE-SC points cost ~20 min each at 800 k frames; SCL32 ~2.5 h.
2. **No 5G-sequence baselines** (`SCL32 - 5G`, `SCL32-CRC11 - 5G`). The 3GPP
   38.212 reliability sequence is not in the repo and the decoder has no CRC
   support. The SCL32 curve here is on the *same* high-SNR code, which is the
   more direct comparison anyway.
3. **No `SC - Optimal SNR`, no FGA, no BP curves.** Per-SNR re-design, factor
   graph permutations and BP are outside what the existing decoder does.
4. **Combiner is the SC path metric, not exact least squares.** [P21]/[P22]
   select the branch minimising ‖y − x̂‖². `AdjustPolarDecoderRelation` selects
   the branch with the smallest SCL path metric, which under min-sum SC is a
   monotone *approximation* of that. It is the one place where this
   reproduction could be losing a little performance relative to the paper, and
   adding a true LS combiner is a ~20-line, opt-in change to
   `AED_relation_check.cpp` (the candidate codeword and `received` are both
   already in scope). Worth doing before quoting absolute gaps against [P22].

## Reproduce

```
cd tools
python3 polar_design.py --n 7 --K 100 --design ga --snr 10.5 \
    --out_dir ../codes/n128_k100_highsnr --id n128_k100_highsnr
for m in utl lta pu random; do
  python3 blta_sample.py --I_file ../codes/n128_k100_highsnr/n128_k100_highsnr_I.txt \
      --n 7 --mode $m --M 32 --seed 1 --out ../codes/n128_k100_highsnr/aut_${m}_M32.txt
done
python3 blta_sample.py --I_file ../codes/n128_k100_highsnr/n128_k100_highsnr_I.txt \
    --n 7 --mode utl --M 4 --seed 1 --out ../codes/n128_k100_highsnr/aut_utl_M4.txt
cd ../codes/n128_k100_highsnr && ./_run.sh          # ~40 min on 8 cores
cd ../../tools && python3 plot_bler.py --dir ../codes/n128_k100_highsnr \
    --runs SC AE32SC_LTA AE32SC_UTL AE32SC_PU AE32SC_RND AE4SCL8_UTL SCL32 \
    --out fig3_repro.png --ymin 5e-5
```

Validation gates (`tools/test_polar_design.py`), all passing:

```
G-STAGE : PASS  (max |Z_python - Z_decoder| = 3.25e-07)
G1      : PASS  (UT admissible = [(1,2)] on [P21] section III-A's (16,7) example)
G2      : PASS  S=(4,1,1,1,3) -> |EC| = 2205 ;  S=(3,5) -> |EC| = 68355
G3      : PASS  (brute-force |Aut_aff| == |BLTA(S)|, n <= 4)
G5      : PASS  (designed information sets are decreasing / UPO-compliant)
```

---

# Second experiment — constant effective list size M × L = 32

Figure: `fig_eff32.png`. Same code, same seeds, same 100-block-error stopping
rule. The question: at a fixed *effective* list size of 32, how does buying
diversity with parallel automorphism branches compare with buying it with SCL
list depth? SCL32 is the `M = 1` endpoint of this family and AE32-SC the
`L = 1` endpoint, so the whole trade-off is one curve family. OSD order 1
(exact-LS re-encoding, `project/OSD`) is included as an independent reference.

| (M, L) | 3.00 | 3.50 | 4.00 | 4.50 | 5.00 |
|---|---|---|---|---|---|
| SC (1, 1)        | 3.06e-01 | 1.56e-01 | 5.98e-02 | 2.12e-02 | 5.23e-03 |
| OSD order 1      | 8.64e-02 | 2.66e-02 | 7.84e-03 | 1.91e-03 | 4.24e-04 |
| SCL32 (1, 32)    | 4.56e-02 | 1.35e-02 | 3.47e-03 | 7.57e-04 | 1.74e-04 |
| AE2-SCL16        | 5.16e-02 | 1.46e-02 | 3.69e-03 | 7.57e-04 | 1.74e-04 |
| AE4-SCL8         | 5.81e-02 | 1.53e-02 | 3.81e-03 | 7.64e-04 | 1.75e-04 |
| AE8-SCL4         | 6.68e-02 | 1.67e-02 | 4.34e-03 | 8.01e-04 | 1.84e-04 |
| AE16-SCL2        | 7.57e-02 | 1.87e-02 | 5.04e-03 | 9.60e-04 | 2.05e-04 |
| AE32-SC (32, 1)  | 7.74e-02 | 2.41e-02 | 7.18e-03 | 1.33e-03 | 3.10e-04 |

## Findings

**1. BLER is monotone in list depth `L` at every SNR point.** At fixed
effective list size, depth beats branch diversity on this code. The penalty for
going all the way to `M = 32, L = 1` is about **0.2 dB** relative to SCL32
(measured at BLER 1e-3: SCL32 at 4.41 dB, AE32-SC at 4.60 dB).

**2. The knee is at M = 4, and up to it parallelism is free.** AE4-SCL8 is
statistically indistinguishable from SCL32 (7.64e-4 vs 7.57e-4 at 4.5 dB;
1.75e-4 vs 1.74e-4 at 5.0 dB, both within the ~10% standard error of a
100-error measurement). So four independent SCL8 decoders — a quarter of the
list depth, no inter-branch path exchange — buy the full SCL32 error rate.

**3. AE2-SCL16 is on the ML bound.** At *both* 4.5 and 5.0 dB it and SCL32
produced byte-identical logs (100 errors in exactly 132109 and 573084 frames).
Two decoders with different internal structure failing on exactly the same
frames means both are decoding to the maximum-likelihood codeword on every
frame that matters; the curve below `M = 4` is the ML bound, not a decoder
property.

**4. OSD order 1 is uniformly worse than every eff-32 configuration,**
including AE32-SC, at all five SNRs (4.24e-4 vs 3.10e-4 at 5.0 dB), while
beating plain SC by roughly an order of magnitude. For `K = 100` an order-1
reprocessing is simply too small a candidate list.

## How this contrasts with PED on eBCH

This is the **opposite** ordering from the PED/POD results on eBCH, where
PED64-SCL4 beat PED4-SCL64 by about 0.5 dB. A plausible explanation: in PED the
BLBC is embedded as a *dynamic-frozen* polar subcode whose reliability profile
is whatever the transformation parameter happens to give, so SC is weak and
transformation diversity recovers a lot. Here the frozen set is *designed*, SCL
is already essentially ML at `L = 16`, and there is nothing left for diversity
to recover. Worth testing directly rather than assuming: the same eff-32 sweep
on a code where SCL32 is *not* near-ML (lower rate, or the `S = (5,3)`
(256,128) code of [P22] Fig. 5) would separate "AE is weaker" from "this code
is too easy".

## Caveat specific to the high-M end

The combiner is the min-sum SC path metric, not the exact `‖y − x̂‖²` used by
[P21]/[P22] (and used by the OSD reference here). An approximate ranking metric
costs more the more candidates it has to rank, so it penalises `M = 32` more
than `M = 2`, and some unknown part of AE32-SC's 0.2 dB may be combiner loss
rather than a diversity limit. The exact-LS combiner is a ~20-line opt-in
addition to `AED_relation_check.cpp` and should be in place before this 0.2 dB
is quoted as a property of the scheme.

## Reproduce (second experiment)

The AE ensembles for M = 2, 8, 16 are nested prefixes of `aut_utl_M32.txt`
(same seed, and acceptance depends only on previously accepted elements), so
the M-sweep is a clean "add more branches" experiment with no ensemble-selection
confound. Verified by `blta_sample.py` output and a prefix check.

```
cd tools
for M in 2 8 16; do
  python3 blta_sample.py --I_file ../codes/n128_k100_highsnr/n128_k100_highsnr_I.txt \
      --n 7 --mode utl --M $M --seed 1 \
      --out ../codes/n128_k100_highsnr/aut_utl_M${M}.txt
done

# the OSD reference needs its own binary (not built by project/PolarAED)
cd ../../OSD && mkdir -p build && cd build && cmake .. && make   # -> OSD_performance

cd ../../PolarAED/codes/n128_k100_highsnr && ./_run3.sh          # ~1 h on 8 cores
cd ../../tools && python3 plot_bler.py --dir ../codes/n128_k100_highsnr --eff32 \
    --runs SC OSD1 SCL32 AE2SCL16_UTL AE4SCL8_UTL AE8SCL4_UTL AE16SCL2_UTL AE32SC_UTL \
    --title "$(printf '(128,100) polar, high-SNR design\nconstant effective list size  M x L = 32,  UTL ensembles')" \
    --out fig_eff32.png --ymin 1e-4
```
