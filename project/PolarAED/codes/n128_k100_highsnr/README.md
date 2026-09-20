# (128, 100) polar code, high-SNR DE/GA design — reproduction of [P21] Fig. 3

Reference: C. Pillet, V. Bioglio, I. Land, *Polar Codes for Automorphism
Ensemble Decoding*, ITW 2021 (arXiv:2102.08250), Fig. 3.

The paper's point: **LTA automorphisms are absorbed by SC and give exactly zero
AE gain**, while UTL automorphisms are not absorbed and let AE-SC approach SCL
performance at SC latency.

## Code

```
cd ../../tools
python3 polar_design.py --n 7 --K 100 --design ga --snr 10.5 \
    --out_dir ../codes/n128_k100_highsnr --id n128_k100_highsnr
python3 monomial_aut.py --n 7 --I_file ../codes/n128_k100_highsnr/n128_k100_highsnr_I.txt
```

Design is DE/GA at 10.5 dB Eb/N0 — the "high-SNR design" of [P21] §IV, chosen
because it naturally admits many UTL automorphisms ([P21] Fig. 2).

```
decreasing / UPO-compliant : True
admissible UT positions    : 11  [(0,1) (0,2) (0,3) (0,4) (1,2) (1,3) (1,4)
                                  (2,3) (2,4) (3,4) (5,6)]
block structure S          : (5, 2)
Aut_aff == BLTA(S) exactly : True
|BLTA(S)|                  : 7863816683520
[1] profile (SC-absorbed)  : (2,1,1,1,1,1)   |[1]| = 805306368
equivalence classes E      : 9765
|A_U| = 2048, |A_P| = 240, pure-UTL elements = 2048
```

`E = 9765` ≫ 32, so an AE-32 ensemble is far from exhausting the code's
diversity.

## Ensembles

```
for m in utl lta pu random; do
  python3 blta_sample.py --I_file .../n128_k100_highsnr_I.txt --n 7 \
      --mode $m --M 32 --seed 1 --out .../aut_${m}_M32.txt
done
python3 blta_sample.py ... --mode utl --M 4 --out .../aut_utl_M4.txt
```

Every file has the identity as row 0. `verify_aut.py` confirms all 32
permutations in each file are genuine automorphisms, and reports the number of
distinct equivalence classes covered:

| ensemble | automorphism check | distinct ECs / branches |
|---|---|---|
| `aut_utl_M32`    | ALL PASS | 32 / 32 |
| `aut_pu_M32`     | ALL PASS | 32 / 32 |
| `aut_random_M32` | ALL PASS | 32 / 32 |
| `aut_lta_M32`    | ALL PASS | **1 / 32** |

The LTA ensemble collapsing to a single equivalence class is the paper's result,
visible *before* any simulation: all 32 LTA branches are SC-absorbed.

## Runs

```
./_run.sh          # 7 configs in parallel; logs land in <run>/log.txt
python3 ../../tools/plot_bler.py --dir . \
    --runs SC SCL32 AE32SC_LTA AE32SC_UTL AE32SC_PU AE32SC_RND AE4SCL8_UTL \
    --title "(128,100) polar, high-SNR DE/GA design" --out fig3_repro.png
```

| run | decoder | ensemble |
|---|---|---|
| `SC`           | SC (L=1)  | — |
| `SCL32`        | SCL L=32  | — |
| `AE32SC_LTA`   | AE-32 SC  | 32 LTA elements (negative control) |
| `AE32SC_UTL`   | AE-32 SC  | 32 UTL elements on admissible positions |
| `AE32SC_PU`    | AE-32 SC  | 32 BLTA EC representatives `A = P·U` |
| `AE32SC_RND`   | AE-32 SC  | 32 uniform draws from BLTA(S) |
| `AE4SCL8_UTL`  | AE-4 SCL8 | 4 UTL elements |

All runs share `seed_string = 111511015` and the same AWGN seed, so identical
decoders produce byte-identical logs — which is exactly how the LTA absorption
shows up.

## Result

`fig3_repro.png`. Headline numbers are in `RESULTS.md`.
