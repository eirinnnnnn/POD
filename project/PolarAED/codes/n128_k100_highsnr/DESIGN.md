# (128, 100) high-SNR DE/GA design — frozen set

Produced by `tools/polar_design.py --n 7 --K 100 --design ga --snr 10.5`
(DE/GA, Eb/N0 = 10.5 dB, the "high-SNR design" of [P21] §IV).

Index convention (PLAN.md §1.2): `F^{⊗7}` in **natural order, no bit reversal**,
`polar_matrix[r][c] = 1 ⟺ c ⊆ r`. Row `r` carries the monomial containing
`x_i` **iff bit i of r is 0**, so `r = 127` is the constant monomial (weight-128
row, most reliable) and `r = 0` is `x₀x₁x₂x₃x₄x₅x₆` (weight 1, least reliable).
Monomial degree = `7 − popcount(r)`.

## Frozen set (28 positions)

```
0 1 2 3 4 5 6 8 9 10 12 16 17 18 20 24 32 33 34 36 40 48 64 65 66 68 72 80
```

Frozen/information map, `1` = frozen, `0` = information, index 0 → 127:

```
  0: 11111110 11101000 11101000 10000000
 32: 11101000 10000000 10000000 00000000
 64: 11101000 10000000 10000000 00000000
 96: 00000000 00000000 00000000 00000000
```

## The design in closed form

```
I = { r : popcount(r) ≥ 3 }  ∪  { 96 }
```

`{ r : popcount(r) ≥ 3 }` is exactly **RM(4, 7)**, the monomials of degree ≤ 4,
of dimension `C(7,0)+C(7,1)+C(7,2)+C(7,3)+C(7,4) = 99`. Index `96 = 1100000₂`
has popcount 2, i.e. degree 5, and is the monomial `x₀x₁x₂x₃x₄`.

> **The high-SNR (128,100) polar code is RM(4,7) plus the single degree-5
> monomial `x₀x₁x₂x₃x₄`.**

Verified:

```
I == {popcount >= 3} | {96}   ->  True
I \ RM(4,7) = [96] ,  RM(4,7) \ I = []
```

This is why a 10.5 dB design SNR was used: pushing the design SNR up drives the
information set toward a Reed–Muller code, and RM codes have the largest
possible affine automorphism group ([P21] Fig. 2 — the maximum of 21 admissible
upper-triangular positions at N = 128 is attained by the eight RM codes).

## Why the automorphism group is BLTA(5, 2)

The one non-RM generator, `x₀x₁x₂x₃x₄`, is **symmetric in `{x₀,…,x₄}` and
involves neither `x₅` nor `x₆`**. Any invertible substitution among `x₀…x₄`
maps it to itself plus lower-degree terms already in RM(4,7), and RM(4,7) itself
is invariant under all of GA(7). The admissible-position map computed by
`tools/monomial_aut.py` shows exactly that shape:

```
    D1111..        rows/cols 0-4: free 5x5 block  ->  s1 = 5
    1D111..
    11D11..
    111D1..
    1111D..
    11111D1        rows/cols 5-6: free 2x2 block  ->  s2 = 2
    111111D
```

(`D` = diagonal, `1` = admissible, `.` = not admissible.)

```
block structure S          : (5, 2)
Aut_aff == BLTA(S) exactly : True
|BLTA(S)|                  : 7 863 816 683 520
[1] (SC-absorbed)          : BLTA(2,1,1,1,1,1),  |[1]| = 805 306 368
equivalence classes E      : 9 765
|A_U| = 2048,  |A_P| = 240,  pure-UTL elements = 2^11 = 2048
admissible UT positions    : 11
```

`E = 9765` is the maximum useful ensemble size; every AE experiment here uses
`M ≤ 32`, so the code's diversity is nowhere near exhausted.

## Contrast: why longer codes need a design algorithm

The same DE/GA rule at (1024, 512) @ 3.0 dB yields `S = (1,1,1,1,1,1,1,1,1,1)`
— **zero** admissible upper-triangular positions, hence no UTL automorphisms and
no possible AE-SC gain. That is [P21] Table III's observation and the reason
phases P3/P4 of PLAN.md need the explicit §III-C / Algorithm 1 designs rather
than a plain reliability sequence.
