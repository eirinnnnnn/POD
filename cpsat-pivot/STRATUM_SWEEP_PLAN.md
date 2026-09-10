# Stratum-preserving reachability sweep — eBCH [64,16] / F^⊗6

## Assignment

Decide whether Trifonov's equation-constructed info set is **Bhattacharyya-optimal
among all reachable pivot profiles that preserve its weight-stratum profile**, and
if not, return the optimal one.

"Reachable" = ∃ coordinate permutation π with `PivotProfile(π(G_b)·G_p) = p_star`.

## Why the family is small (5005, not C(64,16)≈5e14)

Both the equation set and the β (reliability) set have the **same stratum profile**:

| stratum | indices (popcount) | count | in family |
|---|---|---|---|
| weight ≥ 5 | {31,47,55,59,61,62,63} | 7 | **all forced** |
| weight = 4 | {15,23,27,29,30,39,43,45,46,51,53,54,57,58,60} | 15 | **choose 9** |

The RM/degree criterion has no opinion inside a stratum; the equation resolves the
tie algebraically, β resolves it by reliability. So the whole question is:

> **which 9 of the 15 weight-4 indices** → `C(15,9) = 5005` candidate profiles.

## Anchors (Z-sum = Σ Bhattacharyya, AWGN recursion, design SNR 3 dB, R=1/4)

| set | weight-4 choice | Z-sum | p_star[0] | reachable? |
|---|---|---|---|---|
| equation | {15,23,27,29,30,39,43,46,51} | **2.316** | 15 | **YES** (POD permutation) |
| β | {30,45,46,51,53,54,57,58,60} | **0.303** | 30 | unknown — the sweep decides |

**5001 of 5005 candidates beat the equation's Z-sum.** The equation's algebraic
tie-break lands on a near-worst reliability choice. So the sweep is *not* "which few
beat the equation" — it is "of the ~5000 that beat it, which is the most reliable
one that is actually reachable".

## Code facts

- `eBCH_m6_t11.matrix` [64,16], **self-orthogonal**: hull dim h = k = 16
  ⇒ pair(π) ≤ k−h = 0: no complementary pair {j,63−j} may both be pivots.
  *Vacuous within this family* (weight-4 ↔ weight-2, weight≥5 ↔ weight≤1 under
  complement), but it confirms `D_j^⊥ = D_{n-j}` holds for F^⊗6 here — both anchor
  sets are pair-free, consistent with the notes' self-orthogonal theorem.
- Favorable regime: p_star[0] ∈ [15,30] ⇒ the min-weight prefix cuts are **strong**
  (contrast the eGolay tail where p_star[0]≈0 and the cuts were useless).
  A more-reliable (β-like) selection drops the low weight-4 indices ⇒ larger
  p_star[0] ⇒ faster UNSAT proofs. The ordering works *with* us.

## Pipeline — per candidate, in ascending Z-sum order

0. **score** Z-sum at `--snr`; sort ascending; the only candidates that matter are
   those with Z-sum < B_T = Z-sum(equation).

1. **π-free prefilters** (microseconds, sound necessary conditions):
   - **Λ-envelope**: `A(i) ≤ min(Λ(δ_i), k, n−i)` for all i, where
     `A(i) = #{p ∈ p_star : p ≥ i}`, `δ_i = min_{j≥i} 2^popcount(j)`,
     `Λ(d) = max dim of a subspace of C_b whose nonzero words all have weight ≥ d`.
     (plug in the implementation from the parallel `ebch.py`; fall back to skip
      if Λ can't be computed at the needed dimension.)
   - hull pair bound (vacuous here, kept for generality).

2. **span-only prefilter** — `solve_reduced(certificate="span-only",
   min_weight_words=N)`, short cap (~60 s). `span-only` INFEASIBLE ⇒ provably
   unreachable (sound relaxation). Skip.

3. **full SAT** — `solve_reduced(certificate="full", min_weight_words=N)`,
   generous cap. SAT ⇒ record π (verified over F₂). INFEASIBLE ⇒ record.
   UNKNOWN ⇒ escalate to `cms_pivot_reduced_20260903` (native-XOR cuts + Gauss).

4. **stop** at the first full-SAT — by ascending order it is provably optimal
   within the stratum-preserving family.

## Outcomes

- first SAT with Z-sum < B_T → **a strictly better reachable profile**; if it is β,
  β wins the theorem.
- every candidate with Z-sum < B_T is INFEASIBLE and the equation is SAT →
  **the equation is optimal within the stratum-preserving family** (most of the
  theorem; the remaining gap is other stratum profiles).
- residual UNKNOWNs below B_T → not conclusive; list for a longer / CMS pass.

## Compute estimate

- p_star[0] ≥ 15 ⇒ UNSAT proofs in the seconds–minute range (cf. the n=64 5G
  target: INFEASIBLE in 2.8 s).
- Worst case (equation is the answer) ≈ 5000 UNSAT proofs; best case (β reachable)
  = a handful. Realistically hours to ~1 day, resumable via JSON cache.
- Does **not** need n=64 SAT-finding to scale — the sweep produces its first real
  answer while the general encoding work continues.

## Files

- driver: `stratum_sweep_20260904.py`
- solver: `cpsat_pivot_reduced_20260902.py` (`--w-encoding staged`, `--certificate`)
- escalation: `cms_pivot_reduced_20260903.py`
- cache: `stratum_sweep_cache.json`
