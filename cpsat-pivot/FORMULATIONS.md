# Formulations — current version (2026-09-04)

All over `F_2`. `n = 2^m`, `G_p = F^{(x m)}`, `F = [[1,0],[1,1]]`.
Rows of `G_p`: `rho_a` (weight `2^|a|`, support = submasks of `a`).
Columns of `G_p`: `g_b`  (weight `2^(m-|b|)`, support = supermasks of `b`).
`|a| = popcount(a)`, `J : a -> n-1-a`.  `G_p^2 = I`; `G_p^T = J G_p J`.

Base code `C_b = rowspace(G_b)`, `dim k`.  `pi(G_b) := G_b[:, pi]`.
`W := pi(G_b) G_p`.  Target profile `p* = {p_0 < ... < p_{k-1}}`.

**Problem.**  Decide `∃ pi in Sym(n) : ColumnPivotProfile(W) = p*`.


## 1. Exact characterisation (what the solver encodes)

Three layers:

  (1)  W = pi(G_b) G_p                              [definition; F_2-linear in pi]
  (2)  span:  for every non-pivot j <= p_{k-1},
              W[:,j] in span{ W[:,p] : p in p*, p < j }
       (for j < p_0 this degenerates to W[:,j] = 0)
  (3)  independence:  ∃ U in GL(k,F_2) with U W[:,p*] = I_k

**Theorem.** (1) ∧ (2) ∧ (3)  <=>  ColumnPivotProfile(W) = p*.
Each of (2),(3) is necessary; only jointly are they sufficient.
(3) alone is far too weak — it says nothing about non-pivot columns.

Relaxation chain (feasible sets nested), = the `--certificate` levels in
`cpsat_pivot_reduced_20260902.py`:

    cuts-only  = (1) + {W[:,j]=0, j<p_0} + min-weight cuts     (Sec. 3)
    span-only  = (1) + (2)
    full       = (1) + (2) + (3)                               [the exact problem]

INFEASIBLE at ANY level  =>  p* unreachable (sound, one-directional).
SAT is meaningful ONLY at `full`; the relaxations cannot certify reachability.

Cost: (2) is `~k^2 * p_{k-1}` AND-gates, (3) is `k^3` AND-gates.  Both are
degree-2 over F_2 (unknown x unknown) — outside the XOR fragment, so Gaussian
elimination cannot touch them.  Everything else is linear.


## 2. The suffix flag (your D_j) — verified numerically

    D_j := span{ rho_j, ..., rho_{n-1} },   dim D_j = n - j
    E_j := D_j \ D_{j+1} = rho_j + D_{j+1},  |E_j| = 2^(n-1-j)

**Flag identity** (verified for all j at n=64 against the column pivot profile):

    dim( pi(C_b) ∩ D_j )  =  #{ p in p* : p >= j }
    =>   j in p*   <=>   pi(C_b) ∩ E_j != {}

**Duality lemmas** (both verified for all indices at n=64).
With `V_p := span{ g_0, ..., g_{p-1} }` (columns):

    L1.  V_j^perp = D_j
    L2.  V_p = J( D_{n-p} )        [column code = bit-reversed row code]

L2 removes any need for a separate "column formalism": `V_p` is a genuine
decreasing-monomial (polar) code up to the coordinate permutation `J`, so the
standard minimum-distance theorem applies, giving (verified):

    L3.  d_min(V_p) = 2^( m - nu(p) ),      nu(p) := max_{b < p} |b|


## 3. Min-weight cuts (redundant but propagation-friendly)

For `j < p_0` the profile forces `W[:,j] = 0`.  Expanding
`W[r,j] = XOR_a G_b[r,pi(a)] G_p[a,j]` and taking any codeword
`c = XOR_r alpha_r G_b[r,:]`:

    XOR_{a : G_p[a,j] = 1}  c[pi(a)]  =  0        for every c in C_b, j < p_0

**Support is column j of `G_p`.**  (An earlier version used `G_p^{-T}[:,j]`;
that is NOT valid — checked against a known witness permutation it was violated
311/600 times, vs 0/600 for the correct form.  It over-pruned and produced false
INFEASIBLE.  Fixed and re-verified two ways: 0/600 arithmetic violations, and
the full model with pi fixed to the witness returns OPTIMAL at every
certificate level.)

Encoding: `pw[a] = c[pi_a]` via one `AddElement` per position (no n^2
selectors), then one parity per `(c, j)` over `supp(G_p[:,j])`.
`N` low-weight codewords give `N * p_0` cuts; count is independent of `k`.
Logically implied by the `W[:,j]=0` rows already in (2) — these are *cuts* in
the ILP sense: they tighten propagation, not the feasible set.  Low weight is
chosen only for propagation strength.

Codeword generation: exhaustive `2^k` for `k <= 20`; bounded low-order
combination search over randomised RREF bases for `k > 20` (blind mask sampling
missed d_min entirely — reported 16 for eBCH[64,36] whose true d_min is 12).


## 4. pi-free necessary conditions  (`weight_stratum_match_20260904.py`)

**(P) Prefix / span condition.**  "no pivot before p_0"
  <=> pi(C_b) ⊆ D_{p_0}
  <=> V_{p_0} ⊆ pi(C_b^perp)                                   [L1]
  =>  d_min(V_{p_0}) >= d^perp := d_min(C_b^perp)

  With L3 and "least index of popcount r is 2^r - 1", closed form:

      p_0  <=  2^( m - ceil(log2 d^perp) + 1 ) - 1

  eBCH[64,16], m=6, d^perp=6:  p_0 <= 15.  Sharp (d_min(V_15)=8, d_min(V_16)=4).

**(F) Forbidden index.**  supp(A_{E_j}) ∩ supp(A_{C_b}) = {}  =>  j not in p*.
  eBCH[64,16]: forbidden = {0, 32, 60}.

**(C) Stratum capacity.**  Stratum `E_{p_i}` must host exactly `2^(k-1-i)`
  codewords, each of a weight `C_b` actually has:
      sum_{w in supp A_C} A_{E_{p_i}}(w)  >=  2^(k-1-i)

**(T) Transport / exact weight match.**  With
  `x[i,w] = #{c in pi(C_b) ∩ E_{p_i} : wt(c)=w}`:
      sum_i x[i,w] = A_{C_b}(w),  sum_w x[i,w] = 2^(k-1-i),
      0 <= x[i,w] <= A_{E_{p_i}}(w)
  Feasibility necessary; Gale/Hoffman cut over weight subsets S.


## 5. Weight enumerators

  A_{C_b}       : 2^k Gray-code enumeration.
  A_{C_b^perp}  : MacWilliams from A_{C_b} (never touch the 2^(n-k) dual).
  A_{D_j}       : exact Plotkin recursion, O(m n^2)  — verified against brute
                  force for all j = 43..64:

      j >= n/2 :  D_j = {(u,u) : u in D^(m-1)_{j-n/2}}
                  A_{D_j}(2w) = A^(m-1)_{D_{j-n/2}}(w)

      j <  n/2 :  D_j = {(u+v, v) : u in D^(m-1)_j, v free}
                  W_{D_j}(x) = sum_{u} (1+x^2)^(n/2 - wt u) (2x)^(wt u)

  A_{E_j}       = A_{D_j} - A_{D_{j+1}}   (D_j = D_{j+1} disjoint-union E_j)

  This settles the notes' "probably polynomial time (a guess)" affirmatively and
  constructively: every stratum enumerator is computable, including E_31 with
  its 2^32 elements.


## 6. Status on eBCH[64,16] / F^(x6)

Stratum-preserving family: all 22 indices of popcount >= 4 forced... in this
instance the anchor histogram is {4:9, 5:6, 6:1}, so 7 forced (popcount >= 5)
plus choose 9 of the 15 popcount-4 indices  =>  C(15,9) = 5005.
B_T = Z-sum of the equation profile; 5001 candidates below it.

  pi-free kills          : 3718 / 5001   (74.3%)   [ (P) 715 only, (F) 1716 only,
                                                     both 1287 ]
  survivors              : 1283,  Z-sum in [0.9725, 2.312]
  (beta / 5G-optimal, Z-sum 0.294, is killed by (F): index 60)

  solver on survivors    : `full` certificate, ~50 s each, ALL INFEASIBLE so far.

**(C) and (T) never fire**, even with the complete recursion data — the low
strata are far too large for capacity to bind, and transport budgets dwarf
demands.  The weight/stratum program caps at 74% on this instance; the
remaining obstruction is not visible to weight enumerators.

Subcode-embedding condition (correct but not usable):
  for every j there must exist S <= C_b with dim S = #{p >= j} and
  A_S(w) <= A_{D_j}(w).  Exact check is ~ Gauss-binomial(16,8)_2 ≈ 2^64;
  the cheap counting shadow `2^{k_j} - 1 <= sum_{w in supp A_{D_j}} A_{C_b}(w)`
  is vacuous here (pool = all 65535 codewords for every j <= 48).


## 7. Known limits

  * SAT-finding (witness search) works at n=32 (reachable profile, 17 s, no
    hint) but times out at n=64 (>900 s, even with the witness as a hint).
    The model is fine — with pi FIXED to the witness the full model returns
    OPTIMAL in 0.7 s.  The wall is branch-and-search over Sym(n).
  * Consequence: at n=64 the sweep's guarantee is "no candidate below B_T is
    provably reachable".  It can prove optimality when every verdict is
    INFEASIBLE; it cannot discover a reachable one.
  * CP-SAT `INFEASIBLE` is not independently certified (no DRAT).  Key results
    were cross-validated across three encodings instead.


## 8. Single-kernel pruning  (`prune_search_20260904.py`)

Not a SAT problem: with `P` fixed there is no existential quantifier, only
`m N / 2` deterministic linear-algebra evaluations.  Pure GF(2) + numpy.

    L_t          = layer operator, pairs indices differing in bit a_t (MSB=a_1)
    T            = L_m L_{m-1} ... L_1     (unpruned: = F^(x m), asserted)
    prune (t,q)  = delete the single entry L_t[a1, a0]     (rank one)
    M            = G_b P T^{-1},  RREF -> J = info set, non-pivots -> DF relations
    Phi(T, D)    = sum_{i in J} Z_i(T)

Z under pruning: propagate from the channel inward, layer 1 first; unpruned
kernel `(x,y) -> (x+y-xy, xy)`, pruned kernel `(x,y) -> (x,y)`.  Unpruned this
reproduces the standard sequence exactly (regression-tested), and the n=64
baseline `J` reproduces the sweep's anchor profile
`[15,23,27,29,30,31,39,43,46,47,51,55,59,61,62,63]`.

**Conservation law (exact).**  Both kernel maps have output sum `x+y`, and each
layer is a disjoint union of them, so

    sum_i Z_i(T) = N z_0     for EVERY pruning pattern.

Pruning cannot create reliability; it only redistributes a fixed budget between
`J` and `J^c`.  Hence `Phi = N z_0 - sum_{i not in J} Z_i`, and a pruning is
neutral (`delta_e = 0` exactly) whenever its downstream cone lies entirely
inside `J` or entirely inside `J^c` — this accounts for all the exact zeros
(14/80 at n=32, 36/192 at n=64).

**Result (lambda = 0).**  Over 5 codes x 3 SNRs x every single pruning, plus
400 random multi-prunings of 2..12 kernels at n=64:

    forced-J improvements : 0 / all      min delta = +0.000000
    free-J   improvements : 16 / 192 at n=64,  min delta = -0.009907

So the *transform* is not the bottleneck.  With a free information set pruning
can help, barely; with `J` forced by `RREF(G_b P T^{-1})` it never does.  The
gap that matters is alignment, not polarization:

    n=64 k=16, 3 dB:  Phi = 2.3160   vs   free-J optimum 0.2937   (gap 2.02)
    n=32 k=16, 3 dB:  Phi = 0.5397   vs   free-J optimum 0.4277   (gap 0.11)

Note this is NOT explained by majorization: `Z(T_full) > Z(T_pruned)` in the
majorization order fails for 80/192 single prunings at n=64.  The local
de-polarization is real per kernel but does not survive the later layers; the
`delta_e >= 0` observation is empirical, not proved.

**Consequence for the search.**  The standard-kernel `P` is a stationary point
of the basic objective under `R_1` — every single-layer pruning has zero or
positive influence, so greedy/beam search over `R_1` has nothing to descend.
Pruning only becomes an active degree of freedom under the
dynamic-frozen-informativeness objective
`Phi - lambda sum_f Delta_f`, where improvements do appear (e.g. 18/192 at
n=64, 3 dB, lambda=0.05).  But that term is currently unbounded — it drives
`Phi` negative — so those are not yet meaningful optima; the term needs
calibration before a search over it means anything.
