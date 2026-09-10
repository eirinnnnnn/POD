#!/usr/bin/env python3
"""
Stratum-preserving reachability sweep, generalised.  (See STRATUM_SWEEP_PLAN.md.)

Given an eBCH generator, the k-ignorant "equation" permutation P, and a polar
transform F^{x m}:

  * anchor profile  p_t  = PivotProfile( pi(G_b) F^{x m} )  with pi[i]=argmax P[i,:]
    (reachable by construction).
  * stratum family  = all pivot profiles with the SAME popcount(index) histogram
    as p_t: every fully-occupied popcount class is forced, the one partially
    occupied class is the choose-set.
  * B_T = Bhattacharyya sum of p_t (design SNR).  The 5G info set has the same
    histogram, so it is one member of the family.

Walk the family in ASCENDING Bhattacharyya sum; per candidate:
  1. span-only + min-weight cuts, short cap   (span-only UNSAT => unreachable)
  2. full certificate + min-weight cuts, generous cap
  3. UNKNOWN -> flag for cms_pivot_reduced_20260903 escalation
Stop at the first full SAT (= optimal reachable profile within the family), or at
B_T with everything below it INFEASIBLE (=> p_t optimal within the family).

    python stratum_sweep_gen_20260904.py --code m6_t5 --snr 3.0 --workers 8
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from typing import List

import numpy as np

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, gf2_rank, polar_matrix, column_pivot_profile_gf2,
)
from cpsat_pivot_reduced_20260902 import solve_reduced

HERE = os.path.dirname(os.path.abspath(__file__))
POD = os.path.join(HERE, os.pardir, "project", "POD")
PTP = os.path.join(HERE, os.pardir, "project", "PtP")

CODES = {                                   # name -> (generator, P, m)
    "m6_t11": ("eBCH_m6_t11.matrix", "eBCH_m6_t5_P.matrix", 6),   # [64,16]
    "m6_t5":  ("eBCH_m6_t5.matrix",  "eBCH_m6_t5_P.matrix", 6),   # [64,36]
    "m6_t2":  ("eBCH_m6_t2.matrix",  "eBCH_m6_t5_P.matrix", 6),   # [64,51]
    "m5_t3":  ("eBCH_m5_t3.matrix",  "eBCH_m5_t3_P_5g.matrix", 5),
}


def parse_mat(p):
    t = open(p).read().split(); r, c = int(t[0]), int(t[1]); body = t[2:]
    if len(body) == r:
        return np.array([[int(x) for x in row] for row in body], np.uint8)
    return np.array([int(x) for x in body], np.uint8).reshape(r, c)


def bha_seq(m: int, ebn0_db: float, rate: float) -> List[float]:
    z = [math.exp(-rate * 10 ** (ebn0_db / 10))]
    for _ in range(m):
        z = [v for x in z for v in (2 * x - x * x, x * x)]
    return z


def five_g_info_set(n: int, k: int) -> List[int]:
    sys.path.insert(0, PTP)
    from static_ptp_5GNR import NR_5G_POLAR_SEQUENCE_1024
    q = [i for i in NR_5G_POLAR_SEQUENCE_1024 if i < n]
    return sorted(q[-k:])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", choices=sorted(CODES), default="m6_t5")
    ap.add_argument("--snr", type=float, default=3.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cuts-time", type=float, default=90.0)
    ap.add_argument("--span-time", type=float, default=120.0)
    ap.add_argument("--full-time", type=float, default=900.0)
    ap.add_argument("--min-weight-words", type=int, default=40)
    ap.add_argument("--max-cand", type=int, default=10**9)
    ap.add_argument("--max-hours", type=float, default=1e9)
    args = ap.parse_args()

    gen_f, perm_f, m = CODES[args.code]
    Gb = gf2(parse_mat(os.path.join(POD, gen_f)))
    Gp = gf2(polar_matrix(m))
    n = Gb.shape[1]; k = Gb.shape[0]
    assert gf2_rank(Gb) == k
    P = gf2(parse_mat(os.path.join(POD, perm_f)))
    pi_eq = [int(np.argmax(P[i, :])) for i in range(n)]
    p_t = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi_eq], Gp))
    h = k - gf2_rank(gf2_matmul(Gb, Gb.T))
    cache_path = os.path.join(HERE, f"stratum_sweep_{args.code}_cache.json")
    log_prefix = f"stratum_{args.code}"

    # popcount strata
    cls = {}
    for i in range(n):
        cls.setdefault(bin(i).count("1"), []).append(i)
    hist = {}
    for j in p_t:
        hist[bin(j).count("1")] = hist.get(bin(j).count("1"), 0) + 1
    forced = sorted(x for pc, cnt in hist.items()
                    if cnt == len(cls[pc]) for x in cls[pc])
    partial = [pc for pc, cnt in hist.items() if 0 < cnt < len(cls[pc])]
    if not partial:
        print(f"[sweep] {args.code}: every popcount class in p_t is FULLY occupied "
              f"-> the stratum family is a single point.")
        print(f"[sweep] p_t = {p_t}")
        print(f"[sweep] histogram {dict(sorted(hist.items()))}  (class sizes "
              f"{ {q: len(v) for q, v in sorted(cls.items())} })")
        print("[sweep] => nothing to search: p_t is trivially OPTIMAL within its "
              "stratum family.")
        return 0
    assert len(partial) == 1, f"expected exactly one partial stratum, got {partial}"
    pc = partial[0]
    pool = sorted(cls[pc])
    take = hist[pc]
    eq_choice = set(x for x in p_t if bin(x).count("1") == pc)
    p5g = five_g_info_set(n, k)
    fam = math.comb(len(pool), take)

    Z = bha_seq(m, args.snr, k / n)
    B_T = sum(Z[i] for i in p_t)

    print(f"[sweep] {args.code}: eBCH[{n},{k}]  hull dim h={h}  (pair(pi) <= {k-h})")
    print(f"[sweep] anchor p_t stratum histogram {dict(sorted(hist.items()))}")
    print(f"[sweep] forced ({len(forced)}): {forced}")
    print(f"[sweep] choose {take} of {len(pool)} popcount-{pc} indices  -> family {fam}")
    print(f"[sweep] equation choice: {sorted(eq_choice)}")
    print(f"[sweep] 5G info set: {p5g}   (same histogram: "
          f"{sorted(bin(j).count('1') for j in p5g) == sorted(bin(j).count('1') for j in p_t)})")
    print(f"[sweep] B_T (equation Z-sum @ {args.snr}dB) = {B_T:.6g}")
    print(f"[sweep] 5G Z-sum = {sum(Z[i] for i in p5g):.6g}", flush=True)

    cands = []
    for choice in itertools.combinations(pool, take):
        S = sorted(forced + list(choice))
        cands.append((sum(Z[i] for i in S), S, set(choice)))
    cands.sort(key=lambda x: x[0])
    below = sum(1 for s, _, _ in cands if s < B_T - 1e-15)
    print(f"[sweep] {below} / {fam} candidates below B_T; "
          f"most-reliable Z-sum {cands[0][0]:.6g}\n", flush=True)

    # ------------------------------------------------------------------
    # pi-free PREFIX-EMBEDDING THEOREM
    #   W[:,j]=0 for j<p0  <=>  V_p0 := span{g_0..g_{p0-1}}  subset of  pi(C_b^perp)
    #   permutation preserves weight  =>  d_min(V_p0) >= d_min(C_b^perp)
    #   (and more generally A_{V_p0}(w) <= A_{C_b^perp}(w) for all w)
    # Violation => unreachable for EVERY pi.  No solver.
    # ------------------------------------------------------------------
    from math import comb
    import collections as _co

    def _spec(B):
        kk = B.shape[0]; rows = [int("".join(map(str, r)), 2) for r in gf2(B)]
        c = _co.Counter({0: 1}); cur = 0
        for msk in range(1, 1 << kk):
            cur ^= rows[(msk & -msk).bit_length() - 1]
            c[bin(cur).count("1")] += 1
        return c

    A_C = _spec(Gb)                                   # 2^k words, enumerable
    A_D = {}                                          # MacWilliams -> dual spectrum
    for w in range(n + 1):
        A_D[w] = sum(A_C.get(v, 0) *
                     sum((-1) ** j * comb(v, j) * comb(n - v, w - j)
                         for j in range(0, w + 1))
                     for v in A_C) // (1 << k)
    d_dual = min(w for w, a in A_D.items() if w > 0 and a > 0)

    def light_word(p0, thresh, max_order=5):
        """
        Return (weight, order) of the lightest codeword of V_p0 = span{g_0..g_{p0-1}}
        found among XORs of <= max_order basis columns, stopping early if one of
        weight < thresh appears.  An UPPER bound on d_min(V_p0) -- which is the
        sound direction: exhibiting one word of weight < d_min(C^perp) proves
        V_p0 cannot embed into any permutation of C_b^perp.
        """
        cols = [int("".join(map(str, Gp[:, j].tolist())), 2) for j in range(p0)]
        best, bo = n + 1, 0
        for order in range(1, max_order + 1):
            for cb in itertools.combinations(range(p0), order):
                v = 0
                for c in cb:
                    v ^= cols[c]
                if v:
                    w = bin(v).count("1")
                    if w < best:
                        best, bo = w, order
                        if best < thresh:
                            return best, bo
        return best, bo

    p0_vals = sorted({S[0] for _, S, _ in cands})
    dminV = {p0: light_word(p0, d_dual) for p0 in p0_vals}

    from weight_stratum_match_20260904 import WeightStratumMatch
    WSM = WeightStratumMatch(Gb, Gp, cap=20, verbose=True)

    def theorem_kills(p0):
        """(P) prefix: non-None => V_p0 cannot embed in any pi(C_b^perp)."""
        w, o = dminV[p0]
        if w < d_dual:
            return f"(P) wt-{w} word in V_{p0} < d_min(C^perp)={d_dual}"
        return None

    def pi_free_kills(S):
        """All pi-free necessary conditions: prefix (P) + weight stratum match."""
        return theorem_kills(S[0]) or WSM.violation(S)

    print(f"[theorem] d_min(C_b^perp) = {d_dual}   (exact, MacWilliams)")
    for p0 in p0_vals:
        r = theorem_kills(p0)
        w, o = dminV[p0]
        print(f"[theorem] p0={p0:3d}: lightest word found wt {w:3d}  "
              f"{'KILLS  (' + r + ')' if r else '>= d_min(C^perp) -> solver'}")
    n_kill = sum(1 for _, S, _ in cands if pi_free_kills(S))
    print(f"[pi-free] prefix + weight_stratum_match settle {n_kill} / {len(cands)} "
          f"candidates with no solver call\n", flush=True)

    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}

    def key(S): return ",".join(map(str, S))
    def save(): json.dump(cache, open(cache_path, "w"), indent=0)

    def solve(S, cert, tlim):
        return solve_reduced(Gb, Gp, S, w_encoding="staged", certificate=cert,
                             min_weight_words=args.min_weight_words,
                             time_limit_sec=tlim, workers=args.workers, verbose=False)

    t0 = time.time()
    n_inf = n_unk = 0
    for idx, (ssum, S, choice) in enumerate(cands[:args.max_cand], 1):
        if ssum >= B_T - 1e-15:
            print(f"\n[sweep] reached B_T at candidate {idx}: "
                  f"{n_inf} infeasible, {n_unk} unknown below B_T.")
            print("[sweep] => p_t (the equation set) is OPTIMAL within the stratum family"
                  if n_unk == 0 else
                  "[sweep] => NOT conclusive: unknowns remain below B_T (see cache)")
            return 0
        if (time.time() - t0) / 3600 > args.max_hours:
            print(f"\n[sweep] --max-hours reached at candidate {idx}/{below}; "
                  f"{n_inf} infeasible, {n_unk} unknown so far.  resumable via cache.")
            return 0

        prev = cache.get(key(S))
        if prev and prev["status"] in ("infeasible", "sat", "unknown"):
            if prev["status"] == "sat":
                print(f"\n[sweep] cached SAT at {S} -> optimal.")
                return 0
            # 'unknown' entries are escalated separately (cadical_batch);
            # a fast-resume pass must not re-spend solver budget on them.
            continue

        # --- theorem prefilter: zero-cost, sound, no solver ---
        why = pi_free_kills(S)
        if why is not None:
            cache[key(S)] = {"status": "infeasible", "sumZ": ssum,
                             "tier": "pi-free", "why": why}
            save(); n_inf += 1
            continue

        tag = "  [=equation]" if choice == eq_choice else \
              ("  [=5G]" if choice == set(x for x in p5g if bin(x).count("1") == pc) else "")
        print(f"[{idx:5d}/{below}] pc{pc}={sorted(choice)}  Z-sum={ssum:.6g}{tag}", flush=True)

        # ---- single stage: cuts-only (cheapest sound INFEASIBLE filter) ----
        # cuts-only INFEASIBLE  => provably unreachable.
        # cuts-only UNKNOWN     => undecided; span/full/CMS tests run AFTER the
        #                          sweep on just these (and any SAT).
        # Passing the pi-free weight/stratum conditions subsumes the min-weight
        # cuts (same weight-level obstruction), so go straight to the exact model:
        # it is the only stage that can return SAT.
        st, pi = solve(S, "full", args.full_time)
        cache[key(S)] = {"status": st, "sumZ": ssum, "tier": "full",
                         "pi": pi if st == "sat" else None}
        save()
        if st == "infeasible":
            n_inf += 1
            print("        full INFEASIBLE", flush=True)
        elif st == "unknown":
            n_unk += 1
            print("        full UNKNOWN  (-> post-sweep escalation)", flush=True)
        else:
            ok = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp)) == S
            out = os.path.join(HERE, f"{log_prefix}_pi_{ssum:.6f}.matrix")
            with open(out, "w") as fh:
                fh.write(f"{n} {n}\n" + "\n".join(
                    "".join("1" if pi[a] == t else "0" for t in range(n))
                    for a in range(n)) + "\n")
            print(f"\n[sweep] *** SAT ***  {S}")
            print(f"        Z-sum {ssum:.6g} < equation {B_T:.6g}   verify={ok}")
            print(f"        pi -> {out}")
            print("[sweep] OPTIMAL reachable profile within the stratum family.")
            return 0

    print(f"\n[sweep] candidate list exhausted.  {n_inf} infeasible, {n_unk} undecided.")
    print("[sweep] => equation OPTIMAL within the stratum family"
          if n_unk == 0 else
          f"[sweep] => {n_unk} candidates need the post-sweep span/full/CMS tests "
          f"(status 'unknown' in {os.path.basename(cache_path)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
