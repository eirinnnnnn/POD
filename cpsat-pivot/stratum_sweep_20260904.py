#!/usr/bin/env python3
"""
Stratum-preserving reachability sweep  --  eBCH [64,16] / F^{x6}.

See STRATUM_SWEEP_PLAN.md.  Walk the C(15,9)=5005 profiles that share the
equation set's weight-stratum profile, in ASCENDING Bhattacharyya sum, and
for each:

  1. span-only + min-weight cuts,  short cap   (span-only UNSAT => unreachable)
  2. full certificate + min-weight cuts, generous cap
  3. UNKNOWN -> escalate to the CryptoMiniSat-reduced encoding

Stop at the first full SAT: by ascending order it is provably the optimal
reachable profile within the stratum-preserving family.

    python stratum_sweep_20260904.py --snr 3.0 --workers 8 \
        --span-time 60 --full-time 900 --min-weight-words 40
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import time
from typing import Dict, List, Tuple

import numpy as np

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, gf2_rank, polar_matrix, column_pivot_profile_gf2,
)
from cpsat_pivot_reduced_20260902 import solve_reduced, low_weight_codewords

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, os.pardir, "project", "POD", "eBCH_m6_t11.matrix")
CACHE = os.path.join(HERE, "stratum_sweep_cache.json")
N, K, M = 64, 16, 6

W5PLUS = [31, 47, 55, 59, 61, 62, 63]                                   # forced
W4 = [15, 23, 27, 29, 30, 39, 43, 45, 46, 51, 53, 54, 57, 58, 60]      # choose 9
EQ_W4 = {15, 23, 27, 29, 30, 39, 43, 46, 51}


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


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snr", type=float, default=3.0, help="design Eb/N0 (dB)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--span-time", type=float, default=60.0)
    ap.add_argument("--full-time", type=float, default=900.0)
    ap.add_argument("--min-weight-words", type=int, default=40)
    ap.add_argument("--max-cand", type=int, default=5005)
    args = ap.parse_args()

    Gb = gf2(parse_mat(GEN))
    Gp = gf2(polar_matrix(M))
    assert Gb.shape == (K, N) and gf2_rank(Gb) == K
    h = K - gf2_rank(gf2_matmul(Gb, Gb.T))

    Z = bha_seq(M, args.snr, K / N)
    B_T = sum(Z[i] for i in W5PLUS) + sum(Z[i] for i in EQ_W4)

    # enumerate the family, score, sort ascending
    cands = []
    for choice in itertools.combinations(W4, 9):
        S = sorted(W5PLUS + list(choice))
        cands.append((sum(Z[i] for i in S), S))
    cands.sort()
    print(f"[sweep] eBCH[64,16] self-orthogonal (hull dim h={h}); design SNR {args.snr} dB")
    print(f"[sweep] B_T (equation Z-sum) = {B_T:.6g}")
    print(f"[sweep] {sum(1 for s,_ in cands if s < B_T-1e-15)} / {len(cands)} candidates below B_T")
    print(f"[sweep] beta Z-sum = {cands[0][0]:.6g}  ({cands[0][1]})", flush=True)

    cache = {}
    if os.path.exists(CACHE):
        cache = json.load(open(CACHE))

    def key(S):
        return ",".join(map(str, S))

    def save():
        json.dump(cache, open(CACHE, "w"), indent=0)

    def solve(S, cert, tlim, mww):
        return solve_reduced(Gb, Gp, S, w_encoding="staged", certificate=cert,
                             min_weight_words=mww, time_limit_sec=tlim,
                             workers=args.workers, verbose=False)

    n_inf = n_unk = 0
    for idx, (ssum, S) in enumerate(cands[:args.max_cand], 1):
        if ssum >= B_T - 1e-15:
            print(f"\n[sweep] reached B_T at candidate {idx}; every cheaper profile "
                  f"decided.  results: {n_inf} infeasible, {n_unk} unknown.")
            print("[sweep] => the equation set is OPTIMAL within the stratum family"
                  if n_unk == 0 else
                  "[sweep] => NOT conclusive: unknowns remain below B_T (see cache)")
            return 0

        prev = cache.get(key(S))
        if prev and prev["status"] in ("infeasible", "sat"):
            if prev["status"] == "sat":
                print(f"\n[sweep] cached SAT at {S} (Z-sum {ssum:.6g}) -> optimal.")
                return 0
            continue

        chosen4 = sorted(set(S) - set(W5PLUS))
        tag = "  [=equation]" if set(chosen4) == EQ_W4 else ""
        print(f"\n[{idx:4d}/{len(cands)}] w4={chosen4}  Z-sum={ssum:.6g}{tag}", flush=True)

        st, pi = solve(S, "span-only", args.span_time, args.min_weight_words)
        if st == "infeasible":
            cache[key(S)] = {"status": "infeasible", "sumZ": ssum, "tier": "span"}
            save(); n_inf += 1
            print(f"        span-only INFEASIBLE ({args.span_time:.0f}s cap)", flush=True)
            continue

        st, pi = solve(S, "full", args.full_time, args.min_weight_words)
        cache[key(S)] = {"status": st, "sumZ": ssum, "tier": "full",
                         "pi": pi if st == "sat" else None}
        save()
        print(f"        full -> {st.upper()}", flush=True)
        if st == "infeasible":
            n_inf += 1
        elif st == "unknown":
            n_unk += 1
            print("        (escalate later with cms_pivot_reduced_20260903)", flush=True)
        elif st == "sat":
            ok = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp)) == S
            out = os.path.join(HERE, f"stratum_sweep_pi_{ssum:.6f}.matrix")
            with open(out, "w") as fh:
                fh.write(f"{N} {N}\n" + "\n".join(
                    "".join("1" if pi[a] == t else "0" for t in range(N))
                    for a in range(N)) + "\n")
            print(f"\n[sweep] *** SAT ***  profile {S}")
            print(f"        Z-sum {ssum:.6g}  (< equation {B_T:.6g})   verify={ok}")
            print(f"        pi -> {out}")
            print(f"[sweep] this is the OPTIMAL reachable profile in the stratum family"
                  f"{'  == beta' if set(chosen4) == set(x for _,s in [cands[0]] for x in s if bin(x).count('1')==4) else ''}")
            return 0

    print("\n[sweep] candidate list exhausted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
