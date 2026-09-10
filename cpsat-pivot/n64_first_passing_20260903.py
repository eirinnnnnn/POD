#!/usr/bin/env python3
"""
eBCH [64,16] under standard polar F^{x6}.

Walk the k=16 pivot profiles in DESCENDING 5G-reliability-sum order
(best Bhattacharyya sum first).  For each, run a cheap PREFILTER:

    certificate='span-only' + min-weight cuts, 60 s cap.

span-only INFEASIBLE  ==>  provably unreachable, skip.
Otherwise the profile "passes the prefilter": launch the FULL SAT
(certificate='full' + min-weight cuts, no real time cap) and TIME it
until it returns SAT or INFEASIBLE.
"""
from __future__ import annotations
import heapq, sys, time
import numpy as np

sys.path.insert(0, "../project/PtP")
from static_ptp_5GNR import NR_5G_POLAR_SEQUENCE_1024
from cpsat_pivot_reconstruct_duality import polar_matrix, gf2, gf2_rank, column_pivot_profile_gf2, gf2_matmul
from cpsat_pivot_reduced_20260902 import solve_reduced

N, K = 64, 16


def parse_mat(p):
    t = open(p).read().split(); r, c = int(t[0]), int(t[1]); body = t[2:]
    if len(body) == r:
        M = np.array([[int(x) for x in row] for row in body], np.uint8)
    else:
        M = np.array([int(x) for x in body], np.uint8).reshape(r, c)
    return M


def main():
    Gb = gf2(parse_mat("../data/BCH_16_64.matrix"))
    Gp = polar_matrix(6)
    assert Gb.shape == (K, N) and gf2_rank(Gb) == K

    order = [i for i in NR_5G_POLAR_SEQUENCE_1024 if i < N]   # least -> most reliable
    rank = {ch: r for r, ch in enumerate(order)}              # 0..63, higher = better
    members = set(order[-K:])
    nonmembers = set(order[:-K])

    # best-first over k-subsets by total rank (descending), 1-swap neighbours
    start = tuple(sorted(members))
    seen = {start}
    heap = [(-sum(rank[c] for c in start), start)]
    tried = 0

    while heap:
        negsum, S = heapq.heappop(heap)
        Sset = set(S)
        tried += 1
        rsum = -negsum
        p_star = list(S)
        print(f"\n[{tried}] profile {p_star}  rank-sum={rsum}", flush=True)

        st, pi = solve_reduced(Gb, Gp, p_star, w_encoding="staged",
                               certificate="span-only", min_weight_words=48,
                               time_limit_sec=60, workers=8, verbose=False)
        if st == "infeasible":
            print("     prefilter: INFEASIBLE (span-only) -> skip", flush=True)
        else:
            print(f"     prefilter: {st.upper()} -> PASSES.  launching full SAT...", flush=True)
            t0 = time.time()
            st2, pi2 = solve_reduced(Gb, Gp, p_star, w_encoding="staged",
                                     certificate="full", min_weight_words=48,
                                     time_limit_sec=36000, workers=8, verbose=True)
            dt = time.time() - t0
            print(f"\n=== FULL SAT on first prefilter-passing profile ===")
            print(f"    profile      : {p_star}")
            print(f"    status       : {st2}")
            print(f"    wall time    : {dt:.2f} s")
            if st2 == "sat":
                actual = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi2], Gp))
                print(f"    profile check: {'PASS' if actual == p_star else 'FAIL'}")
                print(f"    pi           : {pi2}")
            return

        # expand neighbours: swap one member out, one non-member in
        for out in Sset:
            for inn in nonmembers | (members - Sset):
                if inn in Sset:
                    continue
                T = tuple(sorted((Sset - {out}) | {inn}))
                if T not in seen:
                    seen.add(T)
                    heapq.heappush(heap, (-sum(rank[c] for c in T), T))

        if tried >= 400:
            print("\n[stop] 400 profiles all failed the prefilter.")
            return


if __name__ == "__main__":
    main()
