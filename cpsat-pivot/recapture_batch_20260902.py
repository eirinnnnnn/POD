#!/usr/bin/env python3
"""
Fast consolidated re-test of every info set the --full-search left UNKNOWN.

Two-tier per candidate (both sound for the UNSAT direction):
  1. certificate='span-only'  (no k^3 U layer, tiny model) + min-weight cuts.
     span-only INFEASIBLE  ==>  the info set is provably unreachable.
  2. if span-only is SAT/unknown, escalate to certificate='full' with a
     large budget (a SAT here would be a real reachable info set below m1).

Reads/merges  fullsearch_cache.json  and  fullsearch_unknowns.json,
writes results back into  fullsearch_cache.json  (status upgraded from
'unknown' to 'infeasible' / 'sat').

    python recapture_batch_20260902.py --span-time 120 --full-time 1800 --workers 8
"""
from __future__ import annotations
import argparse, json, os, time
import numpy as np

import cpsat_pivot_egolay_multikernel as M
from cpsat_pivot_reconstruct_duality import gf2, gf2_matmul, gf2_nullspace
from cpsat_pivot_reduced_20260902 import solve_reduced

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "fullsearch_cache.json")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--span-time", type=float, default=120.0)
    ap.add_argument("--full-time", type=float, default=1800.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--min-weight-words", type=int, default=48)
    args = ap.parse_args()

    Gb = M.parse_binary_matrix(M.GOLAY_MATRIX)
    Gp = M.build_multikernel_gp(M.KERNAL_STRING)
    n = 24
    HbT = gf2(gf2_nullspace(Gb))
    J = np.eye(n, dtype=np.uint8)[::-1].copy()
    GpT_rev = gf2_matmul(Gp.T, J)
    Z = M.BHATTACHARYYA_Z
    m1sum = sum(Z[i] for i in M.info_set_of_perm(
        HbT, GpT_rev, M.perm_matrix_to_list(M.parse_binary_matrix(M.PERM_MATRIX))))

    with open(CACHE) as fh:
        cache = json.load(fh)
    # merge in the companion capture's results if present
    comp = os.path.join(HERE, "fullsearch_unknowns.json")
    if os.path.exists(comp):
        for k, v in json.load(open(comp)).items():
            if v["status"] in ("infeasible", "sat"):
                cache[k] = v

    todo = sorted((v["sumZ"], k) for k, v in cache.items() if v["status"] == "unknown")
    print(f"[recapture] m1 sumZ={m1sum:.6g}   {len(todo)} unknown info sets to settle")

    def solve(S, cert, tlim):
        frozen = sorted(set(range(n)) - set(S))
        target = sorted(n - 1 - i for i in frozen)
        return solve_reduced(HbT, GpT_rev, target, kernel_string=M.KERNAL_STRING,
                             w_encoding="staged", certificate=cert,
                             min_weight_words=args.min_weight_words,
                             time_limit_sec=tlim, workers=args.workers, verbose=False)

    n_inf = n_sat = n_unk = 0
    for idx, (ssum, key) in enumerate(todo, 1):
        S = [int(x) for x in key.split(",")]
        t0 = time.time()
        st, pi = solve(S, "span-only", args.span_time)
        tier = "span"
        if st != "infeasible":
            st2, pi2 = solve(S, "full", args.full_time)
            tier = "full"
            st, pi = st2, pi2
        if st == "sat" and M.info_set_of_perm(HbT, GpT_rev, pi) != sorted(S):
            st, pi = "unknown", None

        cache[key] = {"status": st, "sumZ": ssum, "pi": pi, "tier": tier,
                      "secs": round(time.time() - t0, 1)}
        with open(CACHE, "w") as fh:
            json.dump(cache, fh, indent=0)

        if st == "infeasible":
            n_inf += 1
        elif st == "sat":
            n_sat += 1
            out = os.path.join(HERE, f"egolay_multikernel_pi_capture_{ssum:.6f}.matrix")
            M.write_perm_matrix(out, pi)
            print(f"[{idx}/{len(todo)}] *** SAT  {S}  sumZ={ssum:.6g} < m1  -> {out}")
        else:
            n_unk += 1
        print(f"[{idx}/{len(todo)}] {key}  sumZ={ssum:.6g}  {st.upper()} "
              f"({tier}, {cache[key]['secs']}s)   running: {n_inf} inf / {n_sat} sat / {n_unk} unk",
              flush=True)

    print(f"\n[recapture] done.  {n_inf} infeasible, {n_sat} sat, {n_unk} still unknown")
    if n_sat == 0 and n_unk == 0:
        print("[recapture] EVERY info set with sumZ < m1 is provably infeasible => m1 OPTIMAL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
