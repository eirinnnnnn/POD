#!/usr/bin/env python3
"""
Companion to `cpsat_pivot_egolay_multikernel.py --full-search`.

Runs ALONGSIDE the main sweep: repeatedly reads its result cache
(fullsearch_cache.json), picks up every info set the main run left
UNKNOWN (time-limit, no proof), and re-attacks it with a stronger
budget -- larger time limit, more minimum-weight cuts, and a
CryptoMiniSat fallback for the UNSAT direction.

Writes its own cache (fullsearch_unknowns.json) so it never races the
main run's file.  Stops when the main run has finished (its process is
gone) and no unknowns remain.

    python capture_unknowns_20260902.py --time 3600 --min-weight-words 48
"""
from __future__ import annotations

import argparse
import json
import os
import time
from typing import List

import numpy as np

import cpsat_pivot_egolay_multikernel as M
from cpsat_pivot_reconstruct_duality import gf2, gf2_matmul, gf2_nullspace
from cpsat_pivot_reduced_20260902 import solve_reduced

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN_CACHE = os.path.join(HERE, "fullsearch_cache.json")
MY_CACHE = os.path.join(HERE, "fullsearch_unknowns.json")


def load(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def cms_fallback(HbT, GpT_rev, target, time_limit, threads):
    """CryptoMiniSat on the same instance (dense).  Returns 'sat'|'infeasible'|'unknown'."""
    try:
        from cpsat_pivot_crypto_minixor import solve_reconstruction_cms
    except Exception as e:
        print("  (cms unavailable:", e, ")")
        return "unknown", None
    st, pis = solve_reconstruction_cms(
        HbT, GpT_rev, target, time_limit_sec=time_limit, threads=threads,
        w_encoding="dense", dualize="never", prepass=True, prepass_tries=50000,
        verbose=False,
    )
    if st == "UNSAT":
        return "infeasible", None
    if st == "SAT" and pis:
        return "sat", pis[0]
    return "unknown", None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", type=float, default=3600.0)
    ap.add_argument("--workers", type=int, default=4,
                    help="keep < nproc so the main sweep keeps running")
    ap.add_argument("--min-weight-words", type=int, default=48)
    ap.add_argument("--poll", type=float, default=30.0)
    ap.add_argument("--no-cms", action="store_true")
    args = ap.parse_args()

    Gb = M.parse_binary_matrix(M.GOLAY_MATRIX)
    Gp = M.build_multikernel_gp(M.KERNAL_STRING)
    n = 24
    HbT = gf2(gf2_nullspace(Gb))
    J = np.eye(n, dtype=np.uint8)[::-1].copy()
    GpT_rev = gf2_matmul(Gp.T, J)
    Z = M.BHATTACHARYYA_Z
    cur_sum = sum(Z[i] for i in M.info_set_of_perm(HbT, GpT_rev,
                  M.perm_matrix_to_list(M.parse_binary_matrix(M.PERM_MATRIX))))
    print(f"[capture] current m1 sum Z = {cur_sum:.6g}; hunting reachable info sets below it")

    mine = load(MY_CACHE)
    while True:
        main_cache = load(MAIN_CACHE)
        pending = [k for k, v in main_cache.items()
                   if v.get("status") == "unknown" and k not in mine]
        pending.sort(key=lambda k: main_cache[k]["sumZ"])

        if not pending:
            # done?  main run's rc line present and no fresh unknowns
            done = os.path.exists(os.path.join(HERE, "fullsearch_reduced.log")) and any(
                "EXHAUSTED" in L or "GLOBAL OPTIMUM" in L or "NOT CONCLUSIVE" in L
                for L in open(os.path.join(HERE, "fullsearch_reduced.log")))
            if done:
                break
            time.sleep(args.poll)
            continue

        for key in pending:
            S = [int(x) for x in key.split(",")]
            ssum = main_cache[key]["sumZ"]
            frozen = sorted(set(range(n)) - set(S))
            target = sorted(n - 1 - i for i in frozen)
            print(f"\n[capture] {key}  sumZ={ssum:.6g}", flush=True)

            t0 = time.time()
            st, pi = solve_reduced(
                HbT, GpT_rev, target, kernel_string=M.KERNAL_STRING,
                w_encoding="staged", certificate="full",
                min_weight_words=args.min_weight_words,
                time_limit_sec=args.time, workers=args.workers, verbose=False,
            )
            src = "reduced"
            if st == "unknown" and not args.no_cms:
                print(f"  reduced UNKNOWN in {time.time()-t0:.0f}s; trying CryptoMiniSat", flush=True)
                st, pi = cms_fallback(HbT, GpT_rev, target, args.time, args.workers)
                src = "cms"

            if st == "sat" and M.info_set_of_perm(HbT, GpT_rev, pi) != S:
                st, pi = "unknown", None
            mine[key] = {"status": st, "sumZ": ssum, "pi": pi, "src": src,
                         "secs": round(time.time() - t0, 1)}
            with open(MY_CACHE, "w") as fh:
                json.dump(mine, fh, indent=0)
            print(f"  -> {st.upper()} ({src}, {mine[key]['secs']}s)", flush=True)

            if st == "sat":
                out = os.path.join(HERE, f"egolay_multikernel_pi_capture_{ssum:.6f}.matrix")
                M.write_perm_matrix(out, pi)
                print(f"  *** REACHABLE info set BELOW m1: {S}  sumZ={ssum:.6g}")
                print(f"  *** pi = {pi}")
                print(f"  *** wrote {out}")

    sats = [(v["sumZ"], k, v) for k, v in mine.items() if v["status"] == "sat"]
    unk = [k for k, v in mine.items() if v["status"] == "unknown"]
    print("\n[capture] finished.")
    print(f"  infeasible: {sum(1 for v in mine.values() if v['status']=='infeasible')}")
    print(f"  still unknown: {len(unk)}  {unk}")
    if sats:
        best = min(sats)
        print(f"  BEST reachable below m1: sumZ={best[0]:.6g}  set {best[1]}")
    else:
        print("  no reachable info set below m1 found among the recaptured unknowns.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
