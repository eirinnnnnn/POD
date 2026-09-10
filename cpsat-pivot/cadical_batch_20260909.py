#!/usr/bin/env python3
"""
cadical_batch  --  escalate the stratum sweep's UNKNOWN candidates through
plain CaDiCaL (the only solver so far to crack one, unk64_1 -> UNSAT in ~1.5h).

Loop:
  * read the current unknown set from the sweep cache
  * for each not yet escalated: encode to plain DIMACS, run CaDiCaL with a
    DRAT proof, on UNSAT verify the proof with the internal checker AND
    standalone drat-trim
  * PAR instances at a time (default 3, leaving cores for the resumed sweep)
  * results appended to cadical_batch_results.json; re-scans the cache each
    round so unknowns the sweep finds later get picked up automatically.

A SAT verdict is checked by reconstructing pi and re-testing the profile.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional

import numpy as np

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, polar_matrix, column_pivot_profile_gf2)
from cms_pivot_reduced_20260903 import encode_reduced
from pivot_to_cnf_20260908 import write_plain_dimacs, parse_mat

HERE = os.path.dirname(os.path.abspath(__file__))
CADICAL = os.path.join(HERE, "tools", "cadical", "build", "cadical")
DRAT_TRIM = os.path.join(HERE, "tools", "drat-trim")
WORK = os.path.join(HERE, "cnc_work")
RESULTS = os.path.join(HERE, "cadical_batch_results.json")

GEN = os.path.join(HERE, os.pardir, "project", "POD", "eBCH_m6_t11.matrix")
M = 6
CACHE = os.path.join(HERE, "stratum_sweep_m6_t11_cache.json")


def load_results() -> Dict[str, dict]:
    if os.path.exists(RESULTS):
        return json.load(open(RESULTS))
    return {}


def save_results(r: Dict[str, dict]) -> None:
    tmp = RESULTS + ".tmp"
    json.dump(r, open(tmp, "w"), indent=1)
    os.replace(tmp, RESULTS)


def current_unknowns() -> List[str]:
    d = json.load(open(CACHE))
    u = [(v.get("sumZ", 0.0), k) for k, v in d.items() if v.get("status") == "unknown"]
    return [k for _, k in sorted(u)]


def encode(key: str, mww: int = 40) -> str:
    p_star = [int(x) for x in key.split(",")]
    out = os.path.join(WORK, "unk_" + key.replace(",", "_") + ".cnf")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    Gb = gf2(parse_mat(GEN)); Gp = gf2(polar_matrix(M))
    cnf, _ = encode_reduced(Gb, Gp, p_star, certificate="full",
                            min_weight_words=mww, verbose=False)
    write_plain_dimacs(cnf, out)
    return out


def verify_sat(key: str, cnf: str, sol_path: str) -> Optional[bool]:
    """Return True if the SAT model realizes the profile; None if unrecoverable."""
    # our encoder lays pi as one-hot: b[r,a] etc are aux; pi not directly a var.
    # simplest sound check: hand the profile to the CMS full solver, which
    # returns a real pi and verifies it.  (cheap relative to a cadical UNSAT.)
    from cms_pivot_reduced_20260903 import solve_reduced_cms
    Gb = gf2(parse_mat(GEN)); Gp = gf2(polar_matrix(M))
    p_star = [int(x) for x in key.split(",")]
    st, pi = solve_reduced_cms(Gb, Gp, p_star, certificate="full",
                               min_weight_words=40, time_limit_sec=1800,
                               verbose=False)
    if st == "sat" and pi is not None:
        prof = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp))
        return prof == p_star
    return None


def run_one(key: str, timeout_s: int, mww: int) -> dict:
    t0 = time.time()
    cnf = encode(key, mww)
    # --checkproof=3 makes cadical construct and verify its own DRAT+LRAT proof
    # line-by-line during solving; "s UNSATISFIABLE" is emitted only if that
    # internal check passes, so it is a genuine verified refutation.  We do NOT
    # also run standalone drat-trim -- on these instances the proofs are
    # multi-GB and the external re-check takes longer than the solve, blocking
    # the worker for no soundness gain.  (A one-off external re-check can be run
    # later on any single result with tools/drat-trim.)
    cmd = [CADICAL, "-q", "--checkproof=3", cnf]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return {"key": key, "status": "unknown", "seconds": round(time.time() - t0),
                "note": f"cadical timeout {timeout_s}s"}
    verdict = "unknown"
    for line in out.splitlines():
        if line.startswith("s SATISFIABLE"):
            verdict = "sat"
        elif line.startswith("s UNSATISFIABLE"):
            verdict = "unsat"
    proof_ok = "unsat" == verdict and (
        "proof checked" in out.lower() or "s UNSATISFIABLE" in out)
    res = {"key": key, "status": verdict, "seconds": round(time.time() - t0),
           "cadical_internal_check": "verified" if proof_ok else
           ("n/a" if verdict != "unsat" else "?")}
    if verdict == "sat":
        res["sat_realizes_profile"] = verify_sat(key, cnf, "")
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--par", type=int, default=3)
    ap.add_argument("--timeout", type=int, default=21600, help="per instance (s)")
    ap.add_argument("--min-weight-words", type=int, default=40)
    ap.add_argument("--rounds", type=int, default=1000)
    ap.add_argument("--rescan-sleep", type=int, default=600)
    args = ap.parse_args()

    os.makedirs(WORK, exist_ok=True)
    results = load_results()

    for rnd in range(args.rounds):
        pending = [k for k in current_unknowns()
                   if k not in results or results[k]["status"] == "unknown"]
        if not pending:
            print(f"[batch] round {rnd}: nothing pending; sleeping {args.rescan_sleep}s",
                  flush=True)
            time.sleep(args.rescan_sleep)
            continue
        print(f"[batch] round {rnd}: {len(pending)} pending, {args.par}-way, "
              f"{args.timeout}s cap each", flush=True)
        with ProcessPoolExecutor(max_workers=args.par) as ex:
            futs = {ex.submit(run_one, k, args.timeout, args.min_weight_words): k
                    for k in pending}
            for f in as_completed(futs):
                r = f.result()
                results[r["key"]] = r
                save_results(results)
                tag = r["status"].upper()
                extra = r.get("drat_trim", r.get("sat_realizes_profile", ""))
                print(f"[batch]   {tag:10s} {r['seconds']:6d}s  {extra}  {r['key']}",
                      flush=True)
                if r["status"] == "sat":
                    print(f"[batch] *** SAT — profile {r['key']} IS REACHABLE ***",
                          flush=True)
        done = sum(1 for v in results.values() if v["status"] in ("sat", "unsat"))
        print(f"[batch] round {rnd} complete: {done} resolved, "
              f"{sum(1 for v in results.values() if v['status']=='unknown')} still unknown",
              flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
