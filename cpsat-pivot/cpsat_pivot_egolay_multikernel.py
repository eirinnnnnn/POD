#!/usr/bin/env python3
"""
Trial: find a coordinate permutation pi for the extended Golay [24,12] code under
a MULTI-KERNEL polar transform G_p so that the dynamic-frozen construction places
the 12 information positions on the most reliable synthetic channels, i.e. the
pivot / info set minimises the Bhattacharyya sum.

Correspondence with the C++ pipeline (src/ErrorCorrectionCode/PolarMultiKernal.cpp
+ KernalManager.cpp + src/lib/libMath.cpp), verified numerically in the probe below
AND against a full Python replication of the C++ (scratchpad/replicate_cpp.py):

    C++ builds   constraint = G_p * P * H          (P[a, pi[a]] = 1, G_b H = 0,
                                                    H is 24 x 12)
    then         GaussianJordanElimination(constraint^T, corner=1)   -- RIGHT-to-
                 LEFT column sweep -- and marks position j "info" iff j is never
                 the highest-index nonzero of a reduced row.

    That is exactly the RIGHT-to-LEFT greedy independent-row set of `constraint`.
    Folding P and the index reversal J (a -> n-1-a) through the transpose gives

        info set  =  { n-1-q : q in [n] \\ PivotProfile( pi(H^T) * (G_p^T . J) ) }

    where pi(H^T)[r, j] = H^T[r, pi[j]] and (G_p^T . J) is G_p^T with its columns
    reversed.  (The greedy independent-row SET is invariant to the choice of H
    basis, so any nullspace basis of G_b works.)

Therefore the cpsat-pivot instance that targets "info set == p_star" is

        G_b'   = H^T                    (12 x 24; Golay is self-dual)
        G_p'   = G_p^T with columns reversed
        target = sorted{ n-1-i : i in [n] \\ p_star }
        solve  PivotProfile( pi(G_b') G_p' ) = target

and the returned list is exactly the permutation pi (write P[a, pi[a]] = 1).

G_p is not a standard Arikan tensor power, so the solver runs with the dense
W-encoding and no matroid-duality reduction (--dualize never).

Kernel / index conventions (replicating KernalManager::build_Gmatrix):
    kernal_string "657,23,23,23" -> stages [657,23,23,23], radices [3,2,2,2]
    digit s of index i  =  (i // prod(radices[:s])) % radices[s]   (stage 0 = LSB)
    G_p[i,j] = prod_s K_s[digit_s(i)][digit_s(j)]
    K_657 = [[1,1,0],[1,0,1],[1,1,1]]      K_23 = [[1,0],[1,1]]

Usage:
    python cpsat_pivot_egolay_multikernel.py --check-only        # probe + target
    python cpsat_pivot_egolay_multikernel.py --time 900 --workers 12
    python cpsat_pivot_egolay_multikernel.py --hint-current-perm
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List

import numpy as np

from ortools.sat.python import cp_model

from cpsat_pivot_reconstruct_duality import (
    gf2,
    gf2_matmul,
    gf2_rank,
    gf2_nullspace,
    column_pivot_profile_gf2,
    solve_reconstruction,
    build_pivot_reconstruction_model,
)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))

GOLAY_MATRIX = os.path.join(REPO, "data", "extend_Golay_12_24.matrix")
AED_DIR = os.path.join(REPO, "project", "Golay24", "AEDtest_20260108")
PERM_MATRIX = os.path.join(AED_DIR, "657_23_23_23_m1.matrix")
KERNAL_STRING = "657,23,23,23"

# Bhattacharyya values Z (probability domain, exp(-exp(.))) printed by the
# Golay24 binary for kernal_string 657,23,23,23 -- pasted from the user's run.
BHATTACHARYYA_Z: List[float] = [
    0.971283, 0.947432, 0.572901, 0.572477, 0.424491, 0.041475,
    0.423001, 0.284293, 0.013892, 0.032725, 0.016766, 0.000004,
    0.276202, 0.168185, 0.003324, 0.012017, 0.006063, 0.000000,
    0.006258, 0.003144, 0.000000, 0.000005, 0.000002, 0.000000,
]

KERNELS = {
    "23": np.array([[1, 0], [1, 1]], dtype=np.uint8),
    "753": np.array([[1, 1, 1], [1, 0, 1], [0, 1, 1]], dtype=np.uint8),
    "427": np.array([[1, 0, 0], [0, 1, 0], [1, 1, 1]], dtype=np.uint8),
    "657": np.array([[1, 1, 0], [1, 0, 1], [1, 1, 1]], dtype=np.uint8),
}


# ------------------------------------------------------------------
# instance construction
# ------------------------------------------------------------------

def parse_binary_matrix(path: str) -> np.ndarray:
    with open(path, "r") as fh:
        tokens = fh.read().split()
    rows, cols = int(tokens[0]), int(tokens[1])
    body = tokens[2:]
    if len(body) == rows:
        M = np.array([[int(c) for c in r] for r in body], dtype=np.uint8)
    else:
        M = np.array([int(x) for x in body], dtype=np.uint8).reshape(rows, cols)
    if M.shape != (rows, cols):
        raise ValueError(f"{path}: parsed shape {M.shape} != ({rows},{cols})")
    return M


def digit_array(index: int, radices: List[int]) -> List[int]:
    out = []
    for r in radices:
        out.append(index % r)
        index //= r
    return out


def build_multikernel_gp(kernal_string: str) -> np.ndarray:
    stages = kernal_string.split(",")
    radices = [len(s) for s in stages]
    n = int(np.prod(radices))
    Ks = [KERNELS[s] for s in stages]
    G = np.ones((n, n), dtype=np.uint8)
    for i in range(n):
        di = digit_array(i, radices)
        for j in range(n):
            dj = digit_array(j, radices)
            bit = 1
            for s in range(len(stages)):
                bit &= int(Ks[s][di[s], dj[s]])
            G[i, j] = bit
    return G


def perm_matrix_to_list(P: np.ndarray) -> List[int]:
    P = gf2(P)
    if not (P.sum(axis=0) == 1).all() or not (P.sum(axis=1) == 1).all():
        raise ValueError("permutation matrix is not a permutation")
    return [int(np.argmax(P[a])) for a in range(P.shape[0])]


def write_perm_matrix(path: str, pi: List[int]) -> None:
    n = len(pi)
    lines = [f"{n} {n}"]
    for a in range(n):
        row = ["0"] * n
        row[pi[a]] = "1"
        lines.append("".join(row))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def best_bhattacharyya_pivots(z: List[float], k: int) -> List[int]:
    order = sorted(range(len(z)), key=lambda i: (z[i], i))
    return sorted(order[:k])


def info_set_of_perm(HbT: np.ndarray, GpT_rev: np.ndarray, pi: List[int]) -> List[int]:
    """C++ dynamic-frozen info set:
        { n-1-q : q in [n] \\ PivotProfile( pi(H^T) (Gp^T . J) ) }.
    """
    n = HbT.shape[1]
    prof = set(column_pivot_profile_gf2(gf2_matmul(HbT[:, pi], GpT_rev)))
    return sorted(n - 1 - q for q in range(n) if q not in prof)


# ------------------------------------------------------------------
# solve for one exact target info set
# ------------------------------------------------------------------

SOLVER_BACKEND = "cpsat"          # "cpsat" | "reduced"  (set by --reduced)
MIN_WEIGHT_WORDS = 24


def solve_for_info_set(HbT, GpT_rev, want_info, *, time_limit_sec, workers,
                       log=False, hint_pi=None, verbose=False):
    """
    Try to realise C++ info_set(pi) == sorted(want_info).

    Returns (status, pi) with status in {"sat", "infeasible", "unknown"}.
      sat        -> pi is a permutation, verified by info_set_of_perm
      infeasible -> the solver PROVED no permutation realises this info set
      unknown    -> time limit hit before a decision (NOT a proof)
    """
    n = HbT.shape[1]
    want_info = sorted(want_info)
    frozen = sorted(set(range(n)) - set(want_info))
    target = sorted(n - 1 - i for i in frozen)

    if SOLVER_BACKEND == "reduced":
        from cpsat_pivot_reduced_20260902 import solve_reduced
        st, pi_sol = solve_reduced(
            HbT, GpT_rev, target,
            kernel_string=KERNAL_STRING, w_encoding="staged",
            certificate="full", min_weight_words=MIN_WEIGHT_WORDS,
            time_limit_sec=time_limit_sec, workers=workers, log=log, verbose=verbose,
        )
        if st == "sat" and info_set_of_perm(HbT, GpT_rev, pi_sol) != want_info:
            return ("unknown", None)
        return (st, pi_sol)

    model, variables = build_pivot_reconstruction_model(
        HbT, GpT_rev, target,
        w_encoding="dense", certificate="pivot-only", rank_cuts="span",
        row_weight_cuts=True, branch="popcount", verbose=verbose,
    )
    pi_vars = variables["pi"]
    if hint_pi is not None:
        for a in range(n):
            model.AddHint(pi_vars[a], int(hint_pi[a]))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_sec)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.log_search_progress = bool(log)
    st = solver.Solve(model)

    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        pi_sol = [solver.Value(pi_vars[a]) for a in range(n)]
        if info_set_of_perm(HbT, GpT_rev, pi_sol) != want_info:
            return ("unknown", None)   # should not happen; treat as undecided
        return ("sat", pi_sol)
    if st == cp_model.INFEASIBLE:
        return ("infeasible", None)
    return ("unknown", None)


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", type=float, default=1800.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--log", action="store_true")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--search", action="store_true",
                    help="fix the --keep-best most reliable indices, enumerate the "
                         "completions to k positions in ascending Bhattacharyya sum, "
                         "and test each until a reachable info set is found")
    ap.add_argument("--keep-best", type=int, default=10,
                    help="number of most-reliable indices to force into the info set "
                         "(--search mode)")
    ap.add_argument("--max-trials", type=int, default=40,
                    help="max completions to test in --search mode")
    ap.add_argument("--full-search", action="store_true",
                    help="DEFINITIVE: enumerate every k-subset info set with "
                         "Bhattacharyya sum below --beat (ascending) and test each. "
                         "First SAT = global optimum; exhausting the list proves the "
                         "current permutation optimal among all reachable info sets.")
    ap.add_argument("--beat", type=float, default=None,
                    help="Bhattacharyya-sum threshold for --full-search "
                         "(default: the current 657_23_23_23_m1 info-set sum)")
    ap.add_argument("--resume-from", type=int, default=0,
                    help="(deprecated; --full-search now uses fullsearch_cache.json)")
    ap.add_argument("--retry-unknown", action="store_true",
                    help="in --full-search, re-test candidates previously left UNKNOWN")
    ap.add_argument("--reduced", action="store_true",
                    help="use cpsat_pivot_reduced_20260902 (staged W + min-weight cuts) "
                         "as the solver backend")
    ap.add_argument("--hint-current-perm", action="store_true",
                    help="feed the existing 657_23_23_23_m1 permutation as a CP-SAT hint")
    ap.add_argument("--npz-out", default=os.path.join(HERE, "egolay_multikernel_instance.npz"))
    ap.add_argument("--perm-out", default=os.path.join(HERE, "egolay_multikernel_pi.matrix"))
    args = ap.parse_args()

    global SOLVER_BACKEND
    if args.reduced:
        SOLVER_BACKEND = "reduced"

    Gb = parse_binary_matrix(GOLAY_MATRIX)
    k, n = Gb.shape
    if (k, n) != (12, 24):
        raise SystemExit(f"unexpected Golay generator shape {Gb.shape}")

    Gp = build_multikernel_gp(KERNAL_STRING)
    HbT = gf2(gf2_nullspace(Gb))                  # 12 x 24, rows span C_b^perp (= "H^T")
    Jrev = np.eye(n, dtype=np.uint8)[::-1].copy()
    GpT_rev = gf2_matmul(Gp.T, Jrev)             # Gp^T with columns reversed
    if gf2_rank(HbT) != k or gf2_rank(GpT_rev) != n:
        raise SystemExit("rank check failed on H^T / (Gp^T . J)")

    p_star = best_bhattacharyya_pivots(BHATTACHARYYA_Z, k)          # target info set
    frozen = sorted(set(range(n)) - set(p_star))
    target = sorted(n - 1 - i for i in frozen)                     # cpsat pivot target
    print(f"n={n} k={k}   G_p multikernel '{KERNAL_STRING}'  (rank {gf2_rank(Gp)})")
    print(f"target info set p_star (min Bhattacharyya sum) = {p_star}")
    print(f"  sum Z on p_star = {sum(BHATTACHARYYA_Z[i] for i in p_star):.6g}")
    print(f"frozen = [n] \\ p_star                          = {frozen}")
    print(f"cpsat pivot target = sorted(n-1-i for i in frozen) = {target}")

    # ---- probe: reproduce the C++ info set for the existing permutation ----
    pi_cur = perm_matrix_to_list(parse_binary_matrix(PERM_MATRIX))
    info_cur = info_set_of_perm(HbT, GpT_rev, pi_cur)
    print(f"\n[probe] existing 657_23_23_23_m1 permutation:")
    print(f"        info set = {info_cur}")
    print(f"        sum Z    = {sum(BHATTACHARYYA_Z[i] for i in info_cur):.6g}"
          f"   ({'already optimal' if info_cur == p_star else 'not optimal'})")

    if args.check_only:
        return 0

    # ---- full-search: definitive optimality test ------------------------
    if args.full_search:
        import itertools
        import json
        z = BHATTACHARYYA_Z
        beat = args.beat if args.beat is not None else sum(z[i] for i in info_cur)
        cache_path = os.path.join(HERE, "fullsearch_cache.json")

        # sound candidate pool: index i can lie in a k-subset with sum < beat
        # iff z[i] + (sum of the k-1 smallest z over j != i) < beat.
        pool = []
        for i in range(n):
            others = sorted(z[j] for j in range(n) if j != i)
            if z[i] + sum(others[:k - 1]) < beat - 1e-12:
                pool.append(i)
        pool.sort()
        cand = []
        for S in itertools.combinations(pool, k):
            ssum = sum(z[i] for i in S)
            if ssum < beat - 1e-12:
                cand.append((ssum, sorted(S)))
        cand.sort()

        cache = {}
        if os.path.exists(cache_path):
            with open(cache_path) as fh:
                cache = json.load(fh)

        def key(S):
            return ",".join(map(str, S))

        def save():
            with open(cache_path, "w") as fh:
                json.dump(cache, fh, indent=0)

        n_inf = sum(1 for v in cache.values() if v["status"] == "infeasible")
        n_unk = sum(1 for v in cache.values() if v["status"] == "unknown")
        print(f"\n[full-search] threshold beat = {beat:.6g}")
        print(f"[full-search] candidate pool ({len(pool)}): {pool}")
        print(f"[full-search] {len(cand)} info sets; cache has "
              f"{n_inf} infeasible, {n_unk} unknown, "
              f"{len(cache)-n_inf-n_unk} sat.  time limit {args.time:.0f}s/instance\n",
              flush=True)

        for idx, (ssum, S) in enumerate(cand, 1):
            prev = cache.get(key(S))
            if prev and prev["status"] in ("infeasible", "sat"):
                continue                       # already decided
            if prev and prev["status"] == "unknown" and not args.retry_unknown:
                continue
            print(f"[{idx:4d}/{len(cand)}] {S}  sumZ={ssum:.6g}", flush=True)
            status, pi_sol = solve_for_info_set(
                HbT, GpT_rev, S, time_limit_sec=args.time, workers=args.workers,
                log=args.log,
            )
            cache[key(S)] = {"status": status, "sumZ": ssum,
                             "pi": pi_sol if status == "sat" else None}
            save()
            print(f"        {status.upper()}", flush=True)
            if status == "sat":
                out = args.perm_out.replace(".matrix", f"_opt_{ssum:.6f}.matrix")
                write_perm_matrix(out, pi_sol)
                print(f"        pi = {pi_sol}")
                print(f"        wrote {out}")
                print(f"\n[full-search] GLOBAL OPTIMUM among reachable info sets "
                      f"with sumZ < {beat:.6g}:  sumZ = {ssum:.6g}  (current "
                      f"{sum(z[i] for i in info_cur):.6g})")
                return 0

        # ---- verdict ----
        undecided = [(v["sumZ"], k2) for k2, v in cache.items()
                     if v["status"] == "unknown"]
        decided_inf = [k2 for k2, v in cache.items() if v["status"] == "infeasible"]
        tested = {key(S) for _, S in cand}
        missing = [S for _, S in cand if key(S) not in cache]
        if undecided or missing:
            print(f"\n[full-search] NOT CONCLUSIVE.")
            print(f"  {len(decided_inf)} candidates proved infeasible.")
            if missing:
                print(f"  {len(missing)} candidates not yet tested.")
            if undecided:
                print(f"  {len(undecided)} candidates UNKNOWN (time limit hit, "
                      f"no proof). Re-run with a larger --time and --retry-unknown:")
                for sm, k2 in sorted(undecided):
                    print(f"      sumZ={sm:.6g}  [{k2}]")
            return 1

        print(f"\n[full-search] EXHAUSTED: all {len(cand)} candidate info sets with "
              f"Bhattacharyya sum below {beat:.6g} are provably INFEASIBLE.")
        print(f"[full-search] => the 657_23_23_23_m1 permutation "
              f"(sumZ {sum(z[i] for i in info_cur):.6g}) is OPTIMAL among all "
              f"coordinate permutations.")
        return 0

    # ---- search mode: fix the best indices, vary the rest ----------------
    if args.search:
        import itertools
        z = BHATTACHARYYA_Z
        by_rel = sorted(range(n), key=lambda i: (z[i], i))
        keep = sorted(by_rel[:args.keep_best])
        pool = sorted(by_rel[args.keep_best:])          # candidate extra positions
        need = k - len(keep)
        cur_sum = sum(z[i] for i in info_cur)
        print(f"\n[search] forcing info-set superset {keep}")
        print(f"[search] choosing {need} more from {pool}")
        print(f"[search] beating current sum Z = {cur_sum:.6g}\n")

        combos = []
        for extra in itertools.combinations(pool, need):
            S = sorted(keep + list(extra))
            combos.append((sum(z[i] for i in S), S))
        combos.sort()

        found = []
        for rank_i, (ssum, S) in enumerate(combos[:args.max_trials], 1):
            tag = "" if ssum < cur_sum else "  (>= current, stopping)"
            print(f"[{rank_i:2d}] try info set {S}  sum Z = {ssum:.6g}{tag}")
            if ssum >= cur_sum:
                break
            status, pi_sol = solve_for_info_set(
                HbT, GpT_rev, S, time_limit_sec=args.time, workers=args.workers,
                log=args.log, hint_pi=(pi_cur if args.hint_current_perm else None),
            )
            if status != "sat":
                print(f"     -> {status.upper()}"
                      + ("" if status == "infeasible" else "  (time limit, no proof)"))
                continue
            print(f"     -> SAT   pi = {pi_sol}")
            found.append((ssum, S, pi_sol))
            out = args.perm_out.replace(".matrix", f"_search_{ssum:.6f}.matrix")
            write_perm_matrix(out, pi_sol)
            print(f"     -> wrote {out}")
            break   # first SAT is the best by construction (ascending sum)

        if not found:
            print("\n[search] no reachable info set beats the current permutation "
                  "within the tested completions.")
            return 1
        ssum, S, pi_sol = found[0]
        print(f"\n[search] BEST reachable: info set {S}")
        print(f"         sum Z = {ssum:.6g}   (current {cur_sum:.6g})")
        return 0

    np.savez(args.npz_out, Gb=HbT, Gp=GpT_rev,
             p_star=np.array(target, dtype=np.int64))
    print(f"\ninstance (Gb=H^T, Gp=Gp^T.J, p_star=target) -> {args.npz_out}")

    hint = pi_cur if args.hint_current_perm else None
    result = solve_reconstruction(
        HbT, GpT_rev, target,
        time_limit_sec=args.time,
        workers=args.workers,
        log=args.log,
        w_encoding="dense",
        certificate="pivot-only",
        dualize="never",
        hint_pi=hint,
        verbose=True,
    )

    if result is None:
        print("\nNo permutation found: UNSAT (p_star unreachable by any coordinate "
              "permutation) or time-out.")
        return 1

    pi_sol, _ = result
    info_sol = info_set_of_perm(HbT, GpT_rev, pi_sol)
    print("\nRecovered pi:", pi_sol)
    print("achieved info set =", info_sol)
    ok = info_sol == p_star
    print("info set == p_star :", "PASS" if ok else "FAIL")
    print(f"sum Z on achieved info set = {sum(BHATTACHARYYA_Z[i] for i in info_sol):.6g}")

    if ok:
        write_perm_matrix(args.perm_out, pi_sol)
        print(f"permutation matrix (KernalManager format) -> {args.perm_out}")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
