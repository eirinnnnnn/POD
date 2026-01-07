#!/usr/bin/env python3
import os
import json
import random
from pathlib import Path
import numpy as np
import math
import subprocess
import re
from datetime import datetime

# Use your existing group code
from schreier_sims import *
# ( 
    # load_generators_from_json,
    # schreier_sims,
    # compose,
    # inverse,
    # right_coset_key,
# 

# ===================== I/O helpers =====================

def write_perm_matrix(path: Path, perm: np.ndarray):
    """
    Writes an N x N permutation matrix with header 'N N'
    SA_v2-compatible format (no spaces inside rows).
    """
    N = len(perm)
    with path.open("w") as f:
        f.write(f"{N} {N}\n")
        for i in range(N):
            row = ["0"] * N
            row[int(perm[i])] = "1"
            f.write("".join(row) + "\n")

def render_ini(Gmatrix_path, perm_matrix_path, operationArray, bha_value_setting):
    return f"""[AWGN]
end                            = 20.0
seed_string                    = -2
start                          = 2.0
step                           = 0.5
step_type                      = SNR

[AdjustPolarDecoder]
Gmatrix_path                   = {Gmatrix_path}
Hmatrix_path                   =
OnlyInit                       = true
bha_value_setting              = {bha_value_setting}
list_size                      = 4
matrix_src                     = byGmatrix
operationArray                 = {operationArray}
permutation_random_seed        = -1
permutation_src                = {perm_matrix_path}
target_raw_BER                 = 0.01

[Monte_Carlo]
error_max                      = 1
error_min                      = 1
iter_max                       = 1
iter_min                       = 0
monitor_slot_size              = 1
seed_string                    = -1
"""

# ===================== Group logic =====================

def random_perm_np(N, rng: random.Random):
    a = list(range(N))
    rng.shuffle(a)
    return np.array(a, dtype=int)

def canonical_right_coset_rep(P: np.ndarray, chain, n: int):
    """
    Compute canonical representative of right coset P H
    using Schreier–Sims stabilizer chain.

    Returns:
        rep_perm : tuple (length n), canonical permutation
    """

    # identity permutation
    id_perm = np.arange(n, dtype=int)

    base = [level["base_point"] for level in chain]

    best_tuple = None
    best_perm = None

    def dfs(level, h, prefix):
        nonlocal best_tuple, best_perm

        if level == len(chain):
            Ph = compose(P, h)
            tup = tuple(int(Ph[b]) for b in base)
            if best_tuple is None or tup < best_tuple:
                best_tuple = tup
                best_perm = tuple(int(x) for x in Ph)
            return

        level_data = chain[level]
        b = level_data["base_point"]
        orbit = level_data["orbit"]
        trans = level_data["trans"]

        # possible images of b under P h t
        # sort by value to enable pruning
        Ph = compose(P, h)
        candidates = sorted(orbit, key=lambda x: Ph[x])

        # candidates = sorted(orbit, key=lambda x: P[x])

        for x in candidates:
            # val = int(P[x])
            val = int(Ph[x])

            # prefix pruning
            if best_tuple is not None:
                if prefix + (val,) >= best_tuple[: level + 1]:
                    continue

            t = trans[x]          # t(b) = x
            h_next = compose(h, t)
            dfs(level + 1, h_next, prefix + (val,))

    dfs(0, id_perm, ())
    return best_perm

# def canon_right_coset(P: np.ndarray, chain, n: int) -> tuple:
#     """
#     Canonical representative of right coset PH.
#     Safe version matching Schreier–Sims orbit structure.
#     """
#     g = P.copy()

#     for level in chain:
#         b = level["base_point"]
#         orbit = level["orbit"]
#         trans = level["trans"]

#         x = int(g[b])

#         # if x is not in the orbit, we cannot strip further
#         if x not in orbit:
#             break

#         # identity transversal
#         if x == b:
#             continue

#         # safe: x is in orbit and not base
#         t = trans[x]
#         g = compose(g, inverse(t))

#     return tuple(int(v) for v in g)

def random_h(chain, n, rng):
    # make a random element of H by multiplying random strong generators from chain levels
    h = np.arange(n)
    for level in chain:
        gens = level["gens"]
        if not gens:
            continue
        g = rng.choice(gens)
        if rng.random() < 0.5:
            h = compose(h, g)
    return h


# ============== sanity check =================
def sanity_check_A(chain, N, trials=20, seed=0):
    rng = random.Random(seed)

    for i in range(trials):
        P = random_perm_np(N, rng)
        h = random_h(chain, N, rng)

        rep1 = canonical_right_coset_rep(P, chain, N)
        rep2 = canonical_right_coset_rep(compose(P, h), chain, N)

        if rep1 != rep2:
            print("[FAIL A] Coset invariance broken")
            print("P       =", P)
            print("h       =", h)
            print("rep(P)  =", rep1)
            print("rep(Ph) =", rep2)
            return

    print("[PASS A] Coset invariance verified")

def sanity_check_B(chain, N, trials=20, seed=1):
    rng = random.Random(seed)

    for i in range(trials):
        P = random_perm_np(N, rng)
        Q = random_perm_np(N, rng)

        same_coset = sift_membership(
            compose(inverse(Q), P), chain, N
        )

        rep_equal = (
            canonical_right_coset_rep(P, chain, N)
            == canonical_right_coset_rep(Q, chain, N)
        )

        if same_coset != rep_equal:
            print("[FAIL B] Membership disagreement")
            print("P =", P)
            print("Q =", Q)
            print("Q^{-1}P in H =", same_coset)
            print("rep(P) == rep(Q) =", rep_equal)
            return

    print("[PASS B] Membership equivalence verified")

# ===================== Main =====================

def main():

    m = 6
    t = 11
    # ========= config =========
    codetype = f"eBCH_m{m}_t{t}"
    aut_json = f"aut_generators_eBCH_m{m}_t{t}.json"
    aut_json = f"generator_set_eBCH_m{m}_t{t}.json"
    N = 2**m

    max_trials = 10000000
    seed = 111511015 

    # paths
    work_dir = Path(f"./aut_aware_timed_{m}_{t}")
    work_dir.mkdir(exist_ok=True)

    perm_path = work_dir / "tmp.matrix"
    ini_path  = work_dir / "tmp.ini"

    Gmatrix_path_in_ini = f"./eBCH_m{m}_t{t}.matrix"

    m = int(math.log2(N))
    operationArray = "1" * ((N // 2) * m)
    bha_value_setting = "?" * N

    rng = random.Random(seed)


    # Build Schreier–Sims chain (with transversals in RAM)
    n_json, gens, meta = load_generators_from_json(aut_json)
    assert n_json == N
    base = list(range(N))
    chain = schreier_sims(gens, base, N)

    # testing canon check
    sanity_check_A(chain, N)
    sanity_check_B(chain, N)

    # Quick invariance test for right cosets: key(P) == key(P ∘ h)
    # for _ in range(200):
    #     P = random_perm_np(N, rng)
    #     h = gens[rng.randrange(len(gens))]   # pick a random generator in H
    #     if right_coset_key(P, chain, N) != right_coset_key(compose(P, h), chain, N):
    #         print("NOT invariant (should not happen)")
    #         print(P)
    #         print(h)
    #         print(right_coset_key(P, chain, N))
    #         print(right_coset_key(compose(P, h), chain, N))
    #         break
    # else:
    #     print("[ok] invariant under right-mult by H (generators)")



    # Write ini ONCE
    ini_path.write_text(render_ini(
        Gmatrix_path=Gmatrix_path_in_ini,
        perm_matrix_path=str(perm_path),
        operationArray=operationArray,
        bha_value_setting=bha_value_setting,
    ))

    seen = set()          # canonical coset keys
    trials = 0
    accepted = 0

    print("[info] starting coset-distinct random search")

    import re

    re_total = re.compile(r"total_bhattacharyya_value\s+([0-9]*\.[0-9]+)")

    best = {
        "total_bha": float("inf"),
        "perm": None,
        "trial": None,
        "accepted": None,
    }

    best_path = work_dir / "best.json"
    best_log_path = work_dir / "best.log"

    while trials < max_trials:
        trials += 1
        P = random_perm_np(N, rng)

        # key = canon_right_coset(P, chain, N)
        # key = right_coset_key(P, chain, N)

        # # ====== REPRESENTATIVE TRAVERSAL ==========
        key = canonical_right_coset_rep(P, chain, N)
        if key in seen:
            continue

        seen.add(key)
        # # ====== REPRESENTATIVE TRAVERSAL ==========

        accepted += 1
        # overwrite permutation matrix
        write_perm_matrix(perm_path, P)

        # ===== evaluate (placeholder) =====
        out = subprocess.run(
            ["./m7t4_path_metric_AED/AED_20260103", "-ini", str(ini_path)],
            text=True, capture_output=True, check=True
        ).stdout

        m = re_total.search(out)
        if not m:
            raise RuntimeError("Cannot find total_bhattacharyya_value in output")
        total_bha = float(m.group(1))

        if total_bha < best["total_bha"]:
            best["total_bha"] = total_bha
            best["perm"] = P.tolist()
            best["trial"] = trials
            best["accepted"] = accepted

            best_path.write_text(json.dumps(best, indent=2))
            best_log_path.write_text(out)

            print(f"[best] total_bha={total_bha} trials={trials} accepted={accepted}")

    print("[done]")
    print(f"  trials   = {trials}")
    print(f"  accepted = {accepted}")
    print(f"  unique cosets = {len(seen)}")

if __name__ == "__main__":
    main()
