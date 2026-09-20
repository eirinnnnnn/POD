#!/usr/bin/env python3
"""
Reproduces the four feasibility checks recorded in PLAN.md §1.4.

These pin down the conventions the whole PolarAED extension rests on:
  1. a plain polar code drives AdjustPolarDecoder to a purely STATIC frozen set
     (no dynamic-frozen rows), so the existing C++ decoder needs no change;
  2. the admissibility test reproduces [P21] §III-A's worked (16,7) example;
  3. the affine -> permutation direction is  pi(j) = idx(A . bits(j) + b)
     with LSB-first bit order (not A^T, not A^-1);
  4. an LTA automorphism is absorbed by SC bit-for-bit, a UTL one is not.

Checks 1 and 4 shell out to a prebuilt POD binary. Point --binary at one
(e.g. project/POD/m7t10_itw2026/POD) or pass --no-sim to skip them.

Usage:
    python3 _feasibility_check.py --workdir /tmp/polaraed_feas \
        --binary ../../POD/m7t10_itw2026/POD
"""

import argparse
import os
import re
import subprocess
import sys

import numpy as np


# ---------------------------------------------------------------- GF(2) utils

def polar_matrix(n):
    """F^{otimes n} in natural order: M[r][c] = 1 iff c is a submask of r.

    This is exactly what AdjustPolarDecoder builds from an all-ones
    operationArray, and it matches the monomial convention of [P21, Table I]:
    row r <-> monomial containing x_i iff bit i of r is 0.
    """
    N = 1 << n
    return np.array([[1 if (c & r) == c else 0 for c in range(N)]
                     for r in range(N)], dtype=np.uint8)


def rank_mod2(A):
    A = (A.astype(np.uint8) & 1).copy()
    m, nn = A.shape
    r = 0
    for c in range(nn):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(r + 1, m):
            if A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
    return r


def nullspace_mod2(G):
    A = (G.astype(np.uint8) & 1).copy()
    m, nn = A.shape
    r = 0
    pivots = []
    for c in range(nn):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        pivots.append(c)
        r += 1
    pivot_set = set(pivots)
    free = [j for j in range(nn) if j not in pivot_set]
    p2r = {c: i for i, c in enumerate(pivots)}
    rows = []
    for f in free:
        v = np.zeros(nn, dtype=np.uint8)
        v[f] = 1
        for p in pivots:
            v[p] = A[p2r[p], f]
        rows.append(v)
    return np.vstack(rows) if rows else np.zeros((0, nn), dtype=np.uint8)


def save_matrix(path, M):
    """Repo matrix format: '<rows> <cols>' header, then space-separated rows."""
    M = (M.astype(np.uint8) & 1)
    r, c = M.shape
    with open(path, "w") as f:
        f.write("%d %d\n" % (r, c))
        for i in range(r):
            f.write(" ".join(str(int(x)) for x in M[i, :]) + "\n")


def save_perms(path, perms, N):
    with open(path, "w") as f:
        f.write("%d %d\n" % (len(perms), N))
        for p in perms:
            f.write(" ".join(str(int(x)) for x in p) + "\n")


# ------------------------------------------------------- monomial automorphism

def admissible(n, I, i, j):
    """[P21, Thm. 2] in binary-expansion form: may A[i][j] be set to 1?

    Substituting x_i <- x_j removes x_i (set bit i) and adds x_j (clear bit j).
    Admissible iff every affected information index stays in I.
    """
    Iset = set(I)
    for r in I:
        if (r >> i) & 1:      # bit i == 1 means x_i is absent from m_r
            continue
        if ((r | (1 << i)) & ~(1 << j)) not in Iset:
            return False
    return True


def admissible_matrix(n, I):
    return [[admissible(n, I, i, j) if i != j else True for j in range(n)]
            for i in range(n)]


def bits(j, n):
    return np.array([(j >> i) & 1 for i in range(n)], dtype=np.uint8)


def index_of(v, n):
    return int(sum(int(v[i]) << i for i in range(n)))


def perm_from_affine(A, n, b=None):
    """pi(j) = idx(A . bits(j) + b), bits LSB-first. See PLAN.md 1.2."""
    N = 1 << n
    b = np.zeros(n, dtype=np.uint8) if b is None else b
    return [index_of((A @ bits(j, n) + b) % 2, n) for j in range(N)]


def is_automorphism(G, perm):
    return rank_mod2(np.vstack([G, G[:, perm]])) == rank_mod2(G)


# -------------------------------------------------------------------- ini glue

INI = """[AWGN]
start = {snr}
end = {snr}
step = 0.5
step_type = SNR
seed_string = -2

[AdjustPolarDecoder]
matrix_src = byHmatrix
Gmatrix_path =
Hmatrix_path = {h}
operationArray = {op}
bha_value_setting = {q}
target_raw_BER = 0.01
list_size = 1
permutation_src =
permutation_random_seed = -1
OnlyInit = {only_init}

use_AED = {use_aed}
automorphism_src = {aut}
aed_L = {aed_L}

[Monte_Carlo]
iter_min = 0
iter_max = {iters}
error_min = 0
error_max = 100000000
monitor_slot_size = 100000
seed_string = 111511015
"""


def write_ini(path, n, N, h_path, snr, only_init=False, aut="", aed_L=0,
              iters=20000):
    with open(path, "w") as f:
        f.write(INI.format(
            snr=snr, h=h_path, op="1" * (n * (N // 2)), q="?" * N,
            only_init="true" if only_init else "false",
            use_aed="true" if aut else "false",
            aut=aut, aed_L=aed_L, iters=iters))


def run(binary, ini, cwd):
    return subprocess.run([os.path.abspath(binary), "-ini", ini],
                          cwd=cwd, capture_output=True, text=True).stdout


# ----------------------------------------------------------------- the checks

def check1_static_frozen(work, binary):
    """[P21] (16,7) example -> relation string must contain no '2'."""
    n, N = 4, 16
    I = [7, 10, 11, 12, 13, 14, 15]
    G = polar_matrix(n)[I, :]
    H = nullspace_mod2(G)
    save_matrix(os.path.join(work, "p16_7.matrix"), G)
    save_matrix(os.path.join(work, "p16_7_H.matrix"), H.T)
    if binary is None:
        return None
    write_ini(os.path.join(work, "c1.ini"), n, N, "./p16_7_H.matrix", 2.0,
              only_init=True)
    out = run(binary, "./c1.ini", work)
    m = re.search(r"total_bhattacharyya_value\s+[-\d.eE+]+\s+([012]+)", out)
    rel = m.group(1) if m else ""
    info = [i for i, ch in enumerate(rel) if ch == "0"]
    ok = bool(rel) and "2" not in rel and info == I
    print("  relation string : %s" % rel)
    print("  info positions  : %s   (expected %s)" % (info, I))
    print("  no dynamic-frozen rows: %s" % ("2" not in rel if rel else "?"))
    return ok


def check2_admissibility():
    """UT admissible must be exactly {(1,2)}; all LT positions admissible."""
    n = 4
    I = [7, 10, 11, 12, 13, 14, 15]
    ut = [(i, j) for i in range(n) for j in range(i + 1, n)
          if admissible(n, I, i, j)]
    lt = [(i, j) for i in range(n) for j in range(i) if admissible(n, I, i, j)]
    print("  UT admissible : %s   (expected [(1, 2)])" % ut)
    print("  LT admissible : %d of %d (expected all -> decreasing, LTA <= Aut)"
          % (len(lt), n * (n - 1) // 2))
    return ut == [(1, 2)] and len(lt) == n * (n - 1) // 2


def check3_direction():
    """A gives automorphisms where A^T does not, for a non-symmetric case."""
    n = 4
    I = [5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
    G = polar_matrix(n)[I, :]
    discriminating = 0
    for i in range(n):
        for j in range(i):
            if not admissible(n, I, i, j):
                continue
            A = np.eye(n, dtype=np.uint8)
            A[i, j] = 1
            ok_a = is_automorphism(G, perm_from_affine(A, n))
            ok_t = is_automorphism(G, perm_from_affine(A.T.copy(), n))
            if ok_a and not ok_t:
                discriminating += 1
            if ok_t and not ok_a:
                print("  !! transposed convention would be the right one")
                return False
    print("  discriminating lower entries (A aut, A^T not): %d" % discriminating)
    return discriminating > 0


def check4_absorption(work, binary):
    """LTA branch must equal plain SC bit-for-bit; UTL branch must not."""
    n, N = 4, 16
    I = [5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
    G = polar_matrix(n)[I, :]
    H = nullspace_mod2(G)
    save_matrix(os.path.join(work, "p16_10.matrix"), G)
    save_matrix(os.path.join(work, "p16_10_H.matrix"), H.T)

    L = np.eye(n, dtype=np.uint8); L[3, 0] = 1; L[2, 1] = 1     # lower -> LTA
    U = np.eye(n, dtype=np.uint8); U[0, 1] = 1; U[2, 3] = 1     # upper -> UTL
    for A, name in ((L, "lta"), (U, "utl")):
        perm = perm_from_affine(A, n)
        assert is_automorphism(G, perm), "%s element is not an automorphism" % name
        save_perms(os.path.join(work, "aut_%s.txt" % name), [perm], N)

    if binary is None:
        return None

    bler = {}
    for tag, aut, L_ in (("sc", "", 0), ("lta", "./aut_lta.txt", 1),
                         ("utl", "./aut_utl.txt", 1)):
        ini = "c4_%s.ini" % tag
        write_ini(os.path.join(work, ini), n, N, "./p16_10_H.matrix", 1.0,
                  aut=aut, aed_L=L_, iters=20000)
        for stale in ("log.txt",):
            p = os.path.join(work, ini[:-4], stale)
            if os.path.exists(p):
                os.remove(p)
        out = run(binary, "./" + ini, work)
        hits = re.findall(r"BLER =\s+(\d+)/\s*(\d+) = ([\d.]+)", out)
        bler[tag] = hits[-1] if hits else None
        print("  %-3s : %s" % (tag, bler[tag]))

    if not all(bler.values()):
        return None
    absorbed = bler["lta"] == bler["sc"]
    distinct = bler["utl"] != bler["sc"]
    print("  LTA absorbed by SC : %s   (must be True)" % absorbed)
    print("  UTL not absorbed   : %s   (must be True)" % distinct)
    return absorbed and distinct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default="./_feas")
    ap.add_argument("--binary", default=None,
                    help="prebuilt POD/PolarAED binary; omit or use --no-sim "
                         "to run the algebraic checks only")
    ap.add_argument("--no-sim", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    binary = None if args.no_sim else args.binary

    results = {}
    print("[1] plain polar code -> static frozen set")
    results["1"] = check1_static_frozen(args.workdir, binary)
    print("[2] admissibility test vs [P21] (16,7) example")
    results["2"] = check2_admissibility()
    print("[3] affine -> permutation direction")
    results["3"] = check3_direction()
    print("[4] SC absorption of LTA, not of UTL")
    results["4"] = check4_absorption(args.workdir, binary)

    print("\nsummary:")
    failed = False
    for k in sorted(results):
        v = results[k]
        print("  check %s : %s" % (k, "SKIP" if v is None else
                                   ("PASS" if v else "FAIL")))
        failed |= (v is False)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
