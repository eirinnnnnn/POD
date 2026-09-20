#!/usr/bin/env python3
"""
Sample an automorphism ensemble for a polar code and write it in the format
AdjustPolarDecoder::load_automorphism_set() reads:

    <M> <N>
    <permutation of 0..N-1>
    ...

Row 0 is always the identity, so an ensemble file with M rows gives M AE
branches including the un-permuted one.

Modes
  utl     random UTL matrices supported on the admissible upper positions [P21]
  lta     random LTA elements (lower-triangular + translation) -- the negative
          control: every one of these is absorbed by SC
  pu      EC representatives A = P.U, P in A_P, U in A_U            [P22 Thm. 5]
  random  uniform from BLTA(S)                                      [P22 baseline]

Except in "lta" mode, candidates are de-duplicated by EQUIVALENCE CLASS using
the exact test A1 . A2^{-1} in [1]  ([P22, Lem. 6]) -- two automorphisms in the
same EC give bit-identical SC output and would waste a branch.

Usage:
    python3 blta_sample.py --I_file ../codes/<id>/<id>_I.txt --n 7 \
        --mode utl --M 32 --seed 1 --out ../codes/<id>/aut_utl_M32.txt
"""

import argparse
import itertools
import random

import numpy as np

from monomial_aut import (admissible_matrix, absorption_profile,
                          block_structure, inv_mod2, load_I, num_ec, same_ec)


# ------------------------------------------------------- affine -> permutation

def bits(j, n):
    return np.array([(j >> i) & 1 for i in range(n)], dtype=np.uint8)


def index_of(v, n):
    return int(sum(int(v[i]) << i for i in range(n)))


def perm_from_affine(A, n, b=None):
    """pi(j) = idx(A . bits(j) + b), LSB-first. See PLAN.md 1.2."""
    N = 1 << n
    b = np.zeros(n, dtype=np.uint8) if b is None else b
    return [index_of((A @ bits(j, n) + b) % 2, n) for j in range(N)]


# --------------------------------------------------------------- generators

def eye(n):
    return np.eye(n, dtype=np.uint8)


def block_ranges(S):
    out = []
    off = 0
    for s in S:
        out.append((off, off + s))
        off += s
    return out


def rand_utl(n, ut_positions, rng):
    """Upper-triangular unit-diagonal matrix supported on ut_positions."""
    A = eye(n)
    for (i, j) in ut_positions:
        if rng.random() < 0.5:
            A[i, j] = 1
    return A


def rand_lta(n, rng):
    A = eye(n)
    for i in range(n):
        for j in range(i):
            A[i, j] = rng.randint(0, 1)
    b = np.array([rng.randint(0, 1) for _ in range(n)], dtype=np.uint8)
    return A, b


def rand_gl(s, rng):
    """Uniform-ish invertible s x s over F2 by rejection."""
    while True:
        M = np.array([[rng.randint(0, 1) for _ in range(s)] for _ in range(s)],
                     dtype=np.uint8)
        try:
            inv_mod2(M)
            return M
        except ValueError:
            continue


def rand_blta(n, S, rng):
    """Uniform from BLTA(S), linear part only (b is absorbed anyway)."""
    A = np.zeros((n, n), dtype=np.uint8)
    rngs = block_ranges(S)
    for (lo, hi) in rngs:
        A[lo:hi, lo:hi] = rand_gl(hi - lo, rng)
    for bi, (lo, hi) in enumerate(rngs):
        for (lo2, hi2) in rngs[:bi]:
            for i in range(lo, hi):
                for j in range(lo2, hi2):
                    A[i, j] = rng.randint(0, 1)
    return A


def rand_pu(n, S, rng):
    """A = P . U with P in A_P (block-wise permutation), U in A_U (block-wise
    strict upper triangular). [P22, Thm. 5]: every EC has such a representative.
    Also returns the (p, v) descriptor vectors used by the distance heuristic.
    """
    P = np.zeros((n, n), dtype=np.uint8)
    U = eye(n)
    p_desc = []
    v_desc = []
    for (lo, hi) in block_ranges(S):
        s = hi - lo
        perm = list(range(s))
        rng.shuffle(perm)
        for c in range(s):
            P[lo + perm[c], lo + c] = 1
        p_desc += [lo + x for x in perm]
        for i in range(lo, hi):
            for j in range(i + 1, hi):
                bit = rng.randint(0, 1)
                U[i, j] = bit
                v_desc.append(bit)
    return (P @ U) % 2, np.array(p_desc), np.array(v_desc)


def hamming(a, b):
    return int(np.sum(np.asarray(a) != np.asarray(b)))


# ------------------------------------------------------------------- sampling

def sample(n, I, mode, M, seed, D=None, max_tries=200000):
    adm = admissible_matrix(n, I)
    S = block_structure(n, adm)
    S1 = absorption_profile(n, S)
    ut = [(i, j) for i in range(n) for j in range(i + 1, n) if adm[i][j]]
    rng = random.Random(seed)

    accepted_A = [eye(n)]                       # identity is branch 0
    accepted_desc = [(None, None)]
    perms = [list(range(1 << n))]
    seen = {tuple(perms[0])}

    tries = 0
    while len(perms) < M and tries < max_tries:
        tries += 1
        desc = (None, None)
        b = None
        if mode == "utl":
            if not ut:
                raise SystemExit("no admissible UT positions: UTL ensemble is "
                                 "trivial, this code cannot be used for AE-SC")
            A = rand_utl(n, ut, rng)
        elif mode == "lta":
            A, b = rand_lta(n, rng)
        elif mode == "random":
            A = rand_blta(n, S, rng)
        elif mode == "pu":
            A, p_desc, v_desc = rand_pu(n, S, rng)
            desc = (p_desc, v_desc)
        else:
            raise ValueError(mode)

        if mode != "lta":
            # exact EC de-duplication
            if any(same_ec(A, A0, n, S1) for A0 in accepted_A):
                continue
            # optional Hamming-distance heuristic on the (p, v) descriptors
            if D is not None and mode == "pu":
                dP, dU = D
                ok = True
                for (p0, v0) in accepted_desc:
                    if p0 is None:
                        continue
                    if hamming(desc[0], p0) < dP or hamming(desc[1], v0) < dU:
                        ok = False
                        break
                if not ok:
                    continue

        perm = perm_from_affine(A, n, b)
        key = tuple(perm)
        if key in seen:
            continue
        seen.add(key)
        perms.append(perm)
        accepted_A.append(A)
        accepted_desc.append(desc)

    return perms, S, S1, ut, tries


def write_perms(path, perms, N):
    with open(path, "w") as f:
        f.write("%d %d\n" % (len(perms), N))
        for p in perms:
            f.write(" ".join(str(int(x)) for x in p) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I_file", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--mode", required=True,
                    choices=["utl", "lta", "pu", "random"])
    ap.add_argument("--M", type=int, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dP", type=int, default=None)
    ap.add_argument("--dU", type=int, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    n = args.n
    I = load_I(args.I_file)
    D = None if args.dP is None else (args.dP, args.dU or 0)
    perms, S, S1, ut, tries = sample(n, I, args.mode, args.M, args.seed, D)

    write_perms(args.out, perms, 1 << n)
    print("[blta_sample] mode=%s S=%s [1]=%s |EC|=%d admissible-UT=%d"
          % (args.mode, tuple(S), tuple(S1), num_ec(n, S), len(ut)))
    print("  requested M=%d, produced %d branches (identity included) "
          "after %d draws" % (args.M, len(perms), tries))
    if len(perms) < args.M:
        print("  NOTE: ensemble exhausted -- the code has fewer distinct "
              "classes than requested")
    print("  -> %s" % args.out)


if __name__ == "__main__":
    main()
