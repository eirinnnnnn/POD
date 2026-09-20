#!/usr/bin/env python3
"""
Verify an ensemble file against the code it was generated for.

  algebraic : every permutation must satisfy rowspan(G[:, pi]) == rowspan(G)
  structural: report how many DISTINCT equivalence classes the ensemble covers
              (branches in the same EC give bit-identical SC output, so they
              are wasted slots)

Usage:
    python3 verify_aut.py --n 7 --G ../codes/<id>/<id>.matrix \
        --I_file ../codes/<id>/<id>_I.txt --aut ../codes/<id>/aut_utl_M32.txt
"""

import argparse

import numpy as np

from monomial_aut import (absorption_profile, admissible_matrix,
                          block_structure, in_blta, inv_mod2, load_I)
from polar_design import polar_matrix


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


def read_matrix(path):
    with open(path) as f:
        tok = f.read().split()
    r, c = int(tok[0]), int(tok[1])
    return np.array([int(x) for x in tok[2:2 + r * c]],
                    dtype=np.uint8).reshape(r, c)


def read_perms(path):
    with open(path) as f:
        tok = f.read().split()
    m, n = int(tok[0]), int(tok[1])
    return [[int(x) for x in tok[2 + i * n:2 + (i + 1) * n]] for i in range(m)]


def affine_of_perm(perm, n):
    """Recover (A, b) from a permutation, or None if it is not affine.
    b = pi(0); column i of A = bits(pi(2^i)) + b.
    """
    b = np.array([(perm[0] >> i) & 1 for i in range(n)], dtype=np.uint8)
    A = np.zeros((n, n), dtype=np.uint8)
    for i in range(n):
        col = np.array([(perm[1 << i] >> k) & 1 for k in range(n)],
                       dtype=np.uint8)
        A[:, i] = (col + b) % 2
    for j in range(1 << n):
        v = np.array([(j >> i) & 1 for i in range(n)], dtype=np.uint8)
        if int(sum(int(x) << i for i, x in enumerate((A @ v + b) % 2))) != perm[j]:
            return None
    return A, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--G", required=True)
    ap.add_argument("--I_file", required=True)
    ap.add_argument("--aut", required=True)
    args = ap.parse_args()

    n = args.n
    G = read_matrix(args.G)
    I = load_I(args.I_file)
    perms = read_perms(args.aut)
    rG = rank_mod2(G)

    S = block_structure(n, admissible_matrix(n, I))
    S1 = absorption_profile(n, S)

    bad = []
    non_affine = []
    reps = []
    ec_of = []
    for idx, p in enumerate(perms):
        if rank_mod2(np.vstack([G, G[:, p]])) != rG:
            bad.append(idx)
        ab = affine_of_perm(p, n)
        if ab is None:
            non_affine.append(idx)
            ec_of.append(None)
            continue
        A = ab[0]
        hit = None
        for k, A0 in enumerate(reps):
            if in_blta((A @ inv_mod2(A0)) % 2, S1):
                hit = k
                break
        if hit is None:
            reps.append(A)
            hit = len(reps) - 1
        ec_of.append(hit)

    print("[verify_aut] %s" % args.aut)
    print("  branches                : %d" % len(perms))
    print("  automorphism check      : %s"
          % ("ALL PASS" if not bad else "FAIL at %s" % bad))
    if non_affine:
        print("  non-affine permutations : %s" % non_affine)
    print("  distinct equivalence classes covered : %d / %d"
          % (len(reps), len(perms)))
    if len(reps) < len(perms):
        print("  -> %d branch(es) are SC-redundant (same EC as an earlier one)"
              % (len(perms) - len(reps)))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
