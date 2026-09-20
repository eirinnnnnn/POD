#!/usr/bin/env python3
"""
Polar code construction for PolarAED.

Builds an information set I from a reliability rule, then emits the two matrix
files the existing C++ decoder reads:

    <id>.matrix     K x N   generator, rows of F^{otimes n} indexed by I
    <id>_H.matrix   N x (N-K)   parity check, TRANSPOSED (repo convention,
                                same as generator_set_eBCH.py)

Index convention (see PLAN.md 1.2, verified by tools/_feasibility_check.py):
    F^{otimes n}[r][c] = 1  iff  c is a submask of r        (natural order)
    coordinate j  <->  (x_0, ..., x_{n-1}) with x_i = bit i of j   (LSB first)
    row r         <->  monomial containing x_i  iff  bit i of r is 0

Reliability rules:
    ga   Gaussian-approximation density evolution at a design SNR
    rm   Reed-Muller (by row weight / monomial degree), ties broken by GA

Usage:
    python3 polar_design.py --n 7 --K 100 --design ga --snr 10.5 \
        --out_dir ../codes/n128_k100_highsnr --id n128_k100_highsnr
"""

import argparse
import math
import os

import numpy as np


# ------------------------------------------------------------ GF(2) utilities

def polar_matrix(n):
    """F^{otimes n} in natural order: M[r][c] = 1 iff (c & r) == c."""
    N = 1 << n
    return np.array([[1 if (c & r) == c else 0 for c in range(N)]
                     for r in range(N)], dtype=np.uint8)


def nullspace_mod2(G):
    """Basis of {x : G x^T = 0}, returned as rows."""
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
    M = (M.astype(np.uint8) & 1)
    r, c = M.shape
    with open(path, "w") as f:
        f.write("%d %d\n" % (r, c))
        for i in range(r):
            f.write(" ".join(str(int(x)) for x in M[i, :]) + "\n")


# ------------------------------------------------- Gaussian approximation DE/GA

def _phi(x):
    """Chung et al. approximation of the GA phi function."""
    if x <= 0.0:
        return 1.0
    if x < 10.0:
        return math.exp(-0.4527 * (x ** 0.86) + 0.0218)
    return math.sqrt(math.pi / x) * math.exp(-x / 4.0) * (1.0 - 10.0 / (7.0 * x))


def _phi_inv(y):
    """Bisection inverse of _phi on (0, 1)."""
    if y >= 1.0:
        return 0.0
    if y <= 0.0:
        return 1e9
    lo, hi = 0.0, 1.0
    while _phi(hi) > y:
        hi *= 2.0
        if hi > 1e9:
            return 1e9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _phi(mid) > y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ga_means(n, design_snr_db, rate=None, snr_is_ebn0=True):
    """Mean LLRs of the N synthetic channels, natural (non-bit-reversed) order.

    Returned in the SAME index convention as polar_matrix: index N-1 is the
    most reliable (the weight-N row / constant monomial), index 0 the least.
    """
    N = 1 << n
    snr = 10.0 ** (design_snr_db / 10.0)
    if snr_is_ebn0:
        if rate is None:
            raise ValueError("Eb/N0 design SNR needs the code rate")
        esn0 = rate * snr
    else:
        esn0 = snr
    m = [4.0 * esn0] * N            # mean LLR of the physical channel

    # One polarization stage per bit of the index: bit b = 0 picks the "worse"
    # (minus) channel, b = 1 the "better" (plus) one, which is exactly the
    # submask ordering polar_matrix uses.
    #
    # Stage ORDER matters -- (W^-)^+ != (W^+)^-. For u.F^{otimes n} in natural
    # order the outermost butterfly (the one adjacent to the channel) is the
    # MSB and the innermost (adjacent to the u-domain decision) is the LSB, so
    # the transform sequence applied to W runs bit n-1 first down to bit 0 last.
    # This matches the decoder, where do_node_value() branches on bit `level`
    # and level 0 is the decision level. Getting this backwards silently swaps
    # the reliabilities of index pairs such as 1 and 2; see test_polar_design.py.
    for stage in reversed(range(n)):
        step = 1 << stage
        new = list(m)
        for i in range(N):
            if i & step:
                continue
            a, b = m[i], m[i | step]
            new[i] = _phi_inv(1.0 - (1.0 - _phi(a)) * (1.0 - _phi(b)))  # worse
            new[i | step] = a + b                                        # better
        m = new
    return m


def information_set(n, K, design, snr_db=None, rate=None, snr_is_ebn0=True):
    N = 1 << n
    means = ga_means(n, snr_db if snr_db is not None else 0.0,
                     rate=(rate if rate is not None else K / float(N)),
                     snr_is_ebn0=snr_is_ebn0)
    if design == "ga":
        key = [(means[i], i) for i in range(N)]
    elif design == "rm":
        # row weight 2^popcount(i); GA mean breaks ties inside a degree
        key = [(bin(i).count("1") * 1e9 + means[i], i) for i in range(N)]
    else:
        raise ValueError("unknown design %r" % design)
    order = [i for _, i in sorted(key, reverse=True)]   # most reliable first
    return sorted(order[:K]), order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, help="N = 2^n")
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--design", default="ga", choices=["ga", "rm"])
    ap.add_argument("--snr", type=float, default=None, help="design SNR in dB")
    ap.add_argument("--snr_type", default="ebn0", choices=["ebn0", "esn0"])
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--id", default=None)
    args = ap.parse_args()

    n, K, N = args.n, args.K, 1 << args.n
    code_id = args.id or "n%d_k%d" % (N, K)
    os.makedirs(args.out_dir, exist_ok=True)

    I, order = information_set(n, K, args.design, args.snr,
                               rate=K / float(N),
                               snr_is_ebn0=(args.snr_type == "ebn0"))
    G = polar_matrix(n)[I, :]
    H = nullspace_mod2(G)
    assert np.all((G @ H.T) % 2 == 0), "G H^T != 0"

    g_path = os.path.join(args.out_dir, code_id + ".matrix")
    h_path = os.path.join(args.out_dir, code_id + "_H.matrix")
    save_matrix(g_path, G)
    save_matrix(h_path, H.T)
    with open(os.path.join(args.out_dir, code_id + "_I.txt"), "w") as f:
        f.write(" ".join(str(i) for i in I) + "\n")

    dmin_bound = min(1 << bin(i).count("1") for i in I)
    print("[polar_design] (%d, %d) design=%s snr=%s (%s)"
          % (N, K, args.design, args.snr, args.snr_type))
    print("  G -> %s  %s" % (g_path, G.shape))
    print("  H -> %s  %s (transposed)" % (h_path, H.T.shape))
    print("  row-weight lower bound on d_min: %d" % dmin_bound)


if __name__ == "__main__":
    main()
