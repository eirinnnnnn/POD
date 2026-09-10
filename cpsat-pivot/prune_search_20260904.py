#!/usr/bin/env python3
"""
prune_search  --  single-kernel-toggle (single-layer) pruning search.

For a FIXED base code C_b and a FIXED coordinate permutation P, enumerate the
single-kernel prunings

    R_1 = { R_{t,q} : t in [1..m], q in [0..N/2) },     |R_1| = m*N/2

and, for each, evaluate the induced dynamic-frozen structure and objective.

Setup / conventions
-------------------
Index a in [0,N) has bits a_1..a_m, MSB first:  a_t = (a >> (m-t)) & 1.

Layer operator (row-vector convention  x -> x L_t):

    L_t  pairs the indices differing in bit a_t.
    On a pair (a0, a1) with a1 = a0 + 2^(m-t):
         (x[a0], x[a1])  ->  (x[a0] + x[a1], x[a1])          [kernel F]
    Pruning that kernel replaces F by I_2, i.e. deletes the single matrix
    entry L_t[a1, a0].  Hence  L_t^pruned = L_t + e_{a1} e_{a0}^T  (rank one).

Kernel index q at layer t = the m-1 remaining bits of a with bit a_t deleted.

Transform (layer 1 adjacent to the channel, layer m adjacent to u):

    T = L_m L_{m-1} ... L_1          (unpruned:  T = F^{(x m)}, verified)

Message-space representation of the fixed code  C = rowspace(G_b P), where P is
the coordinate permutation as an OPERATOR (G_b P = G_b[:, pi]).  The .matrix
files on disk store the transpose (P[a, pi(a)] = 1), so they are transposed on
load; the RM(2,5) regression test below pins this down.

    M = G_b P T^{-1},   RREF(M) -> pivot set J = information set,
    non-pivot columns of the RREF give the dynamic-frozen relations.

Bhattacharyya under pruning
---------------------------
Z propagates from the channel inward, layer 1 first.  On an unpruned kernel
(general two-channel form, exact on the BEC and the standard upper bound
elsewhere):

    Z[a0] <- Z[a0] + Z[a1] - Z[a0] Z[a1]        (degraded)
    Z[a1] <- Z[a0] Z[a1]                        (upgraded)

On a PRUNED kernel the two channels never combine, so both Z values pass
through unchanged.  With no prunings this reproduces the standard sequence
    z = [v for x in z for v in (2x - x^2, x^2)]   (regression-tested).

Objective
---------
    Phi(T, D) = sum_{i in J(D)} Z_i(T)                       [basic]
    Phi_lambda = Phi - lambda * sum_{f in F(D) cap [m1,m2]} Delta_f(T, D)

Delta_f is the dynamic-frozen informativeness: for a dynamic-frozen position f
whose relation involves r earlier information positions, Delta_f is taken as
the reliability mass that relation makes available, sum over those positions of
(1 - Z).  r = 0 (static frozen) contributes nothing.
"""
from __future__ import annotations

import argparse
import math
import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from cpsat_pivot_reconstruct_duality import gf2, gf2_matmul, gf2_rref, polar_matrix

Prune = Tuple[int, int]  # (t, q),  t in 1..m,  q in [0, N/2)


# ------------------------------------------------------------------
# geometry
# ------------------------------------------------------------------
def pair_of(m: int, t: int, q: int) -> Tuple[int, int]:
    """Indices (a0, a1) of kernel q at layer t.  a1 = a0 + 2^(m-t)."""
    s = m - t
    a0 = ((q >> s) << (s + 1)) | (q & ((1 << s) - 1))
    return a0, a0 | (1 << s)


def all_prunes(m: int) -> List[Prune]:
    return [(t, q) for t in range(1, m + 1) for q in range(1 << (m - 1))]


def layer_matrix(m: int, t: int, pruned: Optional[Set[int]] = None) -> np.ndarray:
    """L_t as an N x N matrix over F_2 (row-vector convention)."""
    n = 1 << m
    L = np.eye(n, dtype=np.uint8)
    for q in range(n >> 1):
        if pruned and q in pruned:
            continue
        a0, a1 = pair_of(m, t, q)
        L[a1, a0] = 1
    return L


def pruned_transform(m: int, prunes: Sequence[Prune] = ()) -> np.ndarray:
    """T = L_m ... L_1 with the given kernels pruned."""
    by_layer: Dict[int, Set[int]] = {}
    for t, q in prunes:
        by_layer.setdefault(t, set()).add(q)
    n = 1 << m
    T = np.eye(n, dtype=np.uint8)
    for t in range(m, 0, -1):                      # T = L_m @ ... @ L_1
        T = gf2_matmul(T, layer_matrix(m, t, by_layer.get(t)))
    return T


# ------------------------------------------------------------------
# reliability
# ------------------------------------------------------------------
def z_pruned(m: int, prunes: Sequence[Prune], z0: float) -> List[float]:
    """Bhattacharyya parameters of the N synthetic channels under pruning."""
    by_layer: Dict[int, Set[int]] = {}
    for t, q in prunes:
        by_layer.setdefault(t, set()).add(q)
    n = 1 << m
    z = [z0] * n
    for t in range(1, m + 1):                      # channel side inward
        skip = by_layer.get(t, ())
        for q in range(n >> 1):
            if q in skip:
                continue
            a0, a1 = pair_of(m, t, q)
            x, y = z[a0], z[a1]
            z[a0] = x + y - x * y
            z[a1] = x * y
    return z


def z0_awgn(ebn0_db: float, rate: float) -> float:
    return math.exp(-rate * 10 ** (ebn0_db / 10))


# ------------------------------------------------------------------
# dynamic-frozen structure
# ------------------------------------------------------------------
class DynFrozen:
    """RREF of M = G_b P T^{-1}:  information set + dynamic-frozen relations."""

    __slots__ = ("J", "rel", "R")

    def __init__(self, M: np.ndarray):
        R, piv = gf2_rref(M)
        k = len(piv)
        self.R = R[:k]
        self.J: List[int] = list(piv)
        # column j not a pivot  ->  u_j is determined by the pivots r with
        # piv[r] < j and R[r, j] = 1  (dynamic frozen; empty => static frozen)
        self.rel: Dict[int, List[int]] = {}
        n = M.shape[1]
        pivset = set(piv)
        for j in range(n):
            if j in pivset:
                continue
            self.rel[j] = [piv[r] for r in range(k) if piv[r] < j and R[r, j]]

    @property
    def static_frozen(self) -> List[int]:
        return [j for j, s in self.rel.items() if not s]

    @property
    def dyn_frozen(self) -> List[int]:
        return [j for j, s in self.rel.items() if s]


def dyn_frozen_of(Gb: np.ndarray, P: np.ndarray, T: np.ndarray) -> DynFrozen:
    Tinv, piv = gf2_rref(np.concatenate([gf2(T), np.eye(T.shape[0], dtype=np.uint8)], 1))
    assert len(piv) == T.shape[0], "transform is singular"
    return DynFrozen(gf2_matmul(gf2_matmul(gf2(Gb), gf2(P)), Tinv[:, T.shape[0]:]))


# ------------------------------------------------------------------
# objective
# ------------------------------------------------------------------
def phi(df: DynFrozen, z: List[float], lam: float = 0.0,
        window: Optional[Tuple[int, int]] = None) -> float:
    """Phi = sum_{i in J} Z_i  -  lam * sum_{f in DF cap window} Delta_f."""
    val = sum(z[i] for i in df.J)
    if lam:
        lo, hi = window if window else (0, len(z))
        val -= lam * sum(sum(1.0 - z[i] for i in s)
                         for f, s in df.rel.items() if s and lo <= f < hi)
    return val


# ------------------------------------------------------------------
# self-test
# ------------------------------------------------------------------
def self_test(m: int = 5) -> None:
    n = 1 << m
    T = pruned_transform(m)
    assert np.array_equal(T, gf2(polar_matrix(m))), "unpruned T != F^(x m)"

    z0 = 0.37
    ref = [z0]
    for _ in range(m):
        ref = [v for x in ref for v in (2 * x - x * x, x * x)]
    got = z_pruned(m, (), z0)
    assert max(abs(a - b) for a, b in zip(ref, got)) < 1e-12, "Z recursion mismatch"

    # every single-kernel pruning is a rank-one change of exactly one layer,
    # keeps T invertible, and strictly lowers total row weight
    w0 = int(T.sum())
    for t, q in all_prunes(m):
        Tp = pruned_transform(m, [(t, q)])
        assert int(Tp.sum()) < w0
        aug = np.concatenate([Tp, np.eye(n, dtype=np.uint8)], 1)
        _, piv = gf2_rref(aug)
        assert len(piv) == n
    # pruning conserves total Z-sum mass ordering: sum over all i is unchanged
    # only for the trivial case; just check values stay in [0,1]
    for t, q in all_prunes(m)[:8]:
        assert all(0.0 <= v <= 1.0 for v in z_pruned(m, [(t, q)], z0))
    # sum conservation: each kernel maps (x,y) -> (x+y-xy, xy), sum x+y; the
    # pruned kernel maps (x,y) -> (x,y), sum x+y.  So sum_i Z_i = N*z0 always.
    import random as _r
    rng = _r.Random(11); R = all_prunes(m)
    for _ in range(50):
        e = rng.sample(R, rng.randint(1, max(2, len(R) // 4)))
        assert abs(sum(z_pruned(m, e, z0)) - n * z0) < 1e-9, "sum(Z) not conserved"

    # orientation pin: eBCH[32,16] under its equation permutation must give
    # exactly RM(2,5) = { a : popcount(a) >= 3 }
    if m == 5:
        root = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                            "project", "POD")
        Gb = gf2(parse_mat(os.path.join(root, "eBCH_m5_t3.matrix")))
        P = np.ascontiguousarray(gf2(parse_mat(
            os.path.join(root, "eBCH_m5_t3_P.matrix"))).T)
        rm = [a for a in range(n) if bin(a).count("1") >= 3]
        assert dyn_frozen_of(Gb, P, T).J == rm, "P orientation / info set mismatch"

    print(f"[self-test] m={m}: transform, Z-recursion, invertibility, "
          f"RM(2,5) anchor  OK ({len(all_prunes(m))} single prunings)")


# ------------------------------------------------------------------
def parse_mat(p: str) -> np.ndarray:
    tok = open(p).read().split()
    r, c = int(tok[0]), int(tok[1]); body = tok[2:]
    if len(body) == r:
        return np.array([[int(x) for x in row] for row in body], np.uint8)
    return np.array([int(x) for x in body], np.uint8).reshape(r, c)


def main() -> int:
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                        "project", "POD")
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default=os.path.join(root, "eBCH_m5_t3.matrix"))
    ap.add_argument("--perm", default=os.path.join(root, "eBCH_m5_t3_P.matrix"))
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--snr", type=float, default=3.0)
    ap.add_argument("--lam", type=float, default=0.0)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test(args.m)

    m, n = args.m, 1 << args.m
    Gb = gf2(parse_mat(args.gen))
    # the .matrix files store P with P[a, pi(a)] = 1, so G_b[:, pi] = G_b P^T
    P = np.ascontiguousarray(gf2(parse_mat(args.perm)).T)
    k = Gb.shape[0]
    assert Gb.shape[1] == n and P.shape == (n, n)
    z0 = z0_awgn(args.snr, k / n)

    T0 = pruned_transform(m)
    z_0 = z_pruned(m, (), z0)
    df0 = dyn_frozen_of(Gb, P, T0)
    phi0 = phi(df0, z_0, args.lam)

    print(f"=== baseline (standard kernel, T = F^(x{m})) ===")
    print(f"  code [{n},{k}]  Eb/N0 {args.snr} dB  z0={z0:.6f}  lambda={args.lam}")
    print(f"  info set J   : {df0.J}")
    print(f"  Phi          : {phi0:.6f}")
    print(f"  static frozen: {len(df0.static_frozen)}   "
          f"dynamic frozen: {len(df0.dyn_frozen)}")
    free0 = sum(sorted(z_0)[:k])
    print(f"  free-J optimum (k most reliable, ignores the code): {free0:.6f}")
    print(f"  alignment gap  Phi - free-J opt : {phi0 - free0:+.6f}")

    rows = []
    for t, q in all_prunes(m):
        a0, a1 = pair_of(m, t, q)
        df = dyn_frozen_of(Gb, P, pruned_transform(m, [(t, q)]))
        zz = z_pruned(m, [(t, q)], z0)
        rows.append(((t, q), (a0, a1), phi(df, zz, args.lam) - phi0,
                     sum(sorted(zz)[:k]) - free0, df, zz))

    rows.sort(key=lambda r: r[2])
    print(f"\n=== {len(rows)} single-kernel prunings, by influence delta_e ===")
    print(f"{'t':>2} {'q':>3}  {'pair':>9}  {'delta_e':>11}  {'delta_free':>11}  "
          f"{'|J diff|':>8}  {'#DF':>4}")
    for (t, q), (a0, a1), d, dfree, df, _ in rows[:args.top]:
        print(f"{t:>2} {q:>3}  {a0:>4},{a1:<4}  {d:>+11.6f}  {dfree:>+11.6f}  "
              f"{len(set(df.J) ^ set(df0.J)):>8}  {len(df.dyn_frozen):>4}")

    imp = [r for r in rows if r[2] < -1e-12]
    freeimp = [r for r in rows if r[3] < -1e-12]
    zero = sum(1 for r in rows if abs(r[2]) < 1e-12)
    print(f"\n  sum(Z) is conserved at N*z0 = {n * z0:.6f} for every pruning")
    print(f"  forced-J improvements : {len(imp)} / {len(rows)}   "
          f"(neutral: {zero})")
    print(f"  free-J   improvements : {len(freeimp)} / {len(rows)}   "
          f"min delta {min(r[3] for r in rows):+.6f}")
    if imp:
        (t, q), _, d, _, df, _ = imp[0]
        print(f"  best forced-J (t={t}, q={q})  delta={d:+.6f}\n    J = {df.J}")
    print(f"  worst delta: {rows[-1][2]:+.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
