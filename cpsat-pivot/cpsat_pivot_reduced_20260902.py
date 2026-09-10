#!/usr/bin/env python3
"""
Reduced CP-SAT encoding for the inverse pivot-profile problem

    PivotProfile( pi(G_b) G_p ) = p_star ,     pi(G_b) = G_b[:, pi]  (over F_2)

Changes vs cpsat_pivot_reconstruct_duality.py  (2026-09-02):

  (2) STAGED MULTI-KERNEL BUTTERFLY W-ENCODING  (--w-encoding staged)
      For G_p a Kronecker product of small kernels
          G_p[i,j] = prod_s  K_s[ d_s(i) ][ d_s(j) ]
      (KernalManager convention: radices R = [m_0, ..., m_{S-1}],
       d_s(i) = (i // prod(R[:s])) % m_s, stage 0 = least significant),
      W = pi(G_b) G_p factors into S local stages, stage s mixing the
      m_s entries that share every digit but digit s:
          x^{(s+1)}[r,i] = XOR_e  K_s[e][d_s(i)] . x^{(s)}[r, i|d_s:=e]
      Cost drops from  O(k n * avg col-support of G_p)  (dense) to
      O(k n * sum_s m_s) = O(k n log n).  Standard Arikan F^{x m} is the
      special case R = [2]*m, K_s = [[1,0],[1,1]].

  (3) SUPPORT-PRUNED / OPTIONAL U-CERTIFICATE  (--certificate ...)
      The pivot-only certificate needs (a) U W[:,p_star] = I_k to certify
      the k pivot columns are jointly independent, and (b) span
      membership for the non-pivots.  (b) is U-free already.  (a) is an
      unknown x unknown product and is inherently ~k^3; a purely local
      U-free encoding of F_2 invertibility does not exist (the "residual
      != 0" trick needs FOR ALL lambda, which the solver does not check).
      What we do:
        * 'span-only'   : keep only (b).  NOT sound on its own -- it
                          allows the k pivot columns to be dependent.
                          Use together with (4) as a fast *filter*:
                          span-only UNSAT  => the instance is UNSAT.
        * 'pruned'      : (a)+(b), but z[i,r,l] is created only when
                          w[r,p_l] is not structurally forced to 0.
        * 'full'        : (a)+(b), every z[i,r,l]  (= duality.py).

  (4) PARTIAL MINIMUM-WEIGHT CUTS  (--min-weight-words N)
      For j < p_star[0] the profile forces W[:,j] = 0, i.e.
          pi(C_b)  _|_  ghat_j ,     ghat_j = (G_p^{-T})[:, j].
      Every codeword of C_b, in particular a low-weight one w_c, must
      satisfy pi(w_c) . ghat_j = 0.  With  pi(w_c)[a] = XOR_t w_c[t] P[a][t]
      this is ONE parity constraint per (w_c, j), linear in the
      permutation selectors P -- a sound necessary condition that CP-SAT
      propagates directly on pi.  N low-weight codewords of C_b are found
      by enumeration (2^k <= 2^20) or random sampling.

Usage:
    python cpsat_pivot_reduced_20260902.py                 # self-tests
    python cpsat_pivot_reduced_20260902.py inst.npz --time 600 --workers 8
    python cpsat_pivot_reduced_20260902.py inst.npz --kernels 657,23,23,23 \
           --w-encoding staged --min-weight-words 24
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from ortools.sat.python import cp_model

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, gf2_rank, gf2_rref, gf2_nullspace,
    column_pivot_profile_gf2, polar_matrix, is_standard_polar_matrix,
    add_xor_eq, add_and_eq, popcount_descending_order,
    solve_reconstruction,
)

# 2x2 / 3x3 kernels, KernalManager names (src/ErrorCorrectionCode/KernalManager.cpp)
KERNELS: Dict[str, np.ndarray] = {
    "23":  np.array([[1, 0], [1, 1]], dtype=np.uint8),
    "753": np.array([[1, 1, 1], [1, 0, 1], [0, 1, 1]], dtype=np.uint8),
    "427": np.array([[1, 0, 0], [0, 1, 0], [1, 1, 1]], dtype=np.uint8),
    "657": np.array([[1, 1, 0], [1, 0, 1], [1, 1, 1]], dtype=np.uint8),
}


# ============================================================
# multi-kernel G_p
# ============================================================

def kernel_radices(kernel_string: str) -> List[int]:
    return [len(tok) for tok in kernel_string.split(",")]


def digit(index: int, radices: List[int], s: int) -> int:
    stride = 1
    for r in radices[:s]:
        stride *= r
    return (index // stride) % radices[s]


def multikernel_from_mats(mats: List[np.ndarray]) -> np.ndarray:
    """G_p[i,j] = prod_s K_s[d_s(i)][d_s(j)], stage 0 = least significant digit."""
    radices = [int(K.shape[0]) for K in mats]
    n = int(np.prod(radices))
    G = np.ones((n, n), dtype=np.uint8)
    for i in range(n):
        di = [digit(i, radices, s) for s in range(len(mats))]
        for j in range(n):
            dj = [digit(j, radices, s) for s in range(len(mats))]
            bit = 1
            for s in range(len(mats)):
                bit &= int(mats[s][di[s], dj[s]])
            G[i, j] = bit
    return G


def build_multikernel_gp(kernel_string: str) -> np.ndarray:
    return multikernel_from_mats([KERNELS[t] for t in kernel_string.split(",")])


def _kernels_match(Gp: np.ndarray, mats: List[np.ndarray]) -> bool:
    rad = [int(K.shape[0]) for K in mats]
    return int(np.prod(rad)) == Gp.shape[0] and np.array_equal(gf2(Gp), multikernel_from_mats(mats))


def infer_kernels(
    Gp: np.ndarray,
    kernel_string: Optional[str],
    kernels_explicit: Optional[List[np.ndarray]] = None,
) -> Optional[Tuple[List[np.ndarray], List[int]]]:
    """
    Return (kernel matrices, radices) reproducing Gp, else None.
    Tries: explicit list, --kernels string, standard polar, and the
    J-conjugates of the string kernels (row/col reversal, transpose) so
    that J.Gp / Gp.J / Gp^T of a multikernel are still recognised.
    """
    n = Gp.shape[0]
    if kernels_explicit is not None:
        if not _kernels_match(Gp, kernels_explicit):
            raise ValueError("explicit kernels do not reproduce Gp")
        return kernels_explicit, [int(K.shape[0]) for K in kernels_explicit]
    if kernel_string:
        base = [KERNELS[t] for t in kernel_string.split(",")]
        for variant in (base,
                        [np.flipud(K).copy() for K in base],     # J . Gp
                        [np.fliplr(K).copy() for K in base],      # Gp . J
                        [K.T.copy() for K in base],               # Gp^T
                        [np.flipud(K).T.copy() for K in base]):   # J . Gp^T
            if _kernels_match(Gp, variant):
                return variant, [int(K.shape[0]) for K in variant]
        raise ValueError("supplied --kernels (and its J-conjugates) do not reproduce Gp")
    if is_standard_polar_matrix(Gp):
        m = int(round(np.log2(n)))
        return [KERNELS["23"]] * m, [2] * m
    return None


# ============================================================
# W encodings
# ============================================================

def build_w_staged_kernel(
    model: cp_model.CpModel,
    b: Dict[Tuple[int, int], cp_model.IntVar],
    kernels: List[np.ndarray],
    radices: List[int],
    k: int,
    n: int,
    constrained_cols: List[int],
    true_lit: cp_model.IntVar,
) -> Tuple[Dict[Tuple[int, int], cp_model.IntVar], Dict[str, int]]:
    """
    W = pi(G_b) G_p by the staged Kronecker recursion, stage 0 = LSB.
        x^{(s+1)}[r,i] = XOR_e  K_s[e][d_s(i)] . x^{(s)}[r, i with digit s = e]
    Verified equal to gf2_matmul(B, build_multikernel_gp(...)).
    """
    x_prev: Dict[Tuple[int, int], cp_model.IntVar] = {
        (r, a): b[r, a] for r in range(k) for a in range(n)
    }
    n_vars = 0
    n_xor = 0
    stride = 1
    for s, K in enumerate(kernels):
        m = radices[s]
        x_next: Dict[Tuple[int, int], cp_model.IntVar] = {}
        for r in range(k):
            for i in range(n):
                d = (i // stride) % m
                base = i - d * stride
                srcs = [x_prev[r, base + e * stride] for e in range(m) if int(K[e, d]) == 1]
                if len(srcs) == 1:
                    x_next[r, i] = srcs[0]
                else:
                    out = model.NewBoolVar(f"xs[{s+1},{r},{i}]")
                    n_vars += 1
                    add_xor_eq(model, srcs + [out], 0, true_lit)
                    n_xor += 1
                    x_next[r, i] = out
        x_prev = x_next
        stride *= m

    w = {(r, j): x_prev[r, j] for r in range(k) for j in constrained_cols}
    return w, {"staged_vars": n_vars, "staged_xor": n_xor}


def build_w_dense(
    model: cp_model.CpModel,
    b: Dict[Tuple[int, int], cp_model.IntVar],
    Gp: np.ndarray,
    k: int,
    n: int,
    constrained_cols: List[int],
    true_lit: cp_model.IntVar,
) -> Tuple[Dict[Tuple[int, int], cp_model.IntVar], Dict[str, int]]:
    w: Dict[Tuple[int, int], cp_model.IntVar] = {}
    n_vars = 0
    for r in range(k):
        for j in constrained_cols:
            w[r, j] = model.NewBoolVar(f"w[{r},{j}]")
            n_vars += 1
            supp = [a for a in range(n) if int(Gp[a, j]) == 1]
            add_xor_eq(model, [b[r, a] for a in supp] + [w[r, j]], 0, true_lit)
    return w, {"dense_vars": n_vars}


# ============================================================
# certificate (span membership + optional U)
# ============================================================

def add_span_membership(
    model, w, k, p_star, true_lit, *, prefix_only: bool = False,
) -> Dict[str, int]:
    """
    (b): every non-pivot column j <= p_star[-1] lies in span of earlier pivots.

    prefix_only=True keeps ONLY the j < p_star[0] rows (W[:,j] = 0), i.e. the
    linear part -- no lambda*w gates.  Still a sound relaxation for the
    infeasible direction ("cuts-only" certificate).
    """
    pivot_set = set(p_star)
    st = {"span_lambda": 0, "span_and": 0, "span_zero": 0}
    for j in range(p_star[-1] + 1):
        if j in pivot_set:
            continue
        prev = [p for p in p_star if p < j]
        if not prev:
            for r in range(k):
                model.Add(w[r, j] == 0)
                st["span_zero"] += 1
            continue
        if prefix_only:
            continue
        lam = [model.NewBoolVar(f"lam[{l},{j}]") for l in range(len(prev))]
        st["span_lambda"] += len(prev)
        for r in range(k):
            ys = []
            for l, p in enumerate(prev):
                y = model.NewBoolVar(f"y[{l},{r},{j}]")
                add_and_eq(model, y, lam[l], w[r, p])
                st["span_and"] += 1
                ys.append(y)
            add_xor_eq(model, ys + [w[r, j]], 0, true_lit)
    return st


def add_u_certificate(
    model, w, u, k, p_star, true_lit, *, pruned: bool,
    always_zero,
) -> Dict[str, int]:
    """(a): U W[:,p_star] = I_k  via z[i,r,l] = u[i,r] AND w[r,p_l]."""
    st = {"z_vars": 0}
    for i in range(k):
        for li, p in enumerate(p_star):
            terms = []
            for r in range(k):
                if pruned and always_zero[r, p]:
                    continue                       # w[r,p] structurally 0
                z = model.NewBoolVar(f"z[{i},{r},{p}]")
                add_and_eq(model, z, u[i, r], w[r, p])
                st["z_vars"] += 1
                terms.append(z)
            add_xor_eq(model, terms, 1 if i == li else 0, true_lit)
    return st


# ============================================================
# partial minimum-weight cuts
# ============================================================

def low_weight_codewords(Gb: np.ndarray, want: int, *, rng_seed: int = 0) -> List[np.ndarray]:
    """
    Up to `want` distinct *minimum-weight* codewords of rowspace(Gb).

    k <= 20 : full 2^k enumeration.
    k >  20 : bounded low-order search -- XOR every combination of <= w_max
              RREF rows (over a few random column orderings of the RREF so the
              search reaches different low-weight words), increasing w_max until
              a stable minimum weight is seen.  This actually finds d_min for
              BCH-scale codes, unlike blind mask sampling.
    """
    import itertools as _it
    Gb = gf2(Gb)
    k, n = Gb.shape
    best = n + 1
    found: Dict[Tuple[int, ...], np.ndarray] = {}

    if k <= 20:
        for mask in range(1, 1 << k):
            c = np.zeros(n, dtype=np.uint8)
            for r in range(k):
                if (mask >> r) & 1:
                    c ^= Gb[r]
            wc = int(c.sum())
            if wc < best:
                best = wc; found = {tuple(c.tolist()): c}
            elif wc == best:
                found[tuple(c.tolist())] = c
        return list(found.values())[:want]

    rng = np.random.default_rng(rng_seed)
    bases = [gf2(Gb).copy()]
    for _ in range(6):                                   # random RREF re-bases
        R, piv = gf2_rref(Gb[:, rng.permutation(n)])
        bases.append(R[:len(piv)])
    w_max = 6
    while True:
        for B in bases:
            kb = B.shape[0]
            for w in range(1, w_max + 1):
                for combo in _it.combinations(range(kb), w):
                    c = np.zeros(n, dtype=np.uint8)
                    for r in combo:
                        c ^= B[r]
                    wc = int(c.sum())
                    if 0 < wc < best:
                        best = wc; found = {tuple(c.tolist()): c}
                    elif wc == best:
                        found[tuple(c.tolist())] = c
                        if len(found) >= max(want, 64):
                            break
        if len(found) >= want or w_max >= 8:
            break
        w_max += 1
    return list(found.values())[:want]


def add_min_weight_prefix_cuts(
    model, pi_perm, Gb, Gp, p_star, num_words, true_lit,
) -> Dict[str, int]:
    """
    W = pi(G_b) G_p, so  W[r,j] = XOR_a G_b[r, pi(a)] G_p[a,j].  For j < p_star[0]
    the profile forces W[:,j] = 0, hence for EVERY codeword c of C_b

        XOR_r alpha_r W[r,j] = XOR_a c[pi(a)] G_p[a,j] = 0
        i.e.   XOR_{a : G_p[a,j] = 1}  c[pi(a)]  =  0 .

    Take c = a low-weight codeword w_c.  c[pi(a)] is one AddElement(pi_a, w_c)
    per position a; then one parity per (w_c, j) over the support of COLUMN j of
    G_p.  Sound necessary condition, linear in pi.  (Earlier versions used
    G_p^{-T}[:,j] here -- that is NOT the valid support and over-prunes.)
    """
    k, n = gf2(Gb).shape
    Gp = gf2(Gp)
    j0 = p_star[0]
    if j0 == 0 or num_words <= 0:
        return {"mw_words": 0, "mw_cuts": 0}

    words = low_weight_codewords(Gb, num_words)
    st = {"mw_words": len(words), "mw_cuts": 0, "mw_min_wt": None}
    if words:
        st["mw_min_wt"] = int(words[0].sum())

    for wi, wc in enumerate(words):
        wlist = [int(v) for v in wc]                       # w_c[t] for t in [n]
        pw = {}                                            # pw[a] = w_c[pi_a]
        for a in range(n):
            v = model.NewBoolVar(f"pw[{wi},{a}]")
            model.AddElement(pi_perm[a], wlist, v)
            pw[a] = v
        for j in range(j0):
            supp_a = [a for a in range(n) if int(Gp[a, j]) == 1]
            if supp_a:
                add_xor_eq(model, [pw[a] for a in supp_a], 0, true_lit)
                st["mw_cuts"] += 1
    return st


# ============================================================
# model
# ============================================================

def build_reduced_model(
    Gb, Gp, p_star, *,
    kernel_string: Optional[str] = None,
    w_encoding: str = "auto",
    certificate: str = "pruned",
    min_weight_words: int = 0,
    branch: str = "popcount",
    verbose: bool = True,
):
    Gb = gf2(Gb); Gp = gf2(Gp); p_star = [int(x) for x in p_star]
    k, n = Gb.shape
    assert Gp.shape == (n, n) and len(p_star) == k
    assert sorted(p_star) == p_star and len(set(p_star)) == k
    assert gf2_rank(Gb) == k and gf2_rank(Gp) == n

    ker = infer_kernels(Gp, kernel_string)
    w_encoding = w_encoding.lower()
    if w_encoding == "auto":
        w_encoding = "staged" if ker is not None else "dense"
    if w_encoding == "staged" and ker is None:
        raise ValueError("--w-encoding staged needs --kernels or standard polar Gp")

    model = cp_model.CpModel()
    true_lit = model.NewBoolVar("TRUE"); model.Add(true_lit == 1)

    pi = [model.NewIntVar(0, n - 1, f"pi[{a}]") for a in range(n)]
    model.AddAllDifferent(pi)

    b: Dict[Tuple[int, int], cp_model.IntVar] = {}
    for r in range(k):
        row = [int(Gb[r, t]) for t in range(n)]
        for a in range(n):
            b[r, a] = model.NewBoolVar(f"b[{r},{a}]")
            if all(v == 0 for v in row):
                model.Add(b[r, a] == 0)
            elif all(v == 1 for v in row):
                model.Add(b[r, a] == 1)
            else:
                model.AddElement(pi[a], row, b[r, a])
        model.Add(sum(b[r, a] for a in range(n)) == int(Gb[r].sum()))  # weight cut

    constrained_cols = list(range(p_star[-1] + 1))
    if w_encoding == "staged":
        w, wst = build_w_staged_kernel(model, b, ker[0], ker[1], k, n, constrained_cols, true_lit)
    else:
        w, wst = build_w_dense(model, b, Gp, k, n, constrained_cols, true_lit)

    span_st = add_span_membership(model, w, k, p_star, true_lit,
                                  prefix_only=(certificate == "cuts-only"))

    cert_st = {"z_vars": 0}
    if certificate in ("pruned", "full"):
        u = {(i, r): model.NewBoolVar(f"u[{i},{r}]") for i in range(k) for r in range(k)}
        # structural zero test for w[r,p]: probe with a random b feed is unsafe;
        # be conservative -> only mark always-zero when Gb row r is all zero.
        always_zero = {(r, p): bool(np.all(Gb[r] == 0)) for r in range(k) for p in p_star}
        cert_st = add_u_certificate(model, w, u, k, p_star, true_lit,
                                    pruned=(certificate == "pruned"),
                                    always_zero=always_zero)
    elif certificate in ("span-only", "cuts-only"):
        u = None
    else:
        raise ValueError("certificate must be cuts-only|span-only|pruned|full")

    mw_st = add_min_weight_prefix_cuts(model, pi, Gb, Gp, p_star, min_weight_words, true_lit)

    order = popcount_descending_order(n) if branch == "popcount" else list(range(n))
    model.AddDecisionStrategy([pi[a] for a in order],
                              cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)

    if verbose:
        print(f"[reduced] n={n} k={k}  W={w_encoding}  certificate={certificate}")
        print(f"          {wst}")
        print(f"          span: {span_st}")
        print(f"          cert: {cert_st}")
        print(f"          min-weight cuts: {mw_st}")

    return model, {"pi": pi, "u": u if certificate in ("pruned","full") else None,
                   "w_encoding": w_encoding, "certificate": certificate}


def solve_reduced(
    Gb, Gp, p_star, *,
    kernel_string=None, w_encoding="auto", certificate="pruned",
    min_weight_words=0, time_limit_sec=300.0, workers=8, log=False, verbose=True,
) -> Tuple[str, Optional[List[int]]]:
    Gb = gf2(Gb); Gp = gf2(Gp); p_star = [int(x) for x in p_star]
    k, n = Gb.shape
    model, vars_ = build_reduced_model(
        Gb, Gp, p_star, kernel_string=kernel_string, w_encoding=w_encoding,
        certificate=certificate, min_weight_words=min_weight_words, verbose=verbose,
    )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_sec)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.log_search_progress = bool(log)
    st = solver.Solve(model)
    name = solver.StatusName(st)
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        pi_sol = [int(solver.Value(vars_["pi"][a])) for a in range(n)]
        actual = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi_sol], Gp))
        if actual != p_star:
            return ("unknown", None)
        return ("sat", pi_sol)
    if st == cp_model.INFEASIBLE:
        return ("infeasible", None)
    return ("unknown", None)


# ============================================================
# self-test / CLI
# ============================================================

def _random_instance(kernel_string: str, k: int, seed: int):
    rng = np.random.default_rng(seed)
    Gp = build_multikernel_gp(kernel_string)
    n = Gp.shape[0]
    while True:
        Gb = rng.integers(0, 2, (k, n), dtype=np.uint8)
        if gf2_rank(Gb) == k:
            break
    pi = list(rng.permutation(n))
    p_star = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp))
    return Gb, Gp, p_star, pi


def _staged_apply_np(B: np.ndarray, kernels, radices) -> np.ndarray:
    x = gf2(B).copy(); stride = 1
    n = B.shape[1]
    for s, K in enumerate(kernels):
        m = radices[s]; xn = np.zeros_like(x)
        for i in range(n):
            d = (i // stride) % m; base = i - d * stride
            acc = np.zeros(x.shape[0], np.uint8)
            for e in range(m):
                if K[e, d]:
                    acc ^= x[:, base + e * stride]
            xn[:, i] = acc
        x = xn; stride *= m
    return x & 1


def self_test() -> None:
    print("[self-test] staged-kernel W == dense matmul")
    for ks in ("23,23,23", "657,23,23,23", "753,23,23", "23,23,23,23,23"):
        Gb, Gp, p_star, pi = _random_instance(ks, k=3, seed=1)
        kers, rad = infer_kernels(Gp, ks)
        B = gf2(Gb)[:, pi]
        assert np.array_equal(_staged_apply_np(B, kers, rad), gf2_matmul(B, Gp)), ks
        print(f"   {ks:16s} OK")

    print("[self-test] solve small reachable instances")
    for ks, k in (("23,23,23", 3), ("657,23,23", 3), ("23,23,23,23", 5)):
        Gb, Gp, p_star, pi = _random_instance(ks, k, seed=7)
        stt, sol = solve_reduced(Gb, Gp, p_star, kernel_string=ks,
                                 w_encoding="staged", certificate="full",
                                 min_weight_words=4,
                                 time_limit_sec=30, workers=4, verbose=False)
        ok = stt == "sat" and column_pivot_profile_gf2(
            gf2_matmul(gf2(Gb)[:, sol], gf2(Gp))) == p_star
        print(f"   {ks:14s} k={k}: {stt}  {'OK' if ok else 'FAIL'}")
        assert ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", nargs="?", default=None)
    ap.add_argument("--kernels", default=None, help="e.g. 657,23,23,23")
    ap.add_argument("--w-encoding", choices=["auto", "staged", "dense"], default="auto")
    ap.add_argument("--certificate", choices=["cuts-only", "span-only", "pruned", "full"], default="pruned")
    ap.add_argument("--min-weight-words", type=int, default=0)
    ap.add_argument("--time", type=float, default=300.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--log", action="store_true")
    args = ap.parse_args()

    if args.npz is None:
        self_test()
        return 0

    d = np.load(args.npz, allow_pickle=False)
    Gb, Gp = gf2(d["Gb"]), gf2(d["Gp"])
    p_star = [int(x) for x in np.asarray(d["p_star"]).reshape(-1)]

    t0 = time.time()
    status, pi_sol = solve_reduced(
        Gb, Gp, p_star,
        kernel_string=args.kernels, w_encoding=args.w_encoding,
        certificate=args.certificate, min_weight_words=args.min_weight_words,
        time_limit_sec=args.time, workers=args.workers, log=args.log, verbose=True,
    )
    print(f"\nstatus: {status}   ({time.time()-t0:.2f}s)")
    if status == "sat":
        print("pi =", pi_sol)
        actual = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi_sol], Gp))
        print("profile check:", "PASS" if actual == p_star else "FAIL")
        return 0
    if status == "infeasible":
        print("p_star provably UNREACHABLE for this (Gb, Gp).")
        return 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
