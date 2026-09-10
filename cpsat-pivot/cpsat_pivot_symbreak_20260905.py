#!/usr/bin/env python3
"""
cpsat_pivot_symbreak  --  two additions aimed at the sweep's UNKNOWN cluster.

Both are OPT-IN; nothing here changes the behaviour of
cpsat_pivot_reduced_20260902.py, which the running sweep imports.

--------------------------------------------------------------------
(A) NARROW PREFIX CUTS FROM V_{p0}          (attacks weak propagation)
--------------------------------------------------------------------
The profile forces W[:,j] = 0 for every j < p0, i.e.

    pi(G_b) g_j = 0     for all j < p0,   g_j := G_p[:,j]

and therefore, by linearity, for EVERY vector of the span

    V_{p0} := span{ g_0, ..., g_{p0-1} } .

So each v in V_{p0} yields a valid cut, at width wt(v):

    XOR_{a : v[a] = 1}  c[pi(a)]  =  0        for every codeword c of C_b.

The 15 basis columns have widths [64, 32x4, 16x6, 8x4] -- a width-64 XOR
implies nothing until 63 of its literals are fixed, so most of what we
currently feed the solver is inert.  V_15 contains 120 vectors of weight 8
(= d_min(V_15) = 2^(m-nu(p0)), the prefix theorem), so taking the low-weight
vectors of the SPAN instead of the basis buys 30x more cuts, all at the
minimum possible width.

--------------------------------------------------------------------
(B) SYMMETRY BREAKING OVER Aut(C_b)         (attacks proof size)
--------------------------------------------------------------------
If tau in Aut(C_b) then G_b[:,tau] = A G_b for some A in GL(k), so with
pi' = tau o pi we get W' = A W.  LEFT multiplication by GL(k) does not change
the column pivot profile, hence pi' is a witness iff pi is.  Refutations
therefore come in orbits and CDCL re-derives each one separately.

We do NOT compute Aut(C_b) (code equivalence is hard).  For an extended
primitive BCH code of length 2^m = 64 (= extended cyclic of length 63) the
generators are given by the construction:

    shift    : i -> 1 + (i mod 63)          on coords 1..63, coord 0 fixed
    frobenius: i -> 2^j * i mod 63

Each is verified by ONE rank computation.  <shift, frob> has order 378 and
acts on {1..63} as i -> a*i + b (a in <2>, b in Z_63): transitive, point
stabiliser of order 6.

Breaking, with sigma := pi^{-1} (the action is sigma' = sigma o tau^{-1},
i.e. precomposition, so tau permutes sigma's ARGUMENTS):

    (1) transitivity => we may demand  sigma(1) = min{ sigma(v) : v = 1..63 }
        -- 62 plain inequalities, kills the factor 63.
    (2) the residual order-6 stabiliser of the point 1 is broken by requiring
        sigma to be lex-least among its 6 images -- 5 lex chains.

Sound for refutation: every orbit keeps its lex-leader, so UNSAT on the broken
model implies UNSAT on the original; and any SAT model is still a true witness.
"""
from __future__ import annotations

import argparse
import itertools
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from ortools.sat.python import cp_model

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, gf2_rank, polar_matrix, column_pivot_profile_gf2)
from cpsat_pivot_reduced_20260902 import (
    add_span_membership, add_u_certificate, add_xor_eq, build_w_dense,
    build_w_staged_kernel, infer_kernels, low_weight_codewords)


# ==================================================================
# (A) narrow cuts
# ==================================================================
def low_weight_span_vectors(Gp: np.ndarray, p0: int, max_vectors: int,
                            max_weight: Optional[int] = None
                            ) -> List[np.ndarray]:
    """Lowest-weight nonzero vectors of V_{p0} = span of columns 0..p0-1."""
    Gp = gf2(Gp); n = Gp.shape[0]
    if p0 <= 0:
        return []
    cols = [int("".join(map(str, Gp[:, j])), 2) for j in range(p0)]
    found: List[Tuple[int, int]] = []
    cur = 0
    for msk in range(1, 1 << p0):                  # 2^p0 Gray code; p0 <= ~20
        cur ^= cols[(msk & -msk).bit_length() - 1]
        w = bin(cur).count("1")
        if max_weight is None or w <= max_weight:
            found.append((w, cur))
    found.sort()
    out = []
    for w, val in found[:max_vectors]:
        out.append(np.array([(val >> (n - 1 - t)) & 1 for t in range(n)], np.uint8))
    return out


def add_narrow_prefix_cuts(model, pi_perm, Gb, Gp, p_star, num_words,
                           num_vectors, true_lit) -> Dict[str, int]:
    """Cuts  XOR_{a: v[a]=1} c[pi(a)] = 0  for low-weight v in V_{p0}."""
    Gb = gf2(Gb); Gp = gf2(Gp)
    k, n = Gb.shape
    j0 = p_star[0]
    st = {"words": 0, "vectors": 0, "cuts": 0, "max_width": 0, "min_width": None}
    if j0 == 0 or num_words <= 0 or num_vectors <= 0:
        return st

    vecs = low_weight_span_vectors(Gp, j0, num_vectors)
    words = low_weight_codewords(Gb, num_words)
    st["words"], st["vectors"] = len(words), len(vecs)
    if not vecs or not words:
        return st
    widths = [int(v.sum()) for v in vecs]
    st["max_width"], st["min_width"] = max(widths), min(widths)

    supports = [np.flatnonzero(v).tolist() for v in vecs]
    for wi, wc in enumerate(words):
        wlist = [int(t) for t in wc]
        pw = {}
        for a in range(n):
            bv = model.NewBoolVar(f"nw[{wi},{a}]")
            model.AddElement(pi_perm[a], wlist, bv)
            pw[a] = bv
        for supp in supports:
            add_xor_eq(model, [pw[a] for a in supp], 0, true_lit)
            st["cuts"] += 1
    return st


# ==================================================================
# (B) automorphisms + symmetry breaking
# ==================================================================
def is_automorphism(Gb: np.ndarray, perm: Sequence[int]) -> bool:
    Gb = gf2(Gb)
    return gf2_rank(np.vstack([Gb, Gb[:, list(perm)]])) == gf2_rank(Gb)


def bch_automorphism_generators(Gb: np.ndarray, m: int = 6
                                ) -> List[Tuple[str, List[int]]]:
    """
    Construction-given automorphisms of an extended primitive BCH code of
    length n = 2^m: the cyclic shift and the Frobenius multipliers, on the
    coordinates 1..n-1 with coordinate 0 fixed.  Each is VERIFIED.
    """
    n = 1 << m
    c = n - 1                                        # cyclic length 63
    gens: List[Tuple[str, List[int]]] = []
    shift = [0] + [1 + ((i + 1) % c) for i in range(c)]
    if is_automorphism(Gb, shift):
        gens.append(("shift", shift))
    for j in range(1, m):
        a = pow(2, j, c)
        idx = [0] + [1 + ((a * i) % c) for i in range(c)]
        if is_automorphism(Gb, idx):
            gens.append((f"frob^{j} (x{a})", idx))
    return gens


def close_group(gens: Sequence[Sequence[int]], n: int, cap: int = 500000
                ) -> List[Tuple[int, ...]]:
    ident = tuple(range(n))
    seen = {ident}
    frontier = [ident]
    gl = [tuple(g) for g in gens]
    while frontier and len(seen) < cap:
        cur = frontier.pop()
        for g in gl:
            nxt = tuple(g[x] for x in cur)
            if nxt not in seen:
                seen.add(nxt); frontier.append(nxt)
    return sorted(seen)


def add_symmetry_breaking(model, pi_perm, sigma_perm, group: Sequence[Tuple[int, ...]],
                          n: int, mode: str = "structured") -> Dict[str, int]:
    """
    sigma = pi^{-1}; the group acts by sigma' = sigma o tau^{-1}.

    "structured": use transitivity on {1..n-1} to demand sigma(1) minimal
                  (n-2 inequalities), then lex-break the point stabiliser.
    "leader":     full lex-leader against every non-identity group element.
    """
    st = {"group_order": len(group), "ineqs": 0, "lex_chains": 0}
    if len(group) <= 1:
        return st

    if mode == "leader":
        for tau in group[1:]:
            _add_lex_le(model, [sigma_perm[v] for v in range(n)],
                        [sigma_perm[tau[v]] for v in range(n)], n)
            st["lex_chains"] += 1
        return st

    # (1) transitive part: sigma(1) is the minimum over the orbit of 1
    orbit = sorted({tau[1] for tau in group})
    for v in orbit:
        if v != 1:
            model.Add(sigma_perm[1] <= sigma_perm[v])
            st["ineqs"] += 1
    # (2) residual point stabiliser of 1
    stab = [tau for tau in group if tau[1] == 1 and tau != tuple(range(n))]
    for tau in stab:
        _add_lex_le(model, [sigma_perm[v] for v in range(n)],
                    [sigma_perm[tau[v]] for v in range(n)], n)
        st["lex_chains"] += 1
    st["stabiliser"] = len(stab) + 1
    return st


def _add_lex_le(model, xs, ys, n) -> None:
    """xs <=_lex ys, via the standard prefix-equality chain."""
    eq = model.NewBoolVar("lexeq0"); model.Add(eq == 1)
    for i in range(n):
        model.Add(xs[i] <= ys[i]).OnlyEnforceIf(eq)
        if i + 1 < n:
            nxt = model.NewBoolVar(f"lexeq{i+1}")
            model.Add(xs[i] == ys[i]).OnlyEnforceIf(nxt)
            model.AddImplication(nxt, eq)
            eq = nxt


# ==================================================================
# model
# ==================================================================
def build_model(Gb, Gp, p_star, *, kernel_string=None, w_encoding="auto",
                certificate="full", min_weight_words=0, narrow_vectors=0,
                symmetry="none", m=6, branch="popcount", verbose=True):
    Gb = gf2(Gb); Gp = gf2(Gp); p_star = [int(x) for x in p_star]
    k, n = Gb.shape
    assert Gp.shape == (n, n) and len(p_star) == k

    ker = infer_kernels(Gp, kernel_string)
    w_encoding = w_encoding.lower()
    if w_encoding == "auto":
        w_encoding = "staged" if ker is not None else "dense"

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
        model.Add(sum(b[r, a] for a in range(n)) == int(Gb[r].sum()))

    constrained_cols = list(range(p_star[-1] + 1))
    if w_encoding == "staged":
        w, wst = build_w_staged_kernel(model, b, ker[0], ker[1], k, n,
                                       constrained_cols, true_lit)
    else:
        w, wst = build_w_dense(model, b, Gp, k, n, constrained_cols, true_lit)

    span_st = add_span_membership(model, w, k, p_star, true_lit, prefix_only=False)

    cert_st = {"z_vars": 0}
    if certificate in ("pruned", "full"):
        u = {(i, r): model.NewBoolVar(f"u[{i},{r}]") for i in range(k) for r in range(k)}
        always_zero = {(r, p): bool(np.all(Gb[r] == 0)) for r in range(k) for p in p_star}
        cert_st = add_u_certificate(model, w, u, k, p_star, true_lit,
                                    pruned=(certificate == "pruned"),
                                    always_zero=always_zero)

    cut_st = add_narrow_prefix_cuts(model, pi, Gb, Gp, p_star,
                                    min_weight_words, narrow_vectors, true_lit)

    sym_st = {"group_order": 1}
    if symmetry != "none":
        sigma = [model.NewIntVar(0, n - 1, f"sig[{v}]") for v in range(n)]
        model.AddInverse(pi, sigma)
        gens = bch_automorphism_generators(Gb, m)
        grp = close_group([g[1] for g in gens], n) if gens else []
        sym_st = add_symmetry_breaking(model, pi, sigma, grp, n, mode=symmetry)
        sym_st["generators"] = [g[0] for g in gens]

    order = sorted(range(n), key=lambda a: -bin(a).count("1")) if branch == "popcount" \
        else list(range(n))
    model.AddDecisionStrategy([pi[a] for a in order],
                              cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)

    if verbose:
        print(f"[symbreak] n={n} k={k} W={w_encoding} cert={certificate}")
        print(f"           narrow cuts : {cut_st}")
        print(f"           symmetry    : {sym_st}")
    return model, {"pi": pi, "stats": {"w": wst, "span": span_st, "cert": cert_st,
                                       "cuts": cut_st, "sym": sym_st}}


def solve(Gb, Gp, p_star, *, time_limit_sec=600.0, workers=8, verbose=True, **kw):
    Gb = gf2(Gb); Gp = gf2(Gp); n = Gb.shape[1]
    model, v = build_model(Gb, Gp, p_star, verbose=verbose, **kw)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(time_limit_sec)
    s.parameters.num_search_workers = int(workers)
    st = s.Solve(model)
    name = s.StatusName(st)
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return "sat", [int(s.Value(v["pi"][a])) for a in range(n)]
    if name == "INFEASIBLE":
        return "infeasible", None
    return "unknown", None
