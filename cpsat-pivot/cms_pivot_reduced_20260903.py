#!/usr/bin/env python3
"""
CryptoMiniSat encoding of the inverse pivot-profile problem, with the
2026-09-02 reductions ported over from cpsat_pivot_reduced_20260902.py.

    PivotProfile( pi(G_b) G_p ) = p_star ,   pi(G_b) = G_b[:, pi]   over F_2

What this adds on top of cpsat_pivot_crypto_minixor.py (which is left
untouched):

  (A) PARTIAL MIN-WEIGHT CUTS AS NATIVE XOR CLAUSES  (--min-weight-words N)
      For j < p_star[0] the profile forces W[:,j] = 0, i.e.
          pi(C_b) _|_ ghat_j ,     ghat_j = (G_p^{-T})[:, j],
      so for every low-weight codeword w_c of C_b

          XOR_{a in supp(ghat_j)}  w_c[pi_a]  =  0 .

      With the permutation selectors P[a][t] already in the encoding,
          w_c[pi_a] = XOR_{t in supp(w_c)} P[a][t]   =: pw[c][a]
      (one XOR clause of width |supp(w_c)|, exactly the same trick the
      channelling b[r,a] already uses), and each cut is then ONE XOR
      clause of width |supp(ghat_j)| over the pw variables.

      This is the ideal shape for CryptoMiniSat: the cuts land directly
      in the Gauss-Jordan matrix, so an inconsistent prefix-orthogonality
      system is refuted by ELIMINATION rather than by branching.  In
      CP-SAT the same cuts took the eGolay optimal-target UNSAT proof
      from 477 s to 0.35 s; here they cost no CNF at all.

  (B) STAGED MULTI-KERNEL BUTTERFLY W  (--w-encoding staged)
      G_p[i,j] = prod_s K_s[d_s(i)][d_s(j)] factors into S local stages,
          x^{(s+1)}[r,i] = XOR_e K_s[e][d_s(i)] . x^{(s)}[r, i|d_s:=e],
      each a short XOR clause.  Standard Arikan F^{x m} is the special
      case radices=[2]*m.  Keeps the parity matrix sparse for multikernel
      G_p, where the dense per-column encoding would produce wide rows.

  (C) CERTIFICATE LEVELS  (--certificate span-only|full)
      'full'      : span membership + U W[:,p_star] = I_k   (as before)
      'span-only' : span membership only -- drops the k^3 AND-gate layer.
                    A RELAXATION: span-only UNSAT  =>  true UNSAT (sound
                    fast filter); span-only SAT proves nothing.

  (D) --prepass-seconds 0  to skip the random-sampling pre-pass when the
      expected answer is UNSAT (sweeps), instead of burning 10 s/instance.

Usage
-----
    python cms_pivot_reduced_20260903.py inst.npz --min-weight-words 48
    python cms_pivot_reduced_20260903.py inst.npz --certificate span-only \
           --prepass-seconds 0 --time 600
    python cms_pivot_reduced_20260903.py inst.npz --kernels 657,23,23,23 \
           --w-encoding staged
    python cms_pivot_reduced_20260903.py inst.npz --dimacs out.cnf   # export
"""
from __future__ import annotations

import argparse
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

from cpsat_pivot_crypto_minixor import (
    CnfXor, gf2, gf2_matmul, gf2_rank, gf2_rref, column_pivot_profile_gf2,
    polar_matrix, is_standard_polar_matrix, polar_log2_size,
    dualize_instance, primal_pi_from_dual_sigma,
    random_sampling_prepass, decode_pi, verify_pi,
    HAVE_PYCRYPTOSAT,
)
try:
    from pycryptosat import Solver as CMSSolver
except ImportError:
    CMSSolver = None

from cpsat_pivot_reduced_20260902 import (
    infer_kernels, multikernel_from_mats, low_weight_codewords, KERNELS,
)


# ============================================================
# encoding
# ============================================================

def encode_reduced(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    kernel_string: Optional[str] = None,
    w_encoding: str = "auto",
    certificate: str = "full",
    min_weight_words: int = 0,
    verbose: bool = True,
) -> Tuple[CnfXor, Dict]:
    Gb, Gp = gf2(Gb), gf2(Gp)
    p_star = [int(x) for x in p_star]
    k, n = Gb.shape
    assert Gp.shape == (n, n) and len(p_star) == k
    assert sorted(p_star) == p_star and len(set(p_star)) == k
    assert gf2_rank(Gb) == k and gf2_rank(Gp) == n
    if certificate not in ("span-only", "full"):
        raise ValueError("certificate must be span-only|full")

    ker = infer_kernels(Gp, kernel_string)
    w_encoding = w_encoding.lower()
    if w_encoding == "auto":
        w_encoding = "staged" if ker is not None else "dense"
    if w_encoding == "staged" and ker is None:
        raise ValueError("--w-encoding staged needs --kernels or a Kronecker G_p")

    cnf = CnfXor()
    stats = {"and_gates": 0, "mw_cuts": 0, "mw_words": 0, "staged_vars": 0}

    # -- 1. permutation selectors, exactly-one per row and column --------
    P = [[cnf.new_var() for _ in range(n)] for _ in range(n)]
    for a in range(n):
        cnf.add_exactly_one([P[a][t] for t in range(n)])
    for t in range(n):
        cnf.add_exactly_one([P[a][t] for a in range(n)])

    # -- 2. channelling  b[r,a] = XOR_{t in supp(Gb[r,:])} P[a][t] -------
    def select(vec: np.ndarray, a: int) -> int:
        """fresh var == vec[pi_a], one native XOR clause."""
        supp = np.flatnonzero(vec).tolist()
        v = cnf.new_var()
        cnf.add_xor([v] + [P[a][t] for t in supp], False)
        return v

    b = {(r, a): select(Gb[r], a) for r in range(k) for a in range(n)}

    # -- 3. W ------------------------------------------------------------
    p_last = p_star[-1]
    cols = list(range(p_last + 1))
    w: Dict[Tuple[int, int], int] = {}

    if w_encoding == "staged":
        mats, radices = ker
        x_prev = {(r, a): b[r, a] for r in range(k) for a in range(n)}
        stride = 1
        for s, K in enumerate(mats):
            m = radices[s]
            x_next = {}
            for r in range(k):
                for i in range(n):
                    d = (i // stride) % m
                    base = i - d * stride
                    srcs = [x_prev[r, base + e * stride]
                            for e in range(m) if int(K[e, d]) == 1]
                    if len(srcs) == 1:
                        x_next[r, i] = srcs[0]
                    else:
                        out = cnf.new_var()
                        stats["staged_vars"] += 1
                        cnf.add_xor([out] + srcs, False)
                        x_next[r, i] = out
            x_prev = x_next
            stride *= m
        w = {(r, j): x_prev[r, j] for r in range(k) for j in cols}
    else:
        for j in cols:
            supp_j = np.flatnonzero(Gp[:, j]).tolist()
            for r in range(k):
                v = cnf.new_var()
                w[r, j] = v
                cnf.add_xor([v] + [b[r, a] for a in supp_j], False)

    # -- 4. span membership for non-pivot columns ------------------------
    pivot_set = set(p_star)
    for j in cols:
        if j in pivot_set:
            continue
        prev = [p for p in p_star if p < j]
        if not prev:
            for r in range(k):
                cnf.add_clause([-w[r, j]])
            continue
        lam = cnf.new_vars(len(prev))
        for r in range(k):
            ys = []
            for l, p in enumerate(prev):
                y = cnf.new_var()
                cnf.add_and_gate(y, lam[l], w[r, p])
                stats["and_gates"] += 1
                ys.append(y)
            cnf.add_xor([w[r, j]] + ys, False)

    # -- 5. certificate  U W[:,p_star] = I_k  (optional) -----------------
    u = None
    if certificate == "full":
        u = [[cnf.new_var() for _ in range(k)] for _ in range(k)]
        for i in range(k):
            for li, p in enumerate(p_star):
                zs = []
                for r in range(k):
                    z = cnf.new_var()
                    cnf.add_and_gate(z, u[i][r], w[r, p])
                    stats["and_gates"] += 1
                    zs.append(z)
                cnf.add_xor(zs, i == li)

    # -- 6. min-weight prefix cuts, NATIVE XOR ---------------------------
    j0 = p_star[0]
    if min_weight_words > 0 and j0 > 0:
        # BUG FIX (2026-09): the cut's support is column j of G_p directly,
        # not of G_p^{-T}.  Derivation: W[r,j] = XOR_a Gb[r,pi(a)] Gp[a,j], so
        # for a codeword c, W[:,j]=0 forces XOR_{a: Gp[a,j]=1} c[pi(a)] = 0.
        # The G_p^{-T} version was checked against a witness permutation and
        # found unsound (311/600 violations); see FORMULATIONS.md sec 3 and
        # cpsat_pivot_reduced_20260902.py's add_min_weight_prefix_cuts.
        words = low_weight_codewords(Gb, min_weight_words)
        stats["mw_words"] = len(words)
        stats["mw_min_wt"] = int(words[0].sum()) if words else None
        for wc in words:
            pw = {a: select(np.asarray(wc, np.uint8), a) for a in range(n)}
            for j in range(j0):
                supp_a = np.flatnonzero(Gp[:, j]).tolist()
                if supp_a:
                    cnf.add_xor([pw[a] for a in supp_a], False)
                    stats["mw_cuts"] += 1

    if verbose:
        print(f"[cms-reduced] n={n} k={k}  W={w_encoding}  cert={certificate}")
        print(f"              vars={cnf.nvars}  CNF={len(cnf.clauses)}  XOR={len(cnf.xors)}")
        print(f"              AND gates={stats['and_gates']}  staged={stats['staged_vars']}")
        print(f"              min-weight: {stats['mw_words']} words "
              f"(min wt {stats.get('mw_min_wt')}), {stats['mw_cuts']} native-XOR cuts")

    return cnf, {"P": P, "u": u, "n": n, "k": k, "stats": stats,
                 "certificate": certificate}


# ============================================================
# solving
# ============================================================

def solve_reduced_cms(
    Gb, Gp, p_star, *,
    kernel_string=None, w_encoding="auto", certificate="full",
    min_weight_words=0, time_limit_sec=3600.0, threads=8,
    prepass_seconds=10.0, prepass_tries=200000,
    dimacs_path=None, sat_verbosity=0, verbose=True,
) -> Tuple[str, Optional[List[int]]]:
    """Returns (status, pi) with status in {'sat','infeasible','unknown'}."""
    Gb, Gp = gf2(Gb), gf2(Gp)
    p_star = [int(x) for x in p_star]
    k, n = Gb.shape

    if prepass_seconds > 0:
        hit = random_sampling_prepass(Gb, Gp, p_star, max_tries=prepass_tries,
                                      max_seconds=prepass_seconds, verbose=verbose)
        if hit is not None:
            return "sat", hit

    cnf, vm = encode_reduced(Gb, Gp, p_star, kernel_string=kernel_string,
                             w_encoding=w_encoding, certificate=certificate,
                             min_weight_words=min_weight_words, verbose=verbose)
    if dimacs_path:
        cnf.write_dimacs(dimacs_path)
        print(f"wrote {dimacs_path}")

    if not HAVE_PYCRYPTOSAT or CMSSolver is None:
        raise RuntimeError("pycryptosat not installed")

    s = CMSSolver(threads=int(threads), time_limit=float(time_limit_sec),
                  verbose=int(sat_verbosity))
    for c in cnf.clauses:
        if not c:
            return "infeasible", None
        s.add_clause(c)
    for vs, rhs in cnf.xors:
        s.add_xor_clause(vs, rhs)

    sat, model = s.solve()
    if sat is True:
        pi = decode_pi(model, vm["P"], n)
        if certificate == "span-only":
            # relaxation: SAT proves nothing -- verify directly.
            if column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp)) != p_star:
                return "unknown", None
        return "sat", pi
    if sat is False:
        return "infeasible", None
    return "unknown", None


# ============================================================
# CLI
# ============================================================

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--kernels", default=None)
    ap.add_argument("--w-encoding", choices=["auto", "staged", "dense"], default="auto")
    ap.add_argument("--certificate", choices=["span-only", "full"], default="full")
    ap.add_argument("--min-weight-words", type=int, default=0)
    ap.add_argument("--time", type=float, default=3600.0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--prepass-seconds", type=float, default=10.0)
    ap.add_argument("--dimacs", default=None)
    ap.add_argument("--sat-verbose", type=int, default=0)
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=False)
    Gb, Gp = gf2(d["Gb"]), gf2(d["Gp"])
    p_star = [int(x) for x in np.asarray(d["p_star"]).reshape(-1)]

    t0 = time.time()
    st, pi = solve_reduced_cms(
        Gb, Gp, p_star, kernel_string=args.kernels, w_encoding=args.w_encoding,
        certificate=args.certificate, min_weight_words=args.min_weight_words,
        time_limit_sec=args.time, threads=args.threads,
        prepass_seconds=args.prepass_seconds, dimacs_path=args.dimacs,
        sat_verbosity=args.sat_verbose, verbose=True,
    )
    dt = time.time() - t0
    print(f"\nstatus: {st}   ({dt:.2f}s)")
    if st == "sat":
        print("pi =", pi)
        print("check:", verify_pi(Gb, Gp, p_star, pi, "cms-reduced"))
        return 0
    if st == "infeasible":
        print("p_star provably UNREACHABLE"
              + ("  (span-only relaxation => sound)" if args.certificate == "span-only" else ""))
        return 1
    return 2


if __name__ == "__main__":
    sys.exit(main())
