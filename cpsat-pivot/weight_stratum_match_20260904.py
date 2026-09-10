#!/usr/bin/env python3
"""
weight_stratum_match  --  pi-free necessary conditions from stratum weight matching.

Flag identity (verified numerically):

    R_j := span{ row_j(G_p), ..., row_{n-1}(G_p) }
    dim( pi(C_b) cap R_j )  =  #{ pivots >= j }
    =>  j is a PIVOT  <=>  pi(C_b) meets the stratum  E_j := row_j(G_p) + R_{j+1}

Every nonzero codeword of pi(C_b) lies in exactly one stratum (its leading
index), and a permutation preserves weight.  So A_{C_b} must distribute over
the strata of the pivot set, giving pi-free necessary conditions:

  (F) FORBIDDEN INDEX
      supp(A_{E_j}) cap supp(A_{C_b}) = {}   =>  j can never be a pivot.
      Any p_star containing j is unreachable.

  (C) STRATUM CAPACITY
      If p_star = {p_0 < ... < p_{k-1}} then stratum E_{p_i} must host exactly
          N_i = 2^{k-1-i}
      codewords of pi(C_b) (the dimension jump), each of a weight in
      supp(A_{C_b}).  Hence
          sum_{w in supp(A_{C_b})} A_{E_{p_i}}(w)  >=  N_i .
      Violation => unreachable.  This bites hardest at the LATE pivots, where
      N_i is small but |E_{p_i}| = 2^{n-1-p_i} is small too.

  (T) GLOBAL TRANSPORT
      sum_i x[i,w] = A_{C_b}(w),  sum_w x[i,w] = N_i,  0 <= x[i,w] <= A_{E_{p_i}}(w).
      Infeasible transportation => unreachable.  Checked by the Gale/Hoffman
      cut  sum_{w in S} A_{C_b}(w) <= sum_i min(N_i, sum_{w in S} A_{E_{p_i}}(w)).

A_{E_j} is computed exactly by enumerating the coset row_j + R_{j+1}, which has
2^{n-1-j} elements -- cheap for LARGE j (exactly where (C) is sharp).  Strata
with n-1-j > --cap are treated as unconstrained (capacity +inf), keeping every
condition sound.
"""
from __future__ import annotations

import argparse
import collections
import itertools
import os
from typing import Dict, List, Optional, Set

import numpy as np

from cpsat_pivot_reconstruct_duality import gf2, polar_matrix


# ------------------------------------------------------------------
def code_weight_enumerator(G: np.ndarray) -> Dict[int, int]:
    """Exact A_C(w) by 2^k enumeration (Gray code)."""
    G = gf2(G); k, n = G.shape
    rows = [int("".join(map(str, r)), 2) for r in G]
    c = collections.Counter({0: 1}); cur = 0
    for msk in range(1, 1 << k):
        cur ^= rows[(msk & -msk).bit_length() - 1]
        c[bin(cur).count("1")] += 1
    return dict(c)


def stratum_weight_enumerator(Gp: np.ndarray, j: int, cap: int = 22
                              ) -> Optional[Dict[int, int]]:
    """
    A_{E_j}(w) for E_j = row_j(Gp) + span{row_{j+1}..row_{n-1}}.
    |E_j| = 2^(n-1-j).  Returns None when that exceeds 2^cap (treat as +inf).
    """
    Gp = gf2(Gp); n = Gp.shape[0]
    d = n - 1 - j
    if d > cap:
        return None
    base = int("".join(map(str, Gp[j])), 2)
    rows = [int("".join(map(str, Gp[i])), 2) for i in range(j + 1, n)]
    c = collections.Counter(); cur = base
    c[bin(cur).count("1")] += 1
    for msk in range(1, 1 << d):
        cur ^= rows[(msk & -msk).bit_length() - 1]
        c[bin(cur).count("1")] += 1
    return dict(c)


# ------------------------------------------------------------------
class WeightStratumMatch:
    """Pre-computes the pi-free stratum data for one (G_b, G_p) pair."""

    def __init__(self, Gb: np.ndarray, Gp: np.ndarray, cap: int = 22,
                 verbose: bool = False):
        self.Gb, self.Gp = gf2(Gb), gf2(Gp)
        self.k, self.n = self.Gb.shape
        self.cap = cap
        self.A_C = code_weight_enumerator(self.Gb)
        self.supp_C: Set[int] = {w for w, a in self.A_C.items() if w > 0 and a > 0}
        self.A_E: Dict[int, Optional[Dict[int, int]]] = {
            j: stratum_weight_enumerator(self.Gp, j, cap) for j in range(self.n)
        }
        # (F) forbidden indices: stratum offers no weight the code has
        self.forbidden: Set[int] = {
            j for j, A in self.A_E.items()
            if A is not None and not (set(A) & self.supp_C)
        }
        # usable capacity of each stratum = # vectors whose weight the code has
        self.cap_ok: Dict[int, Optional[int]] = {
            j: (None if A is None else sum(a for w, a in A.items() if w in self.supp_C))
            for j, A in self.A_E.items()
        }
        if verbose:
            enum = sum(1 for A in self.A_E.values() if A is not None)
            print(f"[wsm] n={self.n} k={self.k}  supp(A_C)={sorted(self.supp_C)}")
            print(f"[wsm] {enum}/{self.n} strata enumerated (cap 2^{cap})")
            print(f"[wsm] forbidden indices ({len(self.forbidden)}): "
                  f"{sorted(self.forbidden)}")

    # --------------------------------------------------------------
    def violation(self, p_star: List[int]) -> Optional[str]:
        """None => passes all conditions.  Else a reason string (unreachable)."""
        p = sorted(int(x) for x in p_star)
        k = len(p)

        # (F) forbidden index used as a pivot
        bad = [j for j in p if j in self.forbidden]
        if bad:
            return f"(F) forbidden index {bad[0]} in p_star (stratum has no weight of C_b)"

        # (C) per-stratum capacity:  E_{p_i} must host N_i = 2^{k-1-i} codewords
        for i, j in enumerate(p):
            N = 1 << (k - 1 - i)
            capj = self.cap_ok[j]
            if capj is not None and capj < N:
                return (f"(C) stratum E_{j} holds only {capj} vectors of a weight in "
                        f"supp(A_C) but must host N={N} codewords")

        # (T) Gale/Hoffman cuts over weight subsets
        for r in (1, 2):
            for S in itertools.combinations(sorted(self.supp_C), r):
                need = sum(self.A_C[w] for w in S)
                have = 0
                for i, j in enumerate(p):
                    N = 1 << (k - 1 - i)
                    A = self.A_E[j]
                    have += N if A is None else min(N, sum(A.get(w, 0) for w in S))
                if need > have:
                    return (f"(T) weights {S}: code needs {need} codewords but the "
                            f"chosen strata can host at most {have}")
        return None


# ------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir,
        "project", "POD", "eBCH_m6_t11.matrix"))
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--cap", type=int, default=22)
    args = ap.parse_args()

    def parse_mat(p):
        t = open(p).read().split(); r, c = int(t[0]), int(t[1]); body = t[2:]
        if len(body) == r:
            return np.array([[int(x) for x in row] for row in body], np.uint8)
        return np.array([int(x) for x in body], np.uint8).reshape(r, c)

    Gb = gf2(parse_mat(args.gen)); Gp = gf2(polar_matrix(args.m))
    W = WeightStratumMatch(Gb, Gp, cap=args.cap, verbose=True)
    print(f"[wsm] A_C = {dict(sorted(W.A_C.items()))}")
    for j in sorted(W.A_E):
        A = W.A_E[j]
        if A is not None:
            print(f"  E_{j:2d}: |E|=2^{W.n-1-j:<2d}  A_E={dict(sorted(A.items()))}"
                  f"   usable={W.cap_ok[j]}"
                  f"{'   <-- FORBIDDEN' if j in W.forbidden else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
