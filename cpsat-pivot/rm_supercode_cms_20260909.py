#!/usr/bin/env python3
"""
rm_supercode_cms  --  STAGE 1 for a GENERAL [2^m, k] code, via CryptoMiniSat.

Problem: find a coordinate permutation P_0 with

        C_b P_0  subseteq  RM(d, m)          i.e.   d(P_0) <= d.

For eBCH the field-expansion permutation solves this in closed form; for an
arbitrary generator matrix there is no such shortcut and we must search.  The
encoding here is PURE permutation + XOR, which is exactly CryptoMiniSat's
Gaussian-elimination regime (unlike the pivot-profile problem, whose U W = I
certificate is bilinear and defeats the XOR engine).

Encoding
--------
one-hot permutation      x[a][v] = 1  iff  pi(a) = v
    exactly-one over v for each a        (ALO clause + pairwise/ladder AMO)
    exactly-one over a for each v
    + native XOR parity row sum_v x[a][v] = 1 and sum_a x[a][v] = 1
      (redundant given ALO+AMO, but puts the constraint in the Gauss matrix)

permuted-code bits       b[r][a] = G_b[r, pi(a)] = XOR_{v : G_b[r,v]=1} x[a][v]
    native XOR:  b[r][a]  XOR  ( XOR_{v in supp(G_b[r])} x[a][v] )  = 0

degree bound             for every column j with  m - wt_2(j) > d  (j != 0):
    supp(Gp[:,j]) = { a : a superset of j }        (size 2^(m - wt_2(j)))
    native XOR:  XOR_{a in supp(Gp[:,j])} b[r][a] = 0    for every row r

j = 0 gives  XOR_a b[r][a] = wt(G_b[r]) mod 2, dropped when every row is even
(reported otherwise).

--enumerate K   add a blocking clause on the permutation after each solution
                and keep going, to seed the stage-2 GL(m,2) walk.
"""
from __future__ import annotations

import argparse
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pycryptosat

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, polar_matrix, column_pivot_profile_gf2)


def _parse_any(p: str) -> np.ndarray:
    tok = open(p).read().split(); r, c = int(tok[0]), int(tok[1]); b = tok[2:]
    if len(b) == r * c:
        return np.array([int(x) for x in b], np.uint8).reshape(r, c)
    return np.array([[int(ch) for ch in row] for row in b], np.uint8)


class Builder:
    def __init__(self) -> None:
        self.nv = 0
        self.clauses: List[List[int]] = []
        self.xors: List[Tuple[List[int], bool]] = []

    def nvar(self) -> int:
        self.nv += 1
        return self.nv

    def clause(self, lits: Sequence[int]) -> None:
        self.clauses.append(list(lits))

    def xor(self, lits: Sequence[int], rhs: bool) -> None:
        self.xors.append((list(lits), bool(rhs)))

    def exactly_one(self, lits: Sequence[int]) -> None:
        lits = list(lits)
        self.clause(lits)                                  # ALO
        for i in range(len(lits)):                         # AMO pairwise
            for j in range(i + 1, len(lits)):
                self.clause([-lits[i], -lits[j]])
        self.xor(lits, True)                               # parity into Gauss

    def feed(self, solver: "pycryptosat.Solver") -> None:
        for c in self.clauses:
            solver.add_clause(c)
        for vs, rhs in self.xors:
            solver.add_xor_clause(vs, rhs)


def build(Gb: np.ndarray, m: int, d: int,
          forbid_cols: Sequence[int] = ()) -> Tuple[Builder, np.ndarray, dict]:
    Gb = gf2(Gb)
    k, n = Gb.shape
    assert n == (1 << m), f"n={n} != 2^{m}"
    Gp = gf2(polar_matrix(m))

    B = Builder()
    # one-hot permutation matrix x[a][v]
    x = [[B.nvar() for _ in range(n)] for _ in range(n)]
    for a in range(n):
        B.exactly_one([x[a][v] for v in range(n)])
    for v in range(n):
        B.exactly_one([x[a][v] for a in range(n)])

    # b[r][a] = XOR_{v in supp(G_b[r])} x[a][v]
    supp_row = [[v for v in range(n) if Gb[r, v]] for r in range(k)]
    b = [[B.nvar() for _ in range(n)] for _ in range(k)]
    for r in range(k):
        for a in range(n):
            B.xor([b[r][a]] + [x[a][v] for v in supp_row[r]], False)

    # degree bound: kill every column of degree > d
    cols = [j for j in range(n) if (m - bin(j).count("1")) > d]
    dropped0 = False
    if 0 in cols:
        cols.remove(0); dropped0 = True
    cols = sorted(set(cols) | set(int(j) for j in forbid_cols))

    even = bool(np.all(Gb.sum(1) % 2 == 0))
    npar = 0
    for j in cols:
        supp_a = [a for a in range(n) if Gp[a, j]]
        for r in range(k):
            B.xor([b[r][aa] for aa in supp_a], False)
            npar += 1
    if dropped0 and not even:
        # j=0 is a real constraint: XOR_a b[r][a] = wt(row r) mod 2
        for r in range(k):
            B.xor([b[r][a] for a in range(n)], bool(Gb[r].sum() % 2))
            npar += 1

    meta = {"k": k, "n": n, "cols": cols, "parities": npar,
            "x_base": x[0][0], "vars": B.nv, "clauses": len(B.clauses),
            "xors": len(B.xors), "even_code": even}
    return B, np.array(x), meta


def _extract_pi(model: Sequence[bool], x: np.ndarray, n: int) -> List[int]:
    pi = [-1] * n
    for a in range(n):
        for v in range(n):
            if model[x[a][v]]:
                pi[a] = v
                break
    return pi


def rm_degree(Gb, Gp, pi, m) -> int:
    W = gf2_matmul(gf2(Gb)[:, list(pi)], gf2(Gp))
    nz = [j for j in range(W.shape[1]) if W[:, j].any()]
    return max(m - bin(j).count("1") for j in nz)


def solve(Gb: np.ndarray, m: int, d: int, *, forbid_cols=(), enumerate_k=1,
          threads=8, max_seconds: Optional[float] = None, verbose=True):
    Gb = gf2(Gb); Gp = gf2(polar_matrix(m))
    k, n = Gb.shape
    B, x, meta = build(Gb, m, d, forbid_cols)
    if verbose:
        print(f"[cms-rm] n={n} k={k} d={d}: {len(meta['cols'])} killed columns, "
              f"{meta['parities']} degree parities")
        print(f"[cms-rm] vars={meta['vars']} clauses={meta['clauses']} "
              f"xors={meta['xors']}  even_code={meta['even_code']}")

    solver = pycryptosat.Solver(threads=threads)
    B.feed(solver)

    found = []
    t0 = time.time()
    for it in range(enumerate_k):
        sat, model = solver.solve()
        dt = time.time() - t0
        if not sat:
            if verbose:
                print(f"[cms-rm] {'UNSAT' if it == 0 else 'no more solutions'} "
                      f"after {it} found  ({dt:.1f}s)")
            break
        pi = _extract_pi(model, x, n)
        dd = rm_degree(Gb, Gp, pi, m)
        prof = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], Gp))
        found.append({"pi": pi, "d": dd, "profile": prof, "seconds": round(dt, 1)})
        if verbose:
            print(f"[cms-rm] solution {it+1}: verified d={dd} "
                  f"(<= {d}: {dd <= d})  profile={prof}  ({dt:.1f}s)")
        # block this exact permutation
        solver.add_clause([-int(x[a][pi[a]]) for a in range(n)])
        if max_seconds and dt > max_seconds:
            break
    return found, meta


# ------------------------------------------------------------------
def main() -> int:
    import os
    from stratum_sweep_gen_20260904 import bha_seq
    HERE = os.path.dirname(os.path.abspath(__file__))
    POD = os.path.join(HERE, os.pardir, "project", "POD")

    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default=os.path.join(POD, "eBCH_m6_t11.matrix"))
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--forbid", default="")
    ap.add_argument("--enumerate", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-seconds", type=float, default=None)
    ap.add_argument("--snr", type=float, default=3.0)
    ap.add_argument("--out-prefix", default=None)
    args = ap.parse_args()

    Gb = gf2(_parse_any(args.gen)); k, n = Gb.shape
    Z = bha_seq(args.m, args.snr, k / n)
    forbid = [int(x) for x in args.forbid.split(",") if x.strip()]

    found, meta = solve(Gb, args.m, args.d, forbid_cols=forbid,
                        enumerate_k=args.enumerate, threads=args.threads,
                        max_seconds=args.max_seconds)
    if not found:
        print("[cms-rm] no realizing P_0.")
        return 1
    for i, sol in enumerate(found):
        zs = sum(Z[j] for j in sol["profile"])
        print(f"  #{i+1}: d={sol['d']}  Z-sum@{args.snr}dB={zs:.6f}  "
              f"profile={sol['profile']}")
        if args.out_prefix:
            P = np.zeros((n, n), np.uint8)
            for a, v in enumerate(sol["pi"]):
                P[a, v] = 1
            fn = f"{args.out_prefix}_{i+1}.matrix"
            with open(fn, "w") as f:
                f.write(f"{n} {n}\n")
                for row in P:
                    f.write("".join(map(str, row)) + "\n")
            print(f"       wrote {fn}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
