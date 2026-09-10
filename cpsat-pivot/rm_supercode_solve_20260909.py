#!/usr/bin/env python3
"""
rm_supercode_solve  --  STAGE 1: find a permutation P_0 realizing a desired
RM degree bound  d(P_0) <= d, i.e.

        C_b P_0  subseteq  RM(d, m).

Convention (matches FORMULATIONS.md): row rho_a of F^{(x m)} has weight 2^|a|,
hence monomial degree m - |a|.  With  W_0 = G_b P_0 F^{(x m)},

        d(P_0) = max{ m - wt_2(j) : W_0[:,j] != 0 }.

So the containment is exactly the vanishing of every "too-high-degree" column:

        W_0[:,j] = 0    for all j with  wt_2(j) < m - d.

Expanding  W_0[r,j] = XOR_{a in supp(Gp[:,j])} G_b[r, pi(a)]  and noting
supp(Gp[:,j]) = supermasks of j (size 2^(m-wt_2(j))), each such column is one
parity per generator row.  The whole model is therefore

        pi in Sym(2^m)                        (AllDifferent)
        XOR_{a : a superset of j} b[r,a] = 0   for wt_2(j) < m-d, all r
        b[r,a] = G_b[r, pi(a)]                 (AddElement)

PURE permutation + XOR.  None of the pivot-profile machinery is present: no
U-certificate (k^3 AND-gates), no span-membership.  Those are what made the
profile problem hard; this is a strictly smaller problem.

Constraint counts for m=6, k=16 (the j=0 column is dropped -- it reduces to
"every generator row has even weight", a tautology for extended BCH):

        d=2 : 41 columns x 16 rows = 656 parities, widths 8..32
        d=3 : 21 columns x 16 rows = 336 parities, widths 16..32
        d=4 :  6 columns x 16 rows =  96 parities, width 32

Feasibility floor:  dim RM(d,m) >= dim C_b.
        [64,16] -> d>=2      [64,36] -> d>=3      [64,51] -> d>=4

A larger d is a strictly WEAKER constraint (RM(d,m) subset RM(d+1,m)), so the
solution set grows with d.  It buys no new reliable positions -- for [64,16]
the 16 most reliable indices all have wt_2 >= 4 already, so they sit inside
RM(2,6) -- but it admits permutations whose NON-pivot support needs degree
> d_min, which the tighter bound would discard.  pi_eq is exactly such a case:
d(pi_eq)=5 yet its pivot profile equals that of the d=2 field permutation.
"""
from __future__ import annotations

import argparse
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from ortools.sat.python import cp_model

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, polar_matrix, column_pivot_profile_gf2)


# ------------------------------------------------------------------
def forced_zero_columns(m: int, d: int, drop_trivial: bool = True) -> List[int]:
    """Columns j whose monomial degree m-wt_2(j) exceeds d."""
    n = 1 << m
    js = [j for j in range(n) if bin(j).count("1") < m - d]
    if drop_trivial and 0 in js:
        js.remove(0)          # XOR over ALL a == row parity == 0 for even codes
    return js


def rm_degree(Gb: np.ndarray, Gp: np.ndarray, pi: Sequence[int], m: int) -> int:
    W = gf2_matmul(gf2(Gb)[:, list(pi)], gf2(Gp))
    nz = [j for j in range(W.shape[1]) if W[:, j].any()]
    return max(m - bin(j).count("1") for j in nz)


def build(Gb: np.ndarray, Gp: np.ndarray, m: int, d: int,
         forbid: Sequence[int] = (), fix_pi: Optional[Sequence[int]] = None,
         verbose: bool = True):
    Gb, Gp = gf2(Gb), gf2(Gp)
    k, n = Gb.shape
    model = cp_model.CpModel()

    pi = [model.NewIntVar(0, n - 1, f"pi[{a}]") for a in range(n)]
    model.AddAllDifferent(pi)
    if fix_pi is not None:
        for a, v in enumerate(fix_pi):
            model.Add(pi[a] == int(v))

    # b[r,a] = G_b[r, pi(a)]
    b: Dict[Tuple[int, int], cp_model.IntVar] = {}
    for r in range(k):
        row = [int(x) for x in Gb[r]]
        for a in range(n):
            b[r, a] = model.NewBoolVar(f"b[{r},{a}]")
            model.AddElement(pi[a], row, b[r, a])
        model.Add(sum(b[r, a] for a in range(n)) == int(Gb[r].sum()))   # weight cut

    cols = forced_zero_columns(m, d)
    cols = sorted(set(cols) | set(int(j) for j in forbid))
    npar = 0
    for j in cols:
        supp = [a for a in range(n) if int(Gp[a, j]) == 1]
        if not supp:
            continue
        for r in range(k):
            _add_xor_zero(model, [b[r, a] for a in supp])
            npar += 1
    if verbose:
        widths = sorted({len([a for a in range(n) if int(Gp[a, j]) == 1]) for j in cols})
        print(f"[rm] m={m} k={k} d={d}: {len(cols)} forced-zero columns "
              f"(widths {widths}) -> {npar} parities")
    return model, {"pi": pi, "cols": cols, "parities": npar}


def _add_xor_zero(model, lits) -> None:
    """XOR of lits == 0, via an integer parity variable (no AND gates)."""
    t = model.NewIntVar(0, len(lits) // 2, "")
    model.Add(sum(lits) == 2 * t)


def solve(Gb, Gp, m, d, *, forbid=(), time_limit=300.0, workers=8,
          verbose=True):
    Gb, Gp = gf2(Gb), gf2(Gp)
    n = Gb.shape[1]
    model, v = build(Gb, Gp, m, d, forbid=forbid, verbose=verbose)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(time_limit)
    s.parameters.num_search_workers = int(workers)
    t0 = time.time()
    st = s.Solve(model)
    dt = time.time() - t0
    name = s.StatusName(st)
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        sol = [int(s.Value(v["pi"][a])) for a in range(n)]
        return "sat", sol, dt
    return ("infeasible" if name == "INFEASIBLE" else "unknown"), None, dt


def report(Gb, Gp, m, pi, Z=None) -> dict:
    Gb, Gp = gf2(Gb), gf2(Gp)
    W = gf2_matmul(Gb[:, list(pi)], Gp)
    prof = column_pivot_profile_gf2(W)
    nz = [j for j in range(W.shape[1]) if W[:, j].any()]
    out = {"profile": prof, "d": max(m - bin(j).count("1") for j in nz),
           "nonzero_cols": len(nz)}
    if Z is not None:
        out["zsum"] = sum(Z[i] for i in prof)
    return out


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
    ap.add_argument("--forbid", default="", help="extra columns forced to zero")
    ap.add_argument("--time", type=float, default=300.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--snr", type=float, default=3.0)
    ap.add_argument("--out", default=None, help="write realizing P_0 .matrix here")
    args = ap.parse_args()

    def parse_any(p):
        tok = open(p).read().split(); r, c = int(tok[0]), int(tok[1]); b = tok[2:]
        return (np.array([int(x) for x in b], np.uint8).reshape(r, c)
                if len(b) == r * c else
                np.array([[int(ch) for ch in row] for row in b], np.uint8))

    Gb = gf2(parse_any(args.gen)); Gp = gf2(polar_matrix(args.m))
    k, n = Gb.shape
    Z = bha_seq(args.m, args.snr, k / n)
    forbid = [int(x) for x in args.forbid.split(",") if x.strip()]

    st, pi, dt = solve(Gb, Gp, args.m, args.d, forbid=forbid,
                       time_limit=args.time, workers=args.workers)
    print(f"\n[rm] status: {st}   ({dt:.2f}s)")
    if st != "sat":
        return 1
    rep = report(Gb, Gp, args.m, pi, Z)
    ok = rep["d"] <= args.d
    print(f"[rm] realizing P_0 found.  verified d(P_0) = {rep['d']}  "
          f"(<= {args.d}: {ok})")
    print(f"[rm] nonzero coefficient columns: {rep['nonzero_cols']}")
    print(f"[rm] pivot profile : {rep['profile']}")
    print(f"[rm] Z-sum @ {args.snr}dB : {rep['zsum']:.6f}")
    if args.out:
        P = np.zeros((n, n), np.uint8)
        for a, v in enumerate(pi):
            P[a, v] = 1
        with open(args.out, "w") as f:
            f.write(f"{n} {n}\n")
            for row in P:
                f.write("".join(map(str, row)) + "\n")
        print(f"[rm] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
