#!/usr/bin/env python3
"""
gl_shear_search  --  STAGE 1 (field expansion) + STAGE 2 (GL(m,2) shear walk)
for the inverse pivot-profile / reliability problem on extended primitive BCH.

    STAGE 1   P_field : cyclic -> monomial basis.  Fixes the RM degree bound
              d(P_0) = max{ m - wt_2(j) : W_0[:,j] != 0 },   W_0 = G_b P_0 F^(x m).
    STAGE 2   walk the Cayley graph of GL(m,2) with transvection generators.
              Q_A is an RM automorphism, so d is INVARIANT along the walk;
              only the pivot reliability changes.

--------------------------------------------------------------------------
TWO CONVENTION TRAPS (both cost real time on 2026-09; do not re-derive them)
--------------------------------------------------------------------------
(1) DIRECTION of the field expansion.  With the house convention
    pi = [argmax(P[a,:])] and the code relabelled as G_b[:, pi], the correct
    map is the DISCRETE LOG

        pi_field(x) = 1 + log_alpha(x)   for x != 0,     pi_field(0) = 0
        (LSB of x = coefficient of alpha^0)

    NOT  pi(j) = lambda(beta_j)  with beta_j = alpha^(j-1), which is its
    inverse.  Both use the same primitive polynomial and the same bit order;
    only the direction differs.  Getting it backwards is silent: at m=6 the
    inverse gives d=5 and at m=7 it gives d=6 -- exactly the degree of the
    repo's P_eq file -- so the result looks like "field expansion doesn't
    help" instead of like a bug.

(2) EQUAL PIVOT PROFILE != EQUAL DEGREE.  P_field and the repo's P_eq have
    the IDENTICAL pivot profile (hence identical Z-sum) at both m=6 and m=7,
    but different coefficient support and different d:

        m=6:  P_eq d=5,  P_field d=2   (both Z-sum 2.316021)
        m=7:  P_eq d=6,  P_field d=4   (both Z-sum 9.124280)

    So the profile does not identify the starting point.  Comparing d across
    two permutations only tells you about the walk if they share a base.

--------------------------------------------------------------------------
COMPOSITION ORDER
--------------------------------------------------------------------------
For A in GL(m,2) with Q_A(z) = A z, the RM-automorphism action on a relabelled
code is RIGHT composition:

        pi'  =  pi o Q_A        i.e.   pi'[a] = pi[Q_A[a]]        (d preserved)
        pi'  =  Q_A o pi                                          (d NOT preserved)

Verified 300/300 random A at m=6.

--------------------------------------------------------------------------
FEASIBILITY FLOOR (dim RM(d,m) >= k, plus d_min compatibility)
--------------------------------------------------------------------------
    eBCH[64,16]   d >= 2   dim RM(2,6)=22
    eBCH[64,36]   d >= 3   dim RM(3,6)=42
    eBCH[128,64]  d >= 4   dim RM(3,7)=64 == k would force C_b = RM(3,7), but
                           d_min 22 != 16, so d=3 is impossible; RM(4,7)=99.
"""
from __future__ import annotations

import argparse, math, random, time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from cpsat_pivot_reconstruct_duality import (
    gf2, gf2_matmul, polar_matrix, column_pivot_profile_gf2)

# Conway polynomials for GF(2^m), as reduction masks (x^m implicit)
CONWAY = {3: 0b011, 4: 0b0011, 5: 0b00101, 6: 0b011011,
          7: 0b0000011, 8: 0b00011101, 9: 0b000010001, 10: 0b0000001001}


def parse_mat(p: str) -> np.ndarray:
    t = open(p).read().split(); r, c = int(t[0]), int(t[1]); b = t[2:]
    if len(b) == r * c:
        return np.array([int(x) for x in b], np.uint8).reshape(r, c)
    return np.array([[int(ch) for ch in row] for row in b], np.uint8)


def bha_seq(m: int, ebn0_db: float, rate: float) -> List[float]:
    z = [math.exp(-rate * 10 ** (ebn0_db / 10))]
    for _ in range(m):
        z = [v for x in z for v in (2 * x - x * x, x * x)]
    return z


# ---------------------------------------------------------------- stage 1
def field_permutation(m: int, red: Optional[int] = None) -> List[int]:
    """pi_field(x) = 1 + log_alpha(x), pi_field(0)=0.  See trap (1)."""
    if red is None:
        red = CONWAY[m]
    n = 1 << m
    pw = [0] * (n - 1); e = 1
    for i in range(n - 1):
        pw[i] = e
        e <<= 1
        if (e >> m) & 1:
            e = (e ^ (1 << m)) ^ red
    if e != 1 or len(set(pw)) != n - 1:
        raise ValueError(f"reduction mask {red:#b} is not primitive for m={m}")
    log = {pw[i]: i for i in range(n - 1)}
    return [0] + [1 + log[x] for x in range(1, n)]


def best_field_permutation(Gb, Gp, m) -> Tuple[List[int], int, int]:
    """Try every primitive polynomial; return the one minimising d."""
    n = 1 << m
    best = None
    for red in range(1, 1 << m):
        try:
            pi = field_permutation(m, red)
        except ValueError:
            continue
        d, _ = rm_degree(Gb, Gp, pi, m)
        if best is None or d < best[1]:
            best = (pi, d, red)
    return best


# ---------------------------------------------------------------- geometry
def rm_degree(Gb, Gp, pi, m) -> Tuple[int, int]:
    W = gf2_matmul(gf2(Gb)[:, list(pi)], gf2(Gp))
    nz = [j for j in range(W.shape[1]) if W[:, j].any()]
    return max(m - bin(j).count("1") for j in nz), len(nz)


def transvections(m: int) -> List[Tuple[int, int]]:
    """T = I + E_{i,j}, i != j.  m(m-1) of them; they generate GL(m,2)."""
    return [(i, j) for i in range(m) for j in range(m) if i != j]


def q_of(m: int, i: int, j: int) -> List[int]:
    """Coordinate permutation Q_A of [0,2^m) for A = I + E_{i,j}: z_i <- z_i + z_j."""
    n = 1 << m
    return [z ^ (1 << i) if (z >> j) & 1 else z for z in range(n)]


def apply_q(pi: Sequence[int], q: Sequence[int]) -> List[int]:
    """pi o Q_A  -- RIGHT composition.  See COMPOSITION ORDER above."""
    return [pi[q[a]] for a in range(len(pi))]


# ---------------------------------------------------------------- objective
class Objective:
    def __init__(self, Gb, Gp, m, snr=3.0):
        self.Gb, self.Gp, self.m = gf2(Gb), gf2(Gp), m
        k, n = self.Gb.shape
        self.Z = bha_seq(m, snr, k / n)
        self.n, self.k = n, k
        self.cache: Dict[Tuple[int, ...], Tuple[float, List[int]]] = {}

    def __call__(self, pi: Sequence[int]) -> Tuple[float, List[int]]:
        key = tuple(pi)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        W = gf2_matmul(self.Gb[:, list(pi)], self.Gp)
        prof = column_pivot_profile_gf2(W)
        val = sum(self.Z[i] for i in prof)
        if len(self.cache) < 400000:
            self.cache[key] = (val, prof)
        return val, prof

    def floor(self) -> float:
        return sum(sorted(self.Z)[:self.k])


# ---------------------------------------------------------------- stage 2
def greedy_walk(obj: Objective, pi0: Sequence[int], m: int,
                max_steps=200, verbose=True):
    """Steepest descent over transvection neighbours."""
    gens = [(i, j, q_of(m, i, j)) for i, j in transvections(m)]
    pi = list(pi0)
    cur, prof = obj(pi)
    path = []
    for step in range(max_steps):
        best = None
        for i, j, q in gens:
            cand = apply_q(pi, q)
            v, pr = obj(cand)
            if best is None or v < best[0]:
                best = (v, cand, (i, j), pr)
        if best[0] >= cur - 1e-12:
            break
        cur, pi, prof = best[0], best[1], best[3]
        path.append((best[2], cur))
        if verbose:
            print(f"  step {step+1:3d}: z{best[2][0]+1}<-z{best[2][0]+1}+z{best[2][1]+1}"
                  f"   Z-sum {cur:.6f}", flush=True)
    return pi, cur, prof, path


def beam_walk(obj: Objective, pi0: Sequence[int], m: int, width=16,
              max_steps=120, verbose=True):
    """Beam search: keep `width` best states, expand all transvections."""
    gens = [(i, j, q_of(m, i, j)) for i, j in transvections(m)]
    v0, p0 = obj(pi0)
    beam = [(v0, list(pi0))]
    seen = {tuple(pi0)}
    best = (v0, list(pi0), p0)
    for step in range(max_steps):
        cands = []
        for val, pi in beam:
            for i, j, q in gens:
                c = apply_q(pi, q)
                t = tuple(c)
                if t in seen:
                    continue
                seen.add(t)
                v, pr = obj(c)
                cands.append((v, c, pr))
        if not cands:
            break
        cands.sort(key=lambda x: x[0])
        beam = [(v, c) for v, c, _ in cands[:width]]
        if cands[0][0] < best[0] - 1e-12:
            best = (cands[0][0], cands[0][1], cands[0][2])
            if verbose:
                print(f"  beam step {step+1:3d}: Z-sum {best[0]:.6f}  "
                      f"(seen {len(seen)})", flush=True)
        elif verbose and (step + 1) % 20 == 0:
            print(f"  beam step {step+1:3d}: no improvement, best {best[0]:.6f}",
                  flush=True)
    return best[1], best[0], best[2]


def random_restart(obj: Objective, pi0, m, restarts=20, kick=6, seed=0,
                   verbose=True):
    """Greedy from pi0, then repeated random-transvection kicks + greedy."""
    rng = random.Random(seed)
    gens = [q_of(m, i, j) for i, j in transvections(m)]
    pi, val, prof, _ = greedy_walk(obj, pi0, m, verbose=False)
    best = (val, pi, prof)
    if verbose:
        print(f"  restart  0: {val:.6f}")
    for r in range(restarts):
        cur = list(best[1])
        for _ in range(kick):
            cur = apply_q(cur, rng.choice(gens))
        pi, val, prof, _ = greedy_walk(obj, cur, m, verbose=False)
        if val < best[0] - 1e-12:
            best = (val, pi, prof)
        if verbose:
            print(f"  restart {r+1:2d}: {val:.6f}   best {best[0]:.6f}", flush=True)
    return best[1], best[0], best[2]


# ---------------------------------------------------------------- driver
def main() -> int:
    import os
    HERE = os.path.dirname(os.path.abspath(__file__))
    POD = os.path.join(HERE, os.pardir, "project", "POD")
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", default=os.path.join(POD, "eBCH_m6_t11.matrix"))
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--snr", type=float, default=3.0)
    ap.add_argument("--mode", choices=["greedy", "beam", "restart"], default="beam")
    ap.add_argument("--width", type=int, default=16)
    ap.add_argument("--steps", type=int, default=120)
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--kick", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--scan-poly", action="store_true",
                    help="try all primitive polys, keep min-d")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    m, n = args.m, 1 << args.m
    Gb = gf2(parse_mat(args.gen)); Gp = gf2(polar_matrix(m))
    k = Gb.shape[0]
    obj = Objective(Gb, Gp, m, args.snr)

    if args.scan_poly:
        pi0, d0, red = best_field_permutation(Gb, Gp, m)
        print(f"[stage1] scanned primitive polys -> mask {red:#0{m+2}b}, d={d0}")
    else:
        red = CONWAY[m]
        pi0 = field_permutation(m, red)
        d0, nnz = rm_degree(Gb, Gp, pi0, m)
        print(f"[stage1] Conway mask {red:#0{m+2}b}: d(P_field) = {d0}, "
              f"{nnz} nonzero columns")
    v0, p0 = obj(pi0)
    print(f"[stage1] Z-sum @ {args.snr}dB = {v0:.6f}")
    print(f"[stage1] free-J floor (k most reliable) = {obj.floor():.6f}")
    print(f"[stage2] |GL({m},2)| = {math.prod((1<<m)-(1<<i) for i in range(m)):,}"
          f"   {len(transvections(m))} transvection generators")

    t0 = time.time()
    if args.mode == "greedy":
        pi, val, prof, _ = greedy_walk(obj, pi0, m)
    elif args.mode == "beam":
        pi, val, prof = beam_walk(obj, pi0, m, args.width, args.steps)
    else:
        pi, val, prof = random_restart(obj, pi0, m, args.restarts, args.kick,
                                       args.seed)
    dt = time.time() - t0

    d1, nnz1 = rm_degree(Gb, Gp, pi, m)
    print(f"\n[stage2] {args.mode}: {v0:.6f} -> {val:.6f}   "
          f"(delta {val-v0:+.6f})   {dt:.1f}s")
    print(f"[stage2] d preserved: {d0} -> {d1}  {'OK' if d0==d1 else 'VIOLATION'}")
    print(f"[stage2] profile = {prof}")
    if args.out:
        P = np.zeros((n, n), np.uint8)
        for a, v in enumerate(pi):
            P[a, v] = 1
        with open(args.out, "w") as f:
            f.write(f"{n} {n}\n")
            for row in P:
                f.write("".join(map(str, row)) + "\n")
        print(f"[stage2] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
