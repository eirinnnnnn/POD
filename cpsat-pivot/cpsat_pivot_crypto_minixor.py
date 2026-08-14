#!/usr/bin/env python3
"""
CryptoMiniSat (CNF + native XOR) reconstruction of a permutation pi with

    PivotProfile( pi(G_b) F^{tensor m} ) = p_star,

convention pi(G_b) = G_b[:, pi], arithmetic over F_2.

Why CryptoMiniSat
-----------------
The problem is an F_2-linear system (permutation -> b -> butterfly -> W ->
identity/span parities) coupled to a small multiplicative certificate
layer.  CP-SAT has no parity reasoning: XOR constraints only propagate
when a single literal remains free.  CryptoMiniSat treats XOR clauses
natively and runs Gauss-Jordan elimination on the parity matrix DURING
search, which is exactly the missing inference.  Everything linear below
is emitted as XOR clauses; only the AND gates are CNF.

Encoding
--------
Permutation:  P[a][t] Boolean, "position a takes source column t".
    - AMO per row and per column (sequential/Sinz encoding),
    - ALO clause per row and per column,
    - native XOR parity  XOR_t P[a][t] = 1  per row and per column
      (redundant with ALO+AMO, but it places the exactly-one structure
      inside the Gauss matrix where it propagates with the rest of the
      linear system).

Channeling:  b[r,a] = G_b[r, pi(a)] = XOR_{t in supp(G_b[r,:])} P[a][t]
    -- a single XOR clause per (r,a).  This identity holds because
    exactly one P[a][t] is true, so the parity over the support equals
    the selected entry.

W:  either the polar butterfly recursion (3-literal XOR clauses, sparse
    parity matrix; default) or dense per-column XORs over the support
    of Gp[:,j] (fewer, longer parity rows -- also fine for Gauss).

Certificate (pivot-only; equivalent to the full profile condition, see
the CP-SAT v2 header for the proof):
    (a) U W[:,p_star] = I_k :
            z[i,r,l] = u[i,r] AND w[r,p_l]        (Tseitin CNF)
            XOR_r z[i,r,l] = delta_{i,l}           (XOR clause)
    (b) span membership for non-pivot j <= p_star[-1]:
            y[l,r,j] = lambda[l,j] AND w[r,p_l]    (Tseitin CNF)
            w[r,j] = XOR_l y[l,r,j]                (XOR clause)
        with the zero-previous-pivots case giving unit clauses w[r,j]=0.

Matroid-duality reduction (--dualize auto|always|never): identical to
the CP-SAT v2 version.  For standard polar Gp,

    PivotProfile(pi(Gb) F) = p_star
        <=>  PivotProfile(sigma(Hb) F) = q_star,

with sigma = pi[::-1], Hb a nullspace basis of Gb, and
q_star = sorted{ n-1-j : j not in p_star }.  Solve whichever of the two
instances has fewer rows.

UNSAT is a *proof of unreachability* of p_star (CryptoMiniSat is a
complete solver), which is often the interesting answer for the inverse
pivot-profile problem.  --enumerate N blocks solutions incrementally
and returns up to N distinct permutations realizing p_star.

Usage
-----
    python cms_pivot_reconstruct.py                       # self-tests
    python cms_pivot_reconstruct.py data.npz --time 14400 --threads 12
    python cms_pivot_reconstruct.py data.npz --enumerate 10
    python cms_pivot_reconstruct.py data.npz --dimacs inst.cnf   # export
                     # extended DIMACS ('x' lines) for cryptominisat5 CLI

Requires:  pip install pycryptosat numpy
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from pycryptosat import Solver as CMSSolver
    HAVE_PYCRYPTOSAT = True
except ImportError:
    HAVE_PYCRYPTOSAT = False


# ============================================================
# F_2 utilities (self-contained; mirror cpsat_pivot_reconstruct_v2)
# ============================================================

def gf2(A: np.ndarray) -> np.ndarray:
    return (np.asarray(A, dtype=np.uint8) & 1)


def gf2_matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = gf2(A)
    B = gf2(B)
    return ((A.astype(np.uint32) @ B.astype(np.uint32)) & 1).astype(np.uint8)


def gf2_rank(A: np.ndarray) -> int:
    A = gf2(A).copy()
    m, n = A.shape
    rank = 0
    for col in range(n):
        pivot = None
        for row in range(rank, m):
            if A[row, col]:
                pivot = row
                break
        if pivot is None:
            continue
        if pivot != rank:
            A[[rank, pivot]] = A[[pivot, rank]]
        for row in range(m):
            if row != rank and A[row, col]:
                A[row, :] ^= A[rank, :]
        rank += 1
        if rank == m:
            break
    return rank


def gf2_rref(A: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    R = gf2(A).copy()
    m, n = R.shape
    pivots: List[int] = []
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, m):
            if R[r, col]:
                piv = r
                break
        if piv is None:
            continue
        if piv != row:
            R[[row, piv]] = R[[piv, row]]
        for r in range(m):
            if r != row and R[r, col]:
                R[r, :] ^= R[row, :]
        pivots.append(col)
        row += 1
        if row == m:
            break
    return R, pivots


def gf2_nullspace(A: np.ndarray) -> np.ndarray:
    A = gf2(A)
    m, n = A.shape
    R, pivots = gf2_rref(A)
    pivot_set = set(pivots)
    free = [c for c in range(n) if c not in pivot_set]
    H = np.zeros((len(free), n), dtype=np.uint8)
    for i, f in enumerate(free):
        H[i, f] = 1
        for r, p in enumerate(pivots):
            H[i, p] = R[r, f]
    return H


def column_pivot_profile_gf2(A: np.ndarray) -> List[int]:
    A = gf2(A)
    k, n = A.shape
    basis = [0 for _ in range(k)]
    pivots: List[int] = []
    for j in range(n):
        v = 0
        for r in range(k):
            if A[r, j]:
                v |= (1 << r)
        for p in reversed(range(k)):
            if ((v >> p) & 1) and basis[p] != 0:
                v ^= basis[p]
        if v != 0:
            p = v.bit_length() - 1
            basis[p] = v
            for q in range(k):
                if q != p and basis[q] != 0 and ((basis[q] >> p) & 1):
                    basis[q] ^= v
            pivots.append(j)
    return pivots


def polar_matrix(m: int) -> np.ndarray:
    F = np.array([[1, 0], [1, 1]], dtype=np.uint8)
    G = np.array([[1]], dtype=np.uint8)
    for _ in range(m):
        G = np.kron(G, F).astype(np.uint8) & 1
    return G


def is_power_of_two(n: int) -> bool:
    return n > 0 and (n & (n - 1)) == 0


def polar_log2_size(n: int) -> int:
    if not is_power_of_two(n):
        raise ValueError(f"n must be a power of two, got n={n}.")
    return int(math.log2(n))


def is_standard_polar_matrix(Gp: np.ndarray) -> bool:
    Gp = gf2(Gp)
    if Gp.ndim != 2 or Gp.shape[0] != Gp.shape[1]:
        return False
    n = Gp.shape[0]
    if not is_power_of_two(n):
        return False
    return np.array_equal(Gp, polar_matrix(polar_log2_size(n)))


def random_full_row_rank_binary_matrix(k: int, n: int, rng: np.random.Generator) -> np.ndarray:
    while True:
        Gb = rng.integers(0, 2, size=(k, n), dtype=np.uint8)
        if gf2_rank(Gb) == k:
            return Gb


# ============================================================
# Matroid-duality reduction
# ============================================================

def dual_target_profile(n: int, p_star: List[int]) -> List[int]:
    pset = set(p_star)
    return sorted(n - 1 - j for j in range(n) if j not in pset)


def dualize_instance(Gb: np.ndarray, p_star: List[int]) -> Tuple[np.ndarray, List[int]]:
    Gb = gf2(Gb)
    k, n = Gb.shape
    Hb = gf2_nullspace(Gb)
    if Hb.shape != (n - k, n):
        raise RuntimeError("Nullspace has unexpected shape; Gb must have full row rank.")
    q_star = dual_target_profile(n, p_star)
    if len(q_star) != n - k:
        raise RuntimeError("Dual target profile has wrong length.")
    return Hb, q_star


def primal_pi_from_dual_sigma(sigma: List[int]) -> List[int]:
    return list(reversed(sigma))


# ============================================================
# CNF + XOR container
# ============================================================

class CnfXor:
    """
    Container for a CNF + XOR-clause formula.

    Clauses:  lists of nonzero signed ints (DIMACS literals).
    XORs:     (vars, rhs) meaning  XOR of the (positive) vars == rhs.
    """

    def __init__(self) -> None:
        self.nvars: int = 0
        self.clauses: List[List[int]] = []
        self.xors: List[Tuple[List[int], bool]] = []

    def new_var(self) -> int:
        self.nvars += 1
        return self.nvars

    def new_vars(self, count: int) -> List[int]:
        return [self.new_var() for _ in range(count)]

    def add_clause(self, lits: List[int]) -> None:
        self.clauses.append([int(l) for l in lits])

    def add_xor(self, vars_: List[int], rhs: bool) -> None:
        vs = [int(v) for v in vars_]
        if len(vs) == 0:
            if rhs:
                # 0 == 1: trivially unsatisfiable formula.
                self.clauses.append([])
            return
        self.xors.append((vs, bool(rhs)))

    # ---------------- gates ----------------

    def add_and_gate(self, z: int, x: int, y: int) -> None:
        """z <-> (x AND y), Tseitin."""
        self.add_clause([-z, x])
        self.add_clause([-z, y])
        self.add_clause([z, -x, -y])

    def add_at_most_one(self, lits: List[int]) -> None:
        """
        AMO over positive literals: pairwise for short lists,
        sequential (Sinz) encoding otherwise.
        """
        n = len(lits)
        if n <= 1:
            return
        if n <= 6:
            for i in range(n):
                for j in range(i + 1, n):
                    self.add_clause([-lits[i], -lits[j]])
            return
        s = self.new_vars(n - 1)
        self.add_clause([-lits[0], s[0]])
        for i in range(1, n - 1):
            self.add_clause([-lits[i], s[i]])
            self.add_clause([-s[i - 1], s[i]])
            self.add_clause([-lits[i], -s[i - 1]])
        self.add_clause([-lits[n - 1], -s[n - 2]])

    def add_exactly_one(self, lits: List[int]) -> None:
        """
        ALO clause + AMO + native XOR parity = 1.  The parity row is
        redundant given ALO+AMO but places the constraint inside the
        Gauss matrix.
        """
        self.add_clause(list(lits))
        self.add_at_most_one(lits)
        self.add_xor(list(lits), True)

    # ---------------- export ----------------

    def write_dimacs(self, path: str) -> None:
        """
        CryptoMiniSat extended DIMACS: XOR clauses on lines starting
        with 'x'; a false RHS is encoded by negating the first literal.
        """
        with open(path, "w") as f:
            f.write(f"p cnf {self.nvars} {len(self.clauses) + len(self.xors)}\n")
            for c in self.clauses:
                f.write(" ".join(map(str, c)) + " 0\n")
            for vs, rhs in self.xors:
                lits = list(vs)
                if not rhs:
                    lits[0] = -lits[0]
                f.write("x " + " ".join(map(str, lits)) + " 0\n")


# ============================================================
# Encoding
# ============================================================

def encode_instance(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    w_encoding: str = "auto",
    verbose: bool = True,
) -> Tuple[CnfXor, Dict]:
    """
    Build the CNF+XOR formula for one instance (pivot-only certificate).
    Returns (formula, varmaps).
    """
    Gb = gf2(Gb)
    Gp = gf2(Gp)
    p_star = [int(x) for x in p_star]
    k, n = Gb.shape

    if Gp.shape != (n, n):
        raise ValueError(f"Gp must have shape ({n},{n}), got {Gp.shape}.")
    if len(p_star) != k:
        raise ValueError(f"p_star must have length k={k}.")
    if sorted(p_star) != p_star or len(set(p_star)) != k:
        raise ValueError("p_star must be sorted and distinct.")
    if min(p_star) < 0 or max(p_star) >= n:
        raise ValueError("p_star entries out of range.")
    if gf2_rank(Gb) != k:
        raise ValueError("Gb must have full row rank k.")
    if gf2_rank(Gp) != n:
        raise ValueError("Gp must be invertible.")

    w_encoding = str(w_encoding).lower()
    if w_encoding not in {"auto", "butterfly", "dense"}:
        raise ValueError("w_encoding must be auto|butterfly|dense.")
    polar_ok = is_standard_polar_matrix(Gp)
    if w_encoding == "auto":
        actual_w = "butterfly" if polar_ok else "dense"
    elif w_encoding == "butterfly":
        if not polar_ok:
            raise ValueError("butterfly encoding requires standard polar Gp.")
        actual_w = "butterfly"
    else:
        actual_w = "dense"

    cnf = CnfXor()

    # ------------------------------------------------------------
    # 1. Permutation matrix P[a][t], exactly-one per row and column.
    # ------------------------------------------------------------
    P = [[cnf.new_var() for _ in range(n)] for _ in range(n)]
    for a in range(n):
        cnf.add_exactly_one([P[a][t] for t in range(n)])
    for t in range(n):
        cnf.add_exactly_one([P[a][t] for a in range(n)])

    # ------------------------------------------------------------
    # 2. Channeling: b[r,a] = XOR_{t in supp(Gb[r,:])} P[a][t].
    # ------------------------------------------------------------
    supports = [np.flatnonzero(Gb[r]).tolist() for r in range(k)]
    b: Dict[Tuple[int, int], int] = {}
    for r in range(k):
        supp = supports[r]
        for a in range(n):
            v = cnf.new_var()
            b[r, a] = v
            # b XOR (sum over support) = 0.  Empty support => unit b=0.
            cnf.add_xor([v] + [P[a][t] for t in supp], False)

    # ------------------------------------------------------------
    # 3. W on columns 0..p_last via butterfly or dense XORs.
    # ------------------------------------------------------------
    p_last = p_star[-1]
    needed_cols = list(range(p_last + 1))

    w: Dict[Tuple[int, int], int] = {}
    n_butterfly = 0

    if actual_w == "butterfly":
        mlog = polar_log2_size(n)
        x_prev: Dict[Tuple[int, int], int] = {(r, a): b[r, a] for r in range(k) for a in range(n)}
        for t in range(mlog):
            s = 1 << t
            x_next: Dict[Tuple[int, int], int] = {}
            for r in range(k):
                for a in range(n):
                    if (a & s) == 0:
                        out = cnf.new_var()
                        n_butterfly += 1
                        cnf.add_xor([out, x_prev[r, a], x_prev[r, a + s]], False)
                        x_next[r, a] = out
                    else:
                        x_next[r, a] = x_prev[r, a]
            x_prev = x_next
        # Final layer: alias W columns (outside the t-loop).
        for r in range(k):
            for j in needed_cols:
                w[r, j] = x_prev[r, j]
    else:
        for j in needed_cols:
            supp_j = np.flatnonzero(Gp[:, j]).tolist()
            for r in range(k):
                v = cnf.new_var()
                w[r, j] = v
                cnf.add_xor([v] + [b[r, a] for a in supp_j], False)

    # ------------------------------------------------------------
    # 4. Span membership for non-pivot columns j <= p_last.
    # ------------------------------------------------------------
    pivot_set = set(p_star)
    n_lambda = 0
    n_and = 0
    for j in needed_cols:
        if j in pivot_set:
            continue
        prev_pivots = [p for p in p_star if p < j]
        q = len(prev_pivots)
        if q == 0:
            for r in range(k):
                cnf.add_clause([-w[r, j]])
            continue
        lambdas = cnf.new_vars(q)
        n_lambda += q
        for r in range(k):
            ys = []
            for ell, p in enumerate(prev_pivots):
                y = cnf.new_var()
                cnf.add_and_gate(y, lambdas[ell], w[r, p])
                n_and += 1
                ys.append(y)
            cnf.add_xor([w[r, j]] + ys, False)

    # ------------------------------------------------------------
    # 5. Certificate  U W[:,p_star] = I_k.
    # ------------------------------------------------------------
    u = [[cnf.new_var() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for ell, p in enumerate(p_star):
            zs = []
            for r in range(k):
                z = cnf.new_var()
                cnf.add_and_gate(z, u[i][r], w[r, p])
                n_and += 1
                zs.append(z)
            cnf.add_xor(zs, i == ell)

    if verbose:
        print("CNF+XOR formula built.")
        print(f"  n = {n}, k = {k}")
        print(f"  W encoding = {actual_w}")
        print(f"  variables = {cnf.nvars}")
        print(f"  CNF clauses = {len(cnf.clauses)}")
        print(f"  XOR clauses = {len(cnf.xors)}")
        print(f"  P vars = {n*n}, b vars = {k*n}, butterfly vars = {n_butterfly}")
        print(f"  lambda vars = {n_lambda}, AND gates = {n_and}")

    varmaps = {"P": P, "u": u, "w": w, "b": b, "n": n, "k": k}
    return cnf, varmaps


# ============================================================
# Solving
# ============================================================

def random_sampling_prepass(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    max_tries: int = 2000,
    max_seconds: float = 10.0,
    seed: int = 0,
    verbose: bool = True,
) -> Optional[List[int]]:
    """
    Try random permutations before invoking SAT, verifying the profile
    directly over F_2.  Profiles realized by many permutations (the
    common case away from reachability boundaries) are found here in
    milliseconds; proof search is reserved for the rare-solution and
    UNSAT regimes, which is where it earns its cost.
    """
    Gb = gf2(Gb)
    Gp = gf2(Gp)
    p_star = [int(x) for x in p_star]
    k, n = Gb.shape
    rng = np.random.default_rng(seed)
    deadline = time.time() + float(max_seconds)
    t = 0
    for t in range(int(max_tries)):
        if time.time() > deadline:
            break
        pi = rng.permutation(n).tolist()
        W = gf2_matmul(Gb[:, pi], Gp)
        if column_pivot_profile_gf2(W) == p_star:
            if verbose:
                print(f"  sampling pre-pass: hit on try {t+1}")
            return pi
    if verbose:
        print(f"  sampling pre-pass: no hit in {t+1} tries; "
              f"falling through to SAT.")
    return None


def decode_pi(solution, P: List[List[int]], n: int) -> List[int]:
    pi: List[int] = []
    for a in range(n):
        hits = [t for t in range(n) if solution[P[a][t]]]
        if len(hits) != 1:
            raise RuntimeError(f"Row {a} of P has {len(hits)} true entries; encoding bug.")
        pi.append(hits[0])
    if sorted(pi) != list(range(n)):
        raise RuntimeError("Decoded pi is not a permutation; encoding bug.")
    return pi


def solve_instance_cms(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    time_limit_sec: float = 3600.0,
    threads: int = 8,
    w_encoding: str = "auto",
    sat_verbosity: int = 0,
    enumerate_count: int = 1,
    dimacs_path: Optional[str] = None,
    prepass: bool = True,
    prepass_tries: int = 2000,
    verbose: bool = True,
) -> Tuple[str, List[Tuple[List[int], np.ndarray]]]:
    """
    Solve one instance.  Returns (status, solutions) where status is
    'SAT' | 'UNSAT' | 'UNKNOWN' and solutions is a list of (pi, U).
    With enumerate_count > 1, previously found permutations are blocked
    incrementally and up to that many distinct pi are returned
    (status 'UNSAT' after >=1 solution means the list is EXHAUSTIVE).
    """
    if prepass and enumerate_count <= 1:
        pi_hit = random_sampling_prepass(
            Gb, Gp, p_star, max_tries=prepass_tries, verbose=verbose)
        if pi_hit is not None:
            k0 = gf2(Gb).shape[0]
            W0 = gf2_matmul(gf2(Gb)[:, pi_hit], gf2(Gp))
            # U = W[:,p_star]^{-1} exists; return identity placeholder is
            # wrong, so compute it by elimination for completeness.
            Wp = W0[:, [int(x) for x in p_star]]
            aug = np.concatenate([Wp, np.eye(k0, dtype=np.uint8)], axis=1)
            R, piv = gf2_rref(aug)
            U0 = R[:, k0:]
            return "SAT", [(pi_hit, U0)]

    cnf, vm = encode_instance(Gb, Gp, p_star, w_encoding=w_encoding, verbose=verbose)
    n, k = vm["n"], vm["k"]
    P, u = vm["P"], vm["u"]

    if dimacs_path is not None:
        cnf.write_dimacs(dimacs_path)
        print(f"Extended DIMACS written to {dimacs_path} "
              f"({cnf.nvars} vars, {len(cnf.clauses)} clauses, {len(cnf.xors)} xors).")

    if not HAVE_PYCRYPTOSAT:
        raise RuntimeError(
            "pycryptosat is not installed (pip install pycryptosat). "
            "Alternatively use --dimacs and run the cryptominisat5 binary."
        )

    solver = CMSSolver(threads=int(threads),
                       time_limit=float(time_limit_sec),
                       verbose=int(sat_verbosity))
    for c in cnf.clauses:
        if len(c) == 0:
            return "UNSAT", []
        solver.add_clause(c)
    for vs, rhs in cnf.xors:
        solver.add_xor_clause(vs, rhs)

    solutions: List[Tuple[List[int], np.ndarray]] = []
    deadline = time.time() + float(time_limit_sec)

    while len(solutions) < max(1, int(enumerate_count)):
        t0 = time.time()
        sat, model = solver.solve()
        dt = time.time() - t0

        if sat is True:
            pi_sol = decode_pi(model, P, n)
            U_sol = np.array([[1 if model[u[i][r]] else 0 for r in range(k)]
                              for i in range(k)], dtype=np.uint8)
            solutions.append((pi_sol, U_sol))
            if verbose:
                print(f"  solution {len(solutions)} found in {dt:.2f}s")
            if len(solutions) >= max(1, int(enumerate_count)):
                return "SAT", solutions
            # Block this permutation: some position must map differently.
            solver.add_clause([-P[a][pi_sol[a]] for a in range(n)])
            if time.time() >= deadline:
                return "SAT", solutions
        elif sat is False:
            if solutions:
                print("  enumeration exhausted: no further permutations "
                      "realize this profile.")
                return "SAT", solutions
            return "UNSAT", []
        else:
            return ("SAT" if solutions else "UNKNOWN"), solutions

    return ("SAT" if solutions else "UNKNOWN"), solutions


def solve_reconstruction_cms(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    time_limit_sec: float = 3600.0,
    threads: int = 8,
    w_encoding: str = "auto",
    dualize: str = "auto",
    sat_verbosity: int = 0,
    enumerate_count: int = 1,
    dimacs_path: Optional[str] = None,
    prepass: bool = True,
    prepass_tries: int = 2000,
    verbose: bool = True,
) -> Tuple[str, List[List[int]]]:
    """
    Top-level solve with optional matroid-duality reduction.
    Returns (status, list_of_pi).  UNSAT means p_star is PROVABLY
    unreachable (the duality reduction is an equivalence, so UNSAT on
    the dual is UNSAT on the primal).
    """
    Gb = gf2(Gb)
    Gp = gf2(Gp)
    p_star = [int(x) for x in p_star]
    k, n = Gb.shape

    dualize = str(dualize).lower()
    if dualize not in {"auto", "always", "never"}:
        raise ValueError("dualize must be auto|always|never.")
    polar_ok = is_standard_polar_matrix(Gp)
    if dualize == "always" and not polar_ok:
        raise ValueError("--dualize always requires standard polar Gp.")

    use_dual = False
    if polar_ok and 0 < n - k < n:
        if dualize == "always":
            use_dual = True
        elif dualize == "auto":
            use_dual = (n - k) < k

    if not use_dual:
        if verbose:
            print(f"Solving PRIMAL instance: n={n}, k={k}.")
        status, sols = solve_instance_cms(
            Gb, Gp, p_star,
            time_limit_sec=time_limit_sec, threads=threads,
            w_encoding=w_encoding, sat_verbosity=sat_verbosity,
            enumerate_count=enumerate_count, dimacs_path=dimacs_path,
            prepass=prepass, prepass_tries=prepass_tries,
            verbose=verbose,
        )
        return status, [pi for (pi, _) in sols]

    Hb, q_star = dualize_instance(Gb, p_star)
    if verbose:
        print(f"Matroid-duality reduction: primal (n={n}, k={k}) "
              f"-> dual (n={n}, k'={n-k}).")
        print(f"  dual target q_star = {q_star}")

    status, sols = solve_instance_cms(
        Hb, Gp, q_star,
        time_limit_sec=time_limit_sec, threads=threads,
        w_encoding=w_encoding, sat_verbosity=sat_verbosity,
        enumerate_count=enumerate_count, dimacs_path=dimacs_path,
        prepass=prepass, prepass_tries=prepass_tries,
        verbose=verbose,
    )

    pis: List[List[int]] = []
    for sigma_sol, U_dual in sols:
        # Verify the dual solution against the dual instance.
        Wd = gf2_matmul(Hb[:, sigma_sol], Gp)
        if column_pivot_profile_gf2(Wd) != q_star:
            print("WARNING: dual solution fails dual profile check; skipping.")
            continue
        pis.append(primal_pi_from_dual_sigma(sigma_sol))

    if status == "SAT" and not pis:
        return "UNKNOWN", []
    return status, pis


# ============================================================
# Verification and self-tests
# ============================================================

def verify_pi(Gb: np.ndarray, Gp: np.ndarray, p_star: List[int],
              pi_sol: List[int], label: str = "") -> bool:
    W = gf2_matmul(gf2(Gb)[:, pi_sol], gf2(Gp))
    actual = column_pivot_profile_gf2(W)
    ok = (actual == [int(x) for x in p_star])
    tag = f" [{label}]" if label else ""
    print(f"Target pivot profile{tag}: {list(p_star)}")
    print(f"Actual pivot profile{tag}: {actual}")
    print(f"Pivot profile check{tag}:", "PASS" if ok else "FAIL")
    return ok


def self_test_duality(trials: int = 40, seed: int = 1) -> None:
    rng = np.random.default_rng(seed)
    checked = 0
    for m in (3, 4, 5):
        n = 2 ** m
        F = polar_matrix(m)
        for _ in range(trials):
            k = int(rng.integers(1, n))
            Gb = random_full_row_rank_binary_matrix(k, n, rng)
            pi = list(rng.permutation(n))
            p = column_pivot_profile_gf2(gf2_matmul(Gb[:, pi], F))
            Hb, q = dualize_instance(Gb, p)
            sigma = list(reversed(pi))
            q_check = column_pivot_profile_gf2(gf2_matmul(Hb[:, sigma], F))
            if q_check != q:
                raise AssertionError(f"duality check failed at n={n}, k={k}")
            checked += 1
    print(f"[self-test] duality identity verified on {checked} random instances.  PASS")


def built_in_self_test() -> Tuple[np.ndarray, np.ndarray, List[int], List[int]]:
    rng = np.random.default_rng(0)
    m = 4
    n = 2 ** m
    k = 10
    Gp = polar_matrix(m)
    Gb = random_full_row_rank_binary_matrix(k, n, rng)
    true_pi = list(rng.permutation(n))
    p_star = column_pivot_profile_gf2(gf2_matmul(Gb[:, true_pi], Gp))
    print("Built-in self-test instance:")
    print(f"  n = {n}, k = {k}  (k > n/2, dual path active under auto)")
    print(f"  p_star = {p_star}")
    return Gb, Gp, p_star, true_pi


def load_npz_instance(path: str) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    data = np.load(path, allow_pickle=False)
    if "Gb" not in data or "Gp" not in data or "p_star" not in data:
        raise ValueError("npz file must contain arrays Gb, Gp, and p_star.")
    Gb = gf2(data["Gb"])
    Gp = gf2(data["Gp"])
    p_star = [int(x) for x in np.asarray(data["p_star"]).reshape(-1)]
    return Gb, Gp, p_star


# ============================================================
# CLI
# ============================================================

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("npz", nargs="?", default=None,
                        help="Optional .npz with Gb, Gp, p_star.")
    parser.add_argument("--time", type=float, default=3600.0,
                        help="Time limit in seconds.")
    parser.add_argument("--threads", "--workers", type=int, default=8,
                        help="CryptoMiniSat threads.")
    parser.add_argument("--w-encoding", choices=["auto", "butterfly", "dense"],
                        default="auto")
    parser.add_argument("--dualize", choices=["auto", "always", "never"],
                        default="auto")
    parser.add_argument("--enumerate", type=int, default=1, dest="enumerate_count",
                        help="Find up to N distinct permutations (default 1).")
    parser.add_argument("--dimacs", type=str, default=None,
                        help="Also export the formula in extended DIMACS "
                             "for the cryptominisat5 binary.")
    parser.add_argument("--sat-verbose", type=int, default=0,
                        help="CryptoMiniSat verbosity (0..15).")
    parser.add_argument("--no-prepass", action="store_true",
                        help="Skip the random-sampling pre-pass.")
    parser.add_argument("--prepass-tries", type=int, default=2000)
    parser.add_argument("--skip-duality-selftest", action="store_true")
    args = parser.parse_args()

    if not args.skip_duality_selftest:
        self_test_duality()
        print()

    if args.npz is None:
        Gb, Gp, p_star, true_pi = built_in_self_test()
    else:
        Gb, Gp, p_star = load_npz_instance(args.npz)
        true_pi = None

    print()
    print("Solving reconstruction problem with CryptoMiniSat...")
    t0 = time.time()
    status, pis = solve_reconstruction_cms(
        Gb, Gp, p_star,
        time_limit_sec=args.time,
        threads=args.threads,
        w_encoding=args.w_encoding,
        dualize=args.dualize,
        sat_verbosity=args.sat_verbose,
        enumerate_count=args.enumerate_count,
        dimacs_path=args.dimacs,
        prepass=not args.no_prepass,
        prepass_tries=args.prepass_tries,
        verbose=True,
    )
    dt = time.time() - t0
    print(f"Status: {status}  ({dt:.2f}s)")

    if status == "UNSAT":
        print("p_star is PROVABLY UNREACHABLE for this (Gb, Gp).")
        return 1
    if status == "UNKNOWN":
        print("No answer within the time limit.")
        return 1

    all_ok = True
    for idx, pi_sol in enumerate(pis):
        print()
        print(f"Recovered pi #{idx+1}:")
        print(pi_sol)
        ok = verify_pi(Gb, Gp, p_star, pi_sol, label=f"pi #{idx+1}")
        all_ok = all_ok and ok

    if true_pi is not None:
        print()
        print("Original true pi used to generate p_star:")
        print([int(x) for x in true_pi])
        print("Note: recovered pi need not equal true pi.")

    return 0 if all_ok else 2


if __name__ == "__main__":
    sys.exit(main())