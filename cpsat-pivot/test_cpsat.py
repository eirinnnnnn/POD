#!/usr/bin/env python3
"""
CP-SAT reconstruction of a permutation pi such that

    PivotProfile( pi(G_b) G_p ) = p_star.

Convention:
    pi(G_b)[r, a] = G_b[r, pi[a]]

So if pi is represented as a list

    pi = [pi[0], ..., pi[n-1]],

then

    pi(G_b) = G_b[:, pi].

All arithmetic for G_b, G_p, W, U, M is over F_2.

Install dependencies:

    pip install ortools numpy

Usage 1: run built-in small self-test

    python cpsat_pivot_reconstruct.py

Usage 2: provide your own matrices in an npz file

    python cpsat_pivot_reconstruct.py data.npz --time 300

where data.npz contains arrays:

    Gb      shape (k,n), entries 0/1
    Gp      shape (n,n), entries 0/1
    p_star  shape (k,), sorted 0-based target pivot positions
"""

import argparse
import sys
from typing import Dict, List, Tuple, Optional

import numpy as np
from ortools.sat.python import cp_model


# ============================================================
# Basic F_2 utilities
# ============================================================

def gf2(A: np.ndarray) -> np.ndarray:
    return (np.asarray(A, dtype=np.uint8) & 1)


def gf2_matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Matrix product over F_2.
    uint8 overflow is harmless for parity because 256 is even.
    """
    A = gf2(A)
    B = gf2(B)
    return ((A @ B) & 1).astype(np.uint8)


def gf2_rank(A: np.ndarray) -> int:
    """
    Rank over F_2 by row elimination.
    """
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


def column_pivot_profile_gf2(A: np.ndarray) -> List[int]:
    """
    Return the column rank profile of A over F_2.

    Scans columns from left to right and records j whenever
    rank(A[:,0:j]) jumps by one.
    """
    A = gf2(A)
    k, n = A.shape

    # basis[p] stores a vector whose highest set bit is p.
    basis = [0 for _ in range(k)]
    pivots: List[int] = []

    for j in range(n):
        v = 0
        for r in range(k):
            if A[r, j]:
                v |= (1 << r)

        # Reduce by current basis.
        for p in reversed(range(k)):
            if ((v >> p) & 1) and basis[p] != 0:
                v ^= basis[p]

        if v != 0:
            p = v.bit_length() - 1
            basis[p] = v

            # Keep basis reduced.
            for q in range(k):
                if q != p and basis[q] != 0 and ((basis[q] >> p) & 1):
                    basis[q] ^= v

            pivots.append(j)

    return pivots


def polar_matrix(m: int) -> np.ndarray:
    """
    Standard Kronecker polar matrix F^{\otimes m} over F_2,
    with F = [[1,0],[1,1]].

    This is only for the built-in self-test.
    Replace by your own G_p in real experiments.
    """
    F = np.array([[1, 0], [1, 1]], dtype=np.uint8)
    G = np.array([[1]], dtype=np.uint8)
    for _ in range(m):
        G = np.kron(G, F).astype(np.uint8) & 1
    return G


def random_full_row_rank_binary_matrix(k: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """
    Generate random G_b in F_2^{k x n} with full row rank k.
    """
    while True:
        Gb = rng.integers(0, 2, size=(k, n), dtype=np.uint8)
        if gf2_rank(Gb) == k:
            return Gb


# ============================================================
# CP-SAT helpers
# ============================================================

def add_xor_eq(
    model: cp_model.CpModel,
    lits: List[cp_model.IntVar],
    rhs: int,
    true_lit: cp_model.IntVar,
) -> None:
    """
    Enforce XOR(lits) == rhs over F_2.

    OR-Tools AddBoolXOr(lits) enforces XOR(lits) == 1.

    Therefore:
        XOR(lits) == 1  -> AddBoolXOr(lits)
        XOR(lits) == 0  -> AddBoolXOr(lits + [true_lit])
    because XOR(lits) XOR 1 == 1 iff XOR(lits) == 0.
    """
    rhs = int(rhs) & 1
    lits = list(lits)

    if len(lits) == 0:
        # Empty XOR is 0.
        if rhs == 0:
            return
        else:
            model.Add(true_lit == 0)  # contradiction
            return

    if rhs == 1:
        model.AddBoolXOr(lits)
    else:
        model.AddBoolXOr(lits + [true_lit])


def add_and_eq(
    model: cp_model.CpModel,
    z: cp_model.IntVar,
    x: cp_model.IntVar,
    y: cp_model.IntVar,
) -> None:
    """
    Enforce z = x AND y for Boolean variables.
    """
    model.Add(z <= x)
    model.Add(z <= y)
    model.Add(z >= x + y - 1)


# ============================================================
# CP-SAT model
# ============================================================

def build_pivot_reconstruction_model(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    verbose: bool = True,
):
    """
    Build the CP-SAT feasibility model.

    Variables:

        pi_a in {0,...,n-1}
            permutation variable, pi_a = pi(a)

        u_{i,r} in {0,1}
            row-operation certificate U

        b_{r,a} in {0,1}
            b_{r,a} = G_b[r, pi_a]

        w_{r,j} in {0,1}
            w_{r,j} = XOR_{a: G_p[a,j]=1} b_{r,a}

        z_{i,r,j} in {0,1}
            z_{i,r,j} = u_{i,r} AND w_{r,j}

    Constraints:

        AllDifferent(pi_0,...,pi_{n-1})

        m_{i,j} := XOR_r z_{i,r,j}

        m_{i,p_l} = delta_{i,l}

        m_{i,j} = 0 for j < p_i
    """
    Gb = gf2(Gb)
    Gp = gf2(Gp)
    p_star = [int(x) for x in p_star]

    k, n = Gb.shape

    if Gp.shape != (n, n):
        raise ValueError(f"Gp must have shape ({n},{n}), got {Gp.shape}.")

    if len(p_star) != k:
        raise ValueError(f"p_star must have length k={k}, got {len(p_star)}.")

    if sorted(p_star) != p_star:
        raise ValueError("p_star must be sorted increasingly.")

    if len(set(p_star)) != len(p_star):
        raise ValueError("p_star must contain distinct entries.")

    if min(p_star) < 0 or max(p_star) >= n:
        raise ValueError("p_star entries must lie in {0,...,n-1}.")

    if gf2_rank(Gb) != k:
        raise ValueError("Gb must have full row rank k over F_2.")

    if gf2_rank(Gp) != n:
        raise ValueError("Gp must be invertible over F_2.")

    model = cp_model.CpModel()

    true_lit = model.NewBoolVar("TRUE")
    model.Add(true_lit == 1)

    # --------------------------------------------------------
    # 1. Permutation variables
    # --------------------------------------------------------

    pi = [
        model.NewIntVar(0, n - 1, f"pi[{a}]")
        for a in range(n)
    ]

    model.AddAllDifferent(pi)

    # --------------------------------------------------------
    # 2. U variables
    # --------------------------------------------------------

    u: Dict[Tuple[int, int], cp_model.IntVar] = {}

    for i in range(k):
        for r in range(k):
            u[i, r] = model.NewBoolVar(f"u[{i},{r}]")

    # --------------------------------------------------------
    # 3. Determine which M columns are constrained
    # --------------------------------------------------------

    # Constraint dictionary:
    #     fixed_m[(i,j)] = rhs
    # meaning:
    #     m_{i,j} = rhs
    fixed_m: Dict[Tuple[int, int], int] = {}

    def fix_m(i: int, j: int, rhs: int) -> None:
        rhs = int(rhs) & 1
        key = (i, j)
        if key in fixed_m and fixed_m[key] != rhs:
            raise ValueError(f"Contradictory constraints on m[{i},{j}].")
        fixed_m[key] = rhs

    # Identity constraints:
    #     M[:, p_star] = I_k
    for i in range(k):
        for ell, p in enumerate(p_star):
            fix_m(i, p, 1 if i == ell else 0)

    # Pre-pivot zero constraints:
    #     M[i,j] = 0 for j < p_i
    for i in range(k):
        for j in range(p_star[i]):
            fix_m(i, j, 0)

    constrained_cols = sorted({j for (_, j) in fixed_m.keys()})

    # --------------------------------------------------------
    # 4. b_{r,a} = G_b[r, pi_a]
    # --------------------------------------------------------

    b: Dict[Tuple[int, int], cp_model.IntVar] = {}

    for r in range(k):
        row_values = [int(Gb[r, t]) for t in range(n)]

        row_all_zero = all(v == 0 for v in row_values)
        row_all_one = all(v == 1 for v in row_values)

        for a in range(n):
            b[r, a] = model.NewBoolVar(f"b[{r},{a}]")

            if row_all_zero:
                model.Add(b[r, a] == 0)
            elif row_all_one:
                model.Add(b[r, a] == 1)
            else:
                # b[r,a] = row_values[pi[a]]
                model.AddElement(pi[a], row_values, b[r, a])

    # --------------------------------------------------------
    # 5. w_{r,j} = XOR_{a: G_p[a,j]=1} b_{r,a}
    # --------------------------------------------------------

    w: Dict[Tuple[int, int], cp_model.IntVar] = {}

    for r in range(k):
        for j in constrained_cols:
            w[r, j] = model.NewBoolVar(f"w[{r},{j}]")

            support = [a for a in range(n) if int(Gp[a, j]) == 1]
            terms = [b[r, a] for a in support]

            # w = XOR(terms)
            # equivalently XOR(terms + [w]) = 0
            add_xor_eq(model, terms + [w[r, j]], 0, true_lit)

    # --------------------------------------------------------
    # 6. z_{i,r,j} = u_{i,r} AND w_{r,j}
    # --------------------------------------------------------

    z: Dict[Tuple[int, int, int], cp_model.IntVar] = {}

    for i in range(k):
        for r in range(k):
            for j in constrained_cols:
                z[i, r, j] = model.NewBoolVar(f"z[{i},{r},{j}]")
                add_and_eq(model, z[i, r, j], u[i, r], w[r, j])

    # --------------------------------------------------------
    # 7. m_{i,j} = XOR_r z_{i,r,j}, and fix constrained m
    # --------------------------------------------------------

    for (i, j), rhs in fixed_m.items():
        terms = [z[i, r, j] for r in range(k)]
        add_xor_eq(model, terms, rhs, true_lit)

    # Optional decision strategy: branch on the permutation first.
    # This often makes sense because U is only a certificate.
    model.AddDecisionStrategy(
        pi,
        cp_model.CHOOSE_FIRST,
        cp_model.SELECT_MIN_VALUE,
    )

    if verbose:
        print("CP-SAT model built.")
        print(f"  n = {n}")
        print(f"  k = {k}")
        print(f"  |p_star| = {len(p_star)}")
        print(f"  constrained M entries = {len(fixed_m)}")
        print(f"  constrained columns = {len(constrained_cols)}")
        print(f"  permutation vars = {len(pi)} integer vars")
        print(f"  U vars = {k*k}")
        print(f"  b vars = {k*n}")
        print(f"  w vars = {k*len(constrained_cols)}")
        print(f"  z vars = {k*k*len(constrained_cols)}")

    variables = {
        "pi": pi,
        "u": u,
        "b": b,
        "w": w,
        "z": z,
        "fixed_m": fixed_m,
        "constrained_cols": constrained_cols,
    }

    return model, variables


# ============================================================
# Solver and verification
# ============================================================

def solve_reconstruction(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    *,
    time_limit_sec: float = 60.0,
    workers: int = 8,
    log: bool = False,
    verbose: bool = True,
) -> Optional[Tuple[List[int], np.ndarray]]:
    """
    Solve the reconstruction problem.

    Returns:
        (pi_solution, U_solution) if feasible,
        None otherwise.
    """
    model, variables = build_pivot_reconstruction_model(
        Gb,
        Gp,
        p_star,
        verbose=verbose,
    )

    pi_vars = variables["pi"]
    u_vars = variables["u"]

    k, n = Gb.shape

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_sec)
    solver.parameters.num_search_workers = int(workers)
    solver.parameters.log_search_progress = bool(log)

    status = solver.Solve(model)

    status_name = solver.StatusName(status)
    print(f"Solver status: {status_name}")

    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        print("No feasible solution found within the given limits.")
        return None

    pi_sol = [solver.Value(pi_vars[a]) for a in range(n)]

    U_sol = np.zeros((k, k), dtype=np.uint8)
    for i in range(k):
        for r in range(k):
            U_sol[i, r] = solver.Value(u_vars[i, r]) & 1

    return pi_sol, U_sol


def verify_solution(
    Gb: np.ndarray,
    Gp: np.ndarray,
    p_star: List[int],
    pi_sol: List[int],
    U_sol: Optional[np.ndarray] = None,
    *,
    print_matrices: bool = False,
) -> bool:
    """
    Verify that pi realizes the target pivot profile.
    Also verifies the U-certificate if U_sol is given.
    """
    Gb = gf2(Gb)
    Gp = gf2(Gp)
    p_star = [int(x) for x in p_star]

    k, n = Gb.shape

    piGb = Gb[:, pi_sol]
    W = gf2_matmul(piGb, Gp)

    actual_pivots = column_pivot_profile_gf2(W)

    print("Target pivot profile: ", p_star)
    print("Actual pivot profile: ", actual_pivots)

    ok_profile = (actual_pivots == p_star)
    print("Pivot profile check:", "PASS" if ok_profile else "FAIL")

    ok_cert = True

    if U_sol is not None:
        U_sol = gf2(U_sol)
        M = gf2_matmul(U_sol, W)

        ok_identity = True
        for i in range(k):
            for ell, p in enumerate(p_star):
                expected = 1 if i == ell else 0
                if int(M[i, p]) != expected:
                    ok_identity = False

        ok_zeros = True
        for i in range(k):
            for j in range(p_star[i]):
                if int(M[i, j]) != 0:
                    ok_zeros = False

        ok_cert = ok_identity and ok_zeros

        print("U-certificate identity check:", "PASS" if ok_identity else "FAIL")
        print("U-certificate zero check:    ", "PASS" if ok_zeros else "FAIL")

        if print_matrices:
            print("pi(Gb) =")
            print(piGb)
            print("W = pi(Gb) Gp =")
            print(W)
            print("U =")
            print(U_sol)
            print("M = U W =")
            print(M)

    return ok_profile and ok_cert


# ============================================================
# Demo / CLI
# ============================================================

def built_in_self_test() -> Tuple[np.ndarray, np.ndarray, List[int], List[int]]:
    """
    Make a small reachable test instance.

    We choose random full-rank G_b, polar G_p, random true pi,
    then define p_star from that true pi. Therefore p_star is reachable.
    """
    rng = np.random.default_rng(0)

    m = 3
    n = 2**m
    k = 3

    Gp = polar_matrix(m)
    Gb = random_full_row_rank_binary_matrix(k, n, rng)

    true_pi = list(rng.permutation(n))

    W_true = gf2_matmul(Gb[:, true_pi], Gp)
    p_star = column_pivot_profile_gf2(W_true)

    if len(p_star) != k:
        raise RuntimeError("Unexpected pivot profile length in self-test.")

    print("Built-in self-test instance:")
    print(f"  n = {n}, k = {k}")
    print(f"  true pi = {true_pi}")
    print(f"  p_star from true pi = {p_star}")

    return Gb, Gp, p_star, true_pi


def load_npz_instance(path: str) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    data = np.load(path, allow_pickle=False)

    if "Gb" not in data or "Gp" not in data or "p_star" not in data:
        raise ValueError("npz file must contain arrays Gb, Gp, and p_star.")

    Gb = gf2(data["Gb"])
    Gp = gf2(data["Gp"])
    p_star = [int(x) for x in np.asarray(data["p_star"]).reshape(-1)]

    return Gb, Gp, p_star


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "npz",
        nargs="?",
        default=None,
        help="Optional .npz file containing Gb, Gp, p_star.",
    )

    parser.add_argument(
        "--time",
        type=float,
        default=60.0,
        help="CP-SAT time limit in seconds.",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of CP-SAT search workers.",
    )

    parser.add_argument(
        "--log",
        action="store_true",
        help="Enable CP-SAT search log.",
    )

    parser.add_argument(
        "--print-matrices",
        action="store_true",
        help="Print pi(Gb), W, U, M after solving.",
    )

    args = parser.parse_args()

    if args.npz is None:
        Gb, Gp, p_star, true_pi = built_in_self_test()
    else:
        Gb, Gp, p_star = load_npz_instance(args.npz)
        true_pi = None

    print()
    print("Solving reconstruction problem...")
    result = solve_reconstruction(
        Gb,
        Gp,
        p_star,
        time_limit_sec=args.time,
        workers=args.workers,
        log=args.log,
        verbose=True,
    )

    if result is None:
        return 1

    pi_sol, U_sol = result

    print()
    print("Recovered pi:")
    print(pi_sol)

    if true_pi is not None:
        print()
        print("Original true pi used to generate p_star:")
        print(true_pi)
        print()
        print("Note: recovered pi need not equal true pi; many pi can share the same pivot profile.")

    print()
    ok = verify_solution(
        Gb,
        Gp,
        p_star,
        pi_sol,
        U_sol,
        print_matrices=args.print_matrices,
    )

    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())