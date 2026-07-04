#!/usr/bin/env python3
"""
Backtracking reconstruction of a permutation pi such that

    PivotSet(Gb @ P_pi @ Gp) = S_star

over F_2, using the lower-triangular/suffix-sum structure of Gp.

Convention:
    P_pi[i, j] = 1 iff pi[j] = i.
Therefore:
    (Gb @ P_pi)[:, j] = Gb[:, pi[j]].

The solver assigns the permuted columns

    b_j = (Gb @ P_pi)[:, j]

from right to left. Since Gp is lower triangular with unit diagonal,

    a_j = (Gb @ P_pi @ Gp)[:, j]
        = b_j + sum_{ell > j} Gp[ell, j] b_ell.

The target prefix-pivot profile is converted to layer constraints and enforced
with sound subspace pruning.
"""

from __future__ import annotations

import argparse
import collections
import dataclasses
import json
import random
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# GF(2) vector/subspace utilities.
# Vectors in F_2^k are represented as Python integers.
# Bit i is the i-th coordinate.
# ---------------------------------------------------------------------------

class XorSpace:
    """Linear subspace of F_2^k stored as a row-reduction XOR basis."""

    __slots__ = ("k", "basis", "dim")

    def __init__(self, k: int, basis: Optional[List[int]] = None, dim: Optional[int] = None):
        self.k = k
        if basis is None:
            self.basis = [0] * k
            self.dim = 0
        else:
            self.basis = list(basis)
            self.dim = sum(1 for x in self.basis if x) if dim is None else int(dim)

    def copy(self) -> "XorSpace":
        return XorSpace(self.k, self.basis, self.dim)

    def contains(self, x: int) -> bool:
        """Return True iff x is in the span."""
        y = int(x)
        for p in range(self.k - 1, -1, -1):
            if (y >> p) & 1:
                bp = self.basis[p]
                if bp:
                    y ^= bp
                else:
                    return False
        return True

    def add(self, x: int) -> bool:
        """
        Add x to the basis.

        Returns True iff the dimension increased.
        """
        y = int(x)
        for p in range(self.k - 1, -1, -1):
            if (y >> p) & 1:
                bp = self.basis[p]
                if bp:
                    y ^= bp
                else:
                    self.basis[p] = y
                    self.dim += 1
                    return True
        return False

    def snapshot(self) -> Tuple[Tuple[int, ...], int]:
        return tuple(self.basis), self.dim


def col_to_int(M: np.ndarray, j: int) -> int:
    """Convert column j of a binary matrix to an integer vector."""
    k = M.shape[0]
    x = 0
    for i in range(k):
        if int(M[i, j]) & 1:
            x |= 1 << i
    return x


def matrix_cols_to_ints(M: np.ndarray) -> List[int]:
    return [col_to_int(M, j) for j in range(M.shape[1])]


def int_to_col(x: int, k: int) -> np.ndarray:
    return np.array([(x >> i) & 1 for i in range(k)], dtype=np.uint8)


def pivot_set_from_col_ints(cols: Sequence[int], k: int) -> List[int]:
    """Column pivot profile by left-to-right greedy rank increments."""
    space = XorSpace(k)
    pivots: List[int] = []
    for j, x in enumerate(cols):
        if not space.contains(x):
            space.add(x)
            pivots.append(j)
    return pivots


def permutation_matrix_from_pi(pi: Sequence[int], n: int) -> np.ndarray:
    """
    Build P with P[pi[j], j] = 1.
    Then (Gb @ P)[:, j] = Gb[:, pi[j]].
    """
    P = np.zeros((n, n), dtype=np.uint8)
    for j, i in enumerate(pi):
        P[int(i), j] = 1
    return P


def pi_from_permutation_matrix(P: np.ndarray) -> List[int]:
    """Return pi satisfying P[pi[j], j] = 1."""
    P = np.asarray(P)
    pi: List[int] = []
    for j in range(P.shape[1]):
        rows = np.flatnonzero(P[:, j] & 1)
        if len(rows) != 1:
            raise ValueError(f"Column {j} of permutation matrix has {len(rows)} ones.")
        pi.append(int(rows[0]))
    return pi


# ---------------------------------------------------------------------------
# Solver result.
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class SearchResult:
    status: str
    pi: Optional[List[int]]
    pivot_set: Optional[List[int]]
    nodes: int
    elapsed_sec: float
    best_remaining_j: int
    message: str


class PivotBacktracker:
    """
    Exact backtracking solver with sound pruning.

    It may still be exponential. Use time/node limits and optionally a hint
    permutation for candidate ordering.
    """

    def __init__(
        self,
        Gb: np.ndarray,
        Gp: np.ndarray,
        target: Sequence[int],
        *,
        seed: int = 0,
        verbose: bool = False,
    ):
        self.Gb = (np.asarray(Gb, dtype=np.uint8) & 1)
        self.Gp = (np.asarray(Gp, dtype=np.uint8) & 1)
        self.k, self.n = self.Gb.shape

        if self.Gp.shape != (self.n, self.n):
            raise ValueError(f"Gp shape must be {(self.n, self.n)}, got {self.Gp.shape}.")

        if not np.all(np.diag(self.Gp) == 1):
            raise ValueError("This backtracker assumes diag(Gp)=1.")

        # Lower-triangular in the convention used here: Gp[ell, j] = 0 for ell < j.
        if not np.all(np.triu(self.Gp, 1) == 0):
            raise ValueError("This backtracker assumes Gp is lower triangular.")

        self.target = sorted(int(x) for x in target)
        if len(self.target) != self.k:
            raise ValueError(
                f"Target pivot set should have length k={self.k}; got {len(self.target)}."
            )

        self.target_set = set(self.target)
        self.rank_target = [
            sum(1 for s in self.target if s <= j) for j in range(self.n)
        ]
        # 1-based pivot layer: s_i has layer i.
        self.pivot_rank: Dict[int, int] = {
            s: i + 1 for i, s in enumerate(self.target)
        }

        # support_gt[j] = {ell > j : Gp[ell,j] = 1}
        self.support_gt: List[List[int]] = [
            [int(x) for x in np.flatnonzero(self.Gp[j + 1 :, j]) + (j + 1)]
            for j in range(self.n)
        ]

        self.col_values = matrix_cols_to_ints(self.Gb)

        self.indices_by_value: Dict[int, List[int]] = collections.defaultdict(list)
        for idx, v in enumerate(self.col_values):
            self.indices_by_value[v].append(idx)

        self.rng = random.Random(seed)
        self.verbose = verbose

    def _recover_pi_from_b(self, b: Sequence[int]) -> List[int]:
        available = {v: list(idxs) for v, idxs in self.indices_by_value.items()}
        pi: List[int] = []
        for v in b:
            if v not in available or not available[v]:
                raise RuntimeError("Internal error: selected column multiset is invalid.")
            pi.append(available[v].pop())
        return pi

    def _compute_A_cols_from_b(self, b: Sequence[int]) -> List[int]:
        A_cols: List[int] = []
        for j in range(self.n):
            x = 0
            # include ell >= j
            for ell in range(j, self.n):
                if int(self.Gp[ell, j]) & 1:
                    x ^= int(b[ell])
            A_cols.append(x)
        return A_cols

    def verify_pi(self, pi: Sequence[int]) -> Tuple[bool, List[int]]:
        b = [self.col_values[int(pi[j])] for j in range(self.n)]
        A_cols = self._compute_A_cols_from_b(b)
        pivots = pivot_set_from_col_ints(A_cols, self.k)
        return pivots == self.target, pivots

    def solve(
        self,
        *,
        time_limit_sec: float = 60.0,
        node_limit: int = 1_000_000,
        hint_pi: Optional[Sequence[int]] = None,
    ) -> SearchResult:
        """
        Run DFS.

        hint_pi, if provided, is used only to order candidates; every step is still
        checked by the mathematical constraints.
        """
        start = time.time()
        sys.setrecursionlimit(max(10000, 10 * self.n))

        K: List[XorSpace] = [XorSpace(self.k) for _ in range(self.k + 1)]
        witness: List[Optional[int]] = [None] * (self.k + 1)
        b: List[Optional[int]] = [None] * self.n
        counts: Dict[int, int] = dict(collections.Counter(self.col_values))

        hint_values: Optional[List[int]] = None
        if hint_pi is not None:
            if len(hint_pi) != self.n:
                raise ValueError(f"hint_pi length must be n={self.n}.")
            hint_values = [self.col_values[int(hint_pi[j])] for j in range(self.n)]

        nodes = 0
        best_remaining_j = self.n

        def suffix_c(j: int) -> int:
            c = 0
            for ell in self.support_gt[j]:
                val = b[ell]
                if val is None:
                    raise RuntimeError(
                        f"Internal error: b[{ell}] needed before assignment at j={j}."
                    )
                c ^= int(val)
            return c

        def candidate_values_for_layer(i: int, c: int) -> Iterable[int]:
            # If r_j^*=0, then a_j must be zero, so b_j=c is forced.
            if i == 0:
                return [c]
            return [v for v, cnt in counts.items() if cnt > 0]

        def dfs(j: int) -> Optional[bool]:
            nonlocal nodes, best_remaining_j

            nodes += 1
            if nodes > node_limit:
                return None
            if time.time() - start > time_limit_sec:
                return None

            if j < 0:
                return True

            if j < best_remaining_j:
                best_remaining_j = j
                if self.verbose:
                    print(
                        f"reached j={j}, nodes={nodes}, elapsed={time.time()-start:.3f}s",
                        flush=True,
                    )

            c = suffix_c(j)
            i = self.rank_target[j]
            h_pivot = self.pivot_rank.get(j)
            hv = hint_values[j] if hint_values is not None else None

            candidates: List[Tuple[int, int, float, int, int]] = []
            for v in candidate_values_for_layer(i, c):
                if counts.get(v, 0) <= 0:
                    continue

                a = int(v) ^ c

                # Capacity pruning: a_j must fit in every K_m, m >= i.
                feasible = True
                expansion_score = 0
                for m in range(i, self.k + 1):
                    if not K[m].contains(a):
                        if K[m].dim >= m:
                            feasible = False
                            break
                        expansion_score += 1
                if not feasible:
                    continue

                # Pivot witness must not already lie in previous layer.
                if h_pivot is not None and K[h_pivot - 1].contains(a):
                    continue

                # Candidate ordering:
                #   1. hint match first, if hint is provided;
                #   2. consume less dimension budget first;
                #   3. randomized tie-break.
                hint_penalty = 0 if (hv is not None and v == hv) else 1
                candidates.append(
                    (hint_penalty, expansion_score, self.rng.random(), int(v), a)
                )

            candidates.sort()

            for _, _, _, v, a in candidates:
                # Snapshot state. This is simple and safe for n up to moderate sizes.
                K_old = [space.copy() for space in K]
                witness_old = list(witness)

                b[j] = v
                counts[v] -= 1

                for m in range(i, self.k + 1):
                    K[m].add(a)

                if h_pivot is not None:
                    witness[h_pivot] = a

                # A processed pivot p_h must stay outside K_{h-1}.
                valid = True
                for h in range(1, self.k + 1):
                    wh = witness[h]
                    if wh is not None and K[h - 1].contains(wh):
                        valid = False
                        break

                if valid:
                    out = dfs(j - 1)
                    if out is True:
                        return True
                    if out is None:
                        return None

                # Restore.
                for idx in range(self.k + 1):
                    K[idx] = K_old[idx]
                witness[:] = witness_old
                counts[v] += 1
                b[j] = None

            return False

        status_bool = dfs(self.n - 1)
        elapsed = time.time() - start

        if status_bool is True:
            assert all(x is not None for x in b)
            b_int = [int(x) for x in b]  # type: ignore[arg-type]
            pi = self._recover_pi_from_b(b_int)
            ok, pivots = self.verify_pi(pi)
            if not ok:
                return SearchResult(
                    status="INTERNAL_ERROR",
                    pi=pi,
                    pivot_set=pivots,
                    nodes=nodes,
                    elapsed_sec=elapsed,
                    best_remaining_j=best_remaining_j,
                    message="DFS reached a leaf, but verification failed.",
                )
            return SearchResult(
                status="FOUND",
                pi=pi,
                pivot_set=pivots,
                nodes=nodes,
                elapsed_sec=elapsed,
                best_remaining_j=best_remaining_j,
                message="Found a permutation whose pivot profile matches the target.",
            )

        if status_bool is None:
            return SearchResult(
                status="LIMIT",
                pi=None,
                pivot_set=None,
                nodes=nodes,
                elapsed_sec=elapsed,
                best_remaining_j=best_remaining_j,
                message="Stopped by time or node limit.",
            )

        return SearchResult(
            status="UNSAT_BY_SEARCH",
            pi=None,
            pivot_set=None,
            nodes=nodes,
            elapsed_sec=elapsed,
            best_remaining_j=best_remaining_j,
            message="Exhausted the search tree under the implemented pruning.",
        )


# ---------------------------------------------------------------------------
# Loading and CLI.
# ---------------------------------------------------------------------------

def load_target(data: np.lib.npyio.NpzFile, key: Optional[str]) -> np.ndarray:
    if key is not None:
        return np.asarray(data[key], dtype=int)
    for candidate in ("desired_p_star", "p_star", "S_star", "target", "target_pivots"):
        if candidate in data.files:
            return np.asarray(data[candidate], dtype=int)
    raise KeyError(
        "Could not find target pivot set. Pass --target-key explicitly."
    )


def load_hint_pi(data: np.lib.npyio.NpzFile, key: Optional[str]) -> Optional[List[int]]:
    if key is None:
        return None
    if key not in data.files:
        raise KeyError(f"hint key {key!r} not found in npz.")
    arr = data[key]
    if arr.ndim == 1:
        return [int(x) for x in arr]
    if arr.ndim == 2:
        return pi_from_permutation_matrix(arr)
    raise ValueError(f"Hint array {key!r} must be pi vector or permutation matrix.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("npz_path", help="Input .npz containing Gb, Gp, and target pivots.")
    parser.add_argument("--gb-key", default="Gb")
    parser.add_argument("--gp-key", default="Gp")
    parser.add_argument("--target-key", default=None)
    parser.add_argument(
        "--hint-key",
        default=None,
        help=(
            "Optional pi/P key used only for candidate ordering, e.g. pi_true. "
            "Do not use this when testing blind reconstruction."
        ),
    )
    parser.add_argument("--time-limit", type=float, default=60.0)
    parser.add_argument("--node-limit", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--save", default=None, help="Optional output .npz path.")
    args = parser.parse_args(argv)

    data = np.load(args.npz_path, allow_pickle=True)
    Gb = np.asarray(data[args.gb_key], dtype=np.uint8) & 1
    Gp = np.asarray(data[args.gp_key], dtype=np.uint8) & 1
    target = load_target(data, args.target_key)

    hint_pi = load_hint_pi(data, args.hint_key)

    solver = PivotBacktracker(Gb, Gp, target, seed=args.seed, verbose=args.verbose)

    if hint_pi is not None:
        ok, piv = solver.verify_pi(hint_pi)
        print(f"Hint verification: ok={ok}, pivot_set={piv}")

    result = solver.solve(
        time_limit_sec=args.time_limit,
        node_limit=args.node_limit,
        hint_pi=hint_pi,
    )

    print(json.dumps(dataclasses.asdict(result), indent=2))

    if args.save is not None:
        if result.pi is not None:
            pi = np.asarray(result.pi, dtype=np.int64)
            P = permutation_matrix_from_pi(pi, solver.n)
            piv = np.asarray(result.pivot_set, dtype=np.int64)
        else:
            pi = np.array([], dtype=np.int64)
            P = np.zeros((solver.n, solver.n), dtype=np.uint8)
            piv = np.array([], dtype=np.int64)

        np.savez(
            args.save,
            status=result.status,
            pi=pi,
            P=P,
            pivot_set=piv,
            target=np.asarray(target, dtype=np.int64),
            nodes=np.asarray(result.nodes, dtype=np.int64),
            elapsed_sec=np.asarray(result.elapsed_sec, dtype=float),
            best_remaining_j=np.asarray(result.best_remaining_j, dtype=np.int64),
            message=result.message,
        )
        print(f"Saved result to {args.save}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
