#!/usr/bin/env python3
"""
pivot_to_cnf  --  export a pivot-profile instance as PLAIN DIMACS CNF for the
cube-and-conquer toolchain (march_cu + iglucose), which do NOT understand the
'x' XOR-clause extension.

Every native XOR of the reduced encoding is Tseitin-expanded to pure CNF via a
chain of 3-variable XORs:

    x1 ^ x2 ^ ... ^ xw = r
    -->  t2 = x1 ^ x2 ;  t3 = t2 ^ x3 ;  ... ;  t_w = t_{w-1} ^ xw ;  t_w = r

each 3-var XOR  a ^ b = c  is the 4 clauses
    (-a -b -c)(-a b c)(a -b c)(a b -c)
so a width-w XOR costs (w-2) aux vars and ~4(w-1) clauses -- linear, no blow-up.

Uses the SAME CnfXor object that cms_pivot_reduced_20260903 builds, so the
encoding (staged W, span, U-certificate, narrow/min-weight cuts) is identical
to what the CMS solver sees; only the file format differs.
"""
from __future__ import annotations

import argparse
import os
from typing import List, Tuple

import numpy as np

from cpsat_pivot_reconstruct_duality import gf2, polar_matrix
from cpsat_pivot_crypto_minixor import CnfXor
from cms_pivot_reduced_20260903 import encode_reduced


def _tseitin_xor_chain(cnf_clauses: List[List[int]], next_var: int,
                       lits: List[int], rhs: bool) -> int:
    """Append CNF for  XOR(lits) == rhs.  Returns the new next_var."""
    def xor3(a: int, b: int, c: int) -> None:            # a ^ b = c
        cnf_clauses.append([-a, -b, -c])
        cnf_clauses.append([-a, b, c])
        cnf_clauses.append([a, -b, c])
        cnf_clauses.append([a, b, -c])

    vs = list(lits)
    if not vs:
        if rhs:
            cnf_clauses.append([])                       # empty clause = UNSAT
        return next_var
    if len(vs) == 1:
        cnf_clauses.append([vs[0] if rhs else -vs[0]])
        return next_var
    acc = vs[0]
    for i in range(1, len(vs) - 1):
        t = next_var; next_var += 1
        xor3(acc, vs[i], t)
        acc = t
    # last: acc ^ vs[-1] = rhs  <=>  acc ^ vs[-1] ^ (rhs as unit) ...
    last = vs[-1]
    if rhs:
        xor3(acc, last, _TRUE_SENTINEL[0])
    else:
        # acc ^ last = 0  ->  acc == last
        cnf_clauses.append([-acc, last])
        cnf_clauses.append([acc, -last])
    return next_var


_TRUE_SENTINEL = [0]   # filled in write_plain_dimacs


def write_plain_dimacs(cnf: CnfXor, path: str) -> Tuple[int, int]:
    """Expand cnf.xors to CNF and write plain DIMACS.  Returns (nvars, nclauses)."""
    clauses: List[List[int]] = [list(c) for c in cnf.clauses]
    nv = cnf.nvars

    # one always-true var for odd-parity (rhs=1) closes
    nv += 1
    true_v = nv
    _TRUE_SENTINEL[0] = true_v
    clauses.append([true_v])

    for vs, rhs in cnf.xors:
        nv = _tseitin_xor_chain(clauses, nv, vs, rhs)

    # drop trivially-satisfied / normalise
    with open(path, "w") as f:
        f.write(f"p cnf {nv} {len(clauses)}\n")
        for c in clauses:
            f.write(" ".join(map(str, c)) + " 0\n")
    return nv, len(clauses)


def parse_mat(p: str) -> np.ndarray:
    tok = open(p).read().split()
    r, c = int(tok[0]), int(tok[1]); body = tok[2:]
    if len(body) == r:
        return np.array([[int(x) for x in row] for row in body], np.uint8)
    return np.array([int(x) for x in body], np.uint8).reshape(r, c)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--p-star", required=True, help="comma-separated pivot indices")
    ap.add_argument("--certificate", choices=["span-only", "full"], default="full")
    ap.add_argument("--min-weight-words", type=int, default=0)
    ap.add_argument("--kernels", default=None)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    Gb = gf2(parse_mat(args.gen))
    Gp = gf2(polar_matrix(args.m))
    p_star = [int(x) for x in args.p_star.split(",")]

    cnf, stats = encode_reduced(
        Gb, Gp, p_star, kernel_string=args.kernels,
        certificate=args.certificate, min_weight_words=args.min_weight_words,
        verbose=True)
    nv, nc = write_plain_dimacs(cnf, args.out)
    print(f"[pivot_to_cnf] xors={len(cnf.xors)}  base clauses={len(cnf.clauses)}")
    print(f"[pivot_to_cnf] wrote {args.out}:  {nv} vars, {nc} clauses")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
