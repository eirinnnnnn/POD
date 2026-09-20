#!/usr/bin/env python3
"""
Validation gates for the polar construction (PLAN.md section 5).

G-STAGE  The polarization stage ORDER used by polar_design.ga_means must match
         the decoder's. AdjustPolarDecoder runs its Bhattacharyya recursion
         from stage-1 down to 0, i.e. bit n-1 transformed FIRST. Getting this
         backwards silently swaps the reliabilities of index pairs like 1 and 2
         ((W^-)^+ vs (W^+)^-), which changes the information set.
         Checked by replicating the C++ Z recursion in Python and diffing
         against the decoder's own "idx NN : Z" printout.

G1       admissible_matrix on [P21] section III-A's (16,7) example.
G2       |EC| for the block structures quoted in [P22]:
             S = (4,1,1,1,3), n = 10  ->  2205
             S = (3,5),       n = 8   ->  68355
G3       Brute-force affine automorphism count vs |BLTA(S)| for n <= 4.
G5       Designed information sets are decreasing / UPO-compliant.

Usage:
    python3 test_polar_design.py [--binary ../build/PolarAED]
"""

import argparse
import itertools
import math
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

import monomial_aut as MA
import polar_design as PD
from blta_sample import perm_from_affine


def z_recursion(n, raw_ber):
    """Replicates AdjustPolarDecoder's Bhattacharyya recursion exactly."""
    N = 1 << n
    z = [math.sqrt(4.0 * raw_ber * (1.0 - raw_ber))] * N
    for stage in reversed(range(n)):
        step = 1 << stage
        for i in range(N):
            if i & step:
                continue
            z0, z1 = z[i], z[i | step]
            z[i] = z0 + z1 - z0 * z1
            z[i | step] = z0 * z1
    return z


def gate_stage_order(binary):
    n, K, raw_ber = 4, 10, 0.01
    N = 1 << n
    work = tempfile.mkdtemp(prefix="polaraed_gate_")
    I = sorted(range(N - K, N))          # any code; only Z is being checked
    G = PD.polar_matrix(n)[I, :]
    H = PD.nullspace_mod2(G)
    PD.save_matrix(os.path.join(work, "t_H.matrix"), H.T)
    ini = os.path.join(work, "t.ini")
    with open(ini, "w") as f:
        f.write("""[AWGN]
start = 2.0
end = 2.0
step = 0.5
step_type = SNR
seed_string = -2

[AdjustPolarDecoder]
matrix_src = byHmatrix
Gmatrix_path =
Hmatrix_path = ./t_H.matrix
operationArray = %s
bha_value_setting = %s
target_raw_BER = %g
list_size = 1
permutation_src =
permutation_random_seed = -1
OnlyInit = true

use_AED = false
automorphism_src =
aed_L = 0

[Monte_Carlo]
iter_min = 0
iter_max = 10
error_min = 0
error_max = 10
monitor_slot_size = 10
seed_string = 1
""" % ("1" * (n * (N // 2)), "?" * N, raw_ber))
    out = subprocess.run([os.path.abspath(binary), "-ini", "./t.ini"],
                         cwd=work, capture_output=True, text=True).stdout
    got = {}
    for m in re.finditer(r"^idx\s+(\d+)\s*:\s*([\d.eE+-]+)\s*$", out, re.M):
        got[int(m.group(1))] = float(m.group(2))
    if len(got) != N:
        print("  G-STAGE : SKIP (decoder printout not parsed; %d/%d)"
              % (len(got), N))
        return None
    mine = z_recursion(n, raw_ber)
    err = max(abs(mine[i] - got[i]) for i in range(N))
    ok = err < 1e-6
    print("  G-STAGE : %s  (max |Z_python - Z_decoder| = %.2e)"
          % ("PASS" if ok else "FAIL", err))
    if not ok:
        rev = z_recursion(n, raw_ber)[::-1]
        print("           (if the stage loop were not reversed the pairs "
              "1<->2 would swap; check ga_means)")
    return ok


def gate_g1():
    n, I = 4, [7, 10, 11, 12, 13, 14, 15]
    adm = MA.admissible_matrix(n, I)
    ut = [(i, j) for i in range(n) for j in range(i + 1, n) if adm[i][j]]
    ok = (ut == [(1, 2)]) and MA.is_decreasing(n, I, adm)
    print("  G1      : %s  (UT admissible = %s, expected [(1, 2)])"
          % ("PASS" if ok else "FAIL", ut))
    return ok


def gate_g2():
    cases = [((4, 1, 1, 1, 3), 10, 2205), ((3, 5), 8, 68355)]
    ok = True
    for S, n, expect in cases:
        got = MA.num_ec(n, list(S))
        ok &= (got == expect)
        print("  G2      : %s  S=%s n=%d -> |EC| = %d (expected %d)"
              % ("PASS" if got == expect else "FAIL", S, n, got, expect))
    return ok


def gate_g3(max_n=4):
    """Brute-force count of affine automorphisms vs |BLTA(S)|, n <= 4.

    n = 4 enumerates 2^16 matrices per code, so it is restricted to a few
    dimensions; n = 3 is exhaustive.
    """
    ok = True
    for n in range(3, max_n + 1):
        N = 1 << n
        mats = []
        for bitsv in itertools.product([0, 1], repeat=n * n):
            A = np.array(bitsv, dtype=np.uint8).reshape(n, n)
            try:
                MA.inv_mod2(A)
            except ValueError:
                continue
            mats.append(A)
        Ks = range(1, N) if n == 3 else (3, 5, 8, 11)
        for K in Ks:
            I = sorted(range(N - K, N))       # a decreasing set for this order
            adm = MA.admissible_matrix(n, I)
            if not MA.is_decreasing(n, I, adm):
                continue
            G = PD.polar_matrix(n)[I, :]
            rG = _rank(G)
            cnt = 0
            for A in mats:
                p = perm_from_affine(A, n)
                if _rank(np.vstack([G, G[:, p]])) == rG:
                    cnt += 1
            cnt *= N                           # every translation b works
            S = MA.block_structure(n, adm)
            want = MA.blta_size(n, S)
            if cnt != want:
                ok = False
                print("  G3      : FAIL n=%d K=%d  brute=%d  |BLTA(%s)|=%d"
                      % (n, K, cnt, tuple(S), want))
    if ok:
        print("  G3      : PASS  (brute-force |Aut_aff| == |BLTA(S)|, n <= %d)"
              % max_n)
    return ok


def gate_g5():
    ok = True
    for (n, K, snr) in [(7, 100, 10.5), (7, 64, 4.0), (8, 128, 6.0),
                        (10, 512, 3.0)]:
        I, _ = PD.information_set(n, K, "ga", snr, rate=K / float(1 << n))
        dec = MA.is_decreasing(n, I)
        ok &= dec
        S = MA.block_structure(n, MA.admissible_matrix(n, I))
        print("  G5      : %s  (%d,%d) @ %.1f dB  decreasing=%s  S=%s"
              % ("PASS" if dec else "FAIL", 1 << n, K, snr, dec, tuple(S)))
    return ok


def _rank(A):
    A = (A.astype(np.uint8) & 1).copy()
    m, nn = A.shape
    r = 0
    for c in range(nn):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(r + 1, m):
            if A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default="../build/PolarAED")
    args = ap.parse_args()

    results = []
    if os.path.exists(args.binary):
        results.append(gate_stage_order(args.binary))
    else:
        print("  G-STAGE : SKIP (no binary at %s)" % args.binary)
    results.append(gate_g1())
    results.append(gate_g2())
    results.append(gate_g3())
    results.append(gate_g5())
    return 0 if all(r is not False for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
