#!/usr/bin/env python3
"""galois-free extended primitive narrow-sense BCH generator, [2^m, k], coordinates beta_0=0, beta_j=alpha^(j-1).
Validated against project/POD/eBCH_m6_t11.matrix and eBCH_m7_t10.matrix (identical row space).
Conway polynomials: gl_shear_search_20260919.CONWAY."""
import numpy as np
import gl_shear_search_20260919 as G
from cpsat_pivot_reconstruct_duality import gf2_nullspace

def ebch(m: int, t: int) -> np.ndarray:
    n = 1 << m; red = G.CONWAY[m]
    pw = [0] * (n - 1); e = 1
    for i in range(n - 1):
        pw[i] = e; e <<= 1
        if (e >> m) & 1: e = (e ^ (1 << m)) ^ red
    log = {pw[i]: i for i in range(n - 1)}
    beta = [0] + pw
    def powf(x, s): return 0 if x == 0 else pw[(log[x] * s) % (n - 1)]
    rows = []
    for s in range(1, 2 * t, 2):            # odd s cover all conjugates of alpha^1..alpha^{2t}
        for b in range(m):
            rows.append([(powf(beta[j], s) >> b) & 1 for j in range(n)])
    rows.append([1] * n)                     # extension parity
    return G.gf2(gf2_nullspace(np.array(rows, np.uint8)))

def dump_cp(Gb, m, path):
    """CP[x] = column x of G_b[:, pi_field], k bits, for gl_flag_enum."""
    n = 1 << m; k = Gb.shape[0]; Cp = Gb[:, G.field_permutation(m)]
    with open(path, 'w') as f:
        f.write(f'{m} {k}\n')
        for x in range(n):
            v = sum(int(Cp[r, x]) << r for r in range(k))
            f.write(f'{v >> 64:016x} {v & ((1 << 64) - 1):016x}\n')
