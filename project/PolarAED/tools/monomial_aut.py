#!/usr/bin/env python3
"""
Affine automorphism structure of a polar code, from its information set.

Implements, in binary-expansion form:
  * [P21, Thm. 2]  admissibility of a variable change x_i <- x_j
  * [P22, Thm. 3]  the block structure S with Aut_aff = BLTA(S)
  * [P22, Lem. 14] |BLTA(S)|
  * [P22, Thm. 4 / Eq. 48]  the SC-absorption group [1] and the EC count
  * [P21, Eq. 4] / [P22, Alg. 1]  the auxiliary matrix A_G

Conventions are those of PLAN.md 1.2 (LSB-first; x_i present in m_r iff bit i
of r is 0), verified by tools/_feasibility_check.py.

Usage:
    python3 monomial_aut.py --I_file ../codes/<id>/<id>_I.txt --n 7
"""

import argparse

import numpy as np


# ------------------------------------------------------------- admissibility

def affected(n, I, i, j):
    """Indices produced by the variable change x_i <- x_j, for r in I."""
    return [((r | (1 << i)) & ~(1 << j)) for r in I if not ((r >> i) & 1)]


def admissible(n, I, i, j):
    """May A[i][j] be set to 1? ([P21, Thm. 2])"""
    if i == j:
        return True
    Iset = set(I)
    return all(t in Iset for t in affected(n, I, i, j))


def admissible_matrix(n, I):
    return [[admissible(n, I, i, j) for j in range(n)] for i in range(n)]


def is_decreasing(n, I, adm=None):
    """All lower-triangular positions admissible <=> UPO-compliant / decreasing."""
    adm = adm or admissible_matrix(n, I)
    return all(adm[i][j] for i in range(n) for j in range(i))


def AG_matrix(n, I):
    """A_G[i][j] = indices that must be ADDED to I to free position (i, j)."""
    Iset = set(I)
    return [[sorted(set(t for t in affected(n, I, i, j) if t not in Iset))
             if i != j else []
             for j in range(n)] for i in range(n)]


# --------------------------------------------------------- BLTA block structure

def block_structure(n, adm):
    """Partition [0, n) into the diagonal blocks implied by admissible upper
    entries: a block must contain every column reachable from its rows.
    """
    blocks = []
    start = 0
    while start < n:
        end = start
        i = start
        while i <= end:
            for j in range(i + 1, n):
                if adm[i][j] and j > end:
                    end = j
            i += 1
        blocks.append(end - start + 1)
        start = end + 1
    return blocks


def blta_is_exact(n, adm, S):
    """True if EVERY upper position inside a block is admissible, i.e. the
    automorphism group is exactly BLTA(S) rather than a proper subgroup.
    [P22, Thm. 3] says this holds for decreasing codes; check it anyway.
    """
    off = 0
    for s in S:
        for i in range(off, off + s):
            for j in range(i + 1, off + s):
                if not adm[i][j]:
                    return False
        off += s
    return True


def blta_size(n, S):
    """|BLTA(S)| = 2^{n(n+1)/2} * prod_i prod_{j=2..s_i} (2^j - 1)  [P22, Lem. 14]"""
    total = 1 << (n * (n + 1) // 2)
    for s in S:
        for j in range(2, s + 1):
            total *= (1 << j) - 1
    return total


def absorption_profile(n, S):
    """[1] as a BLTA profile. Working assumption of [P22]: BLTA(2,1,...,1) when
    the first block has size > 1, else plain LTA(n) = BLTA(1,...,1).
    """
    if S and S[0] > 1:
        return [2] + [1] * (n - 2)
    return [1] * n


def num_ec(n, S):
    """Number of equivalence classes E = |BLTA(S)| / |[1]|  ([P22, Eq. 48])."""
    return blta_size(n, S) // blta_size(n, absorption_profile(n, S))


def au_ap_sizes(S):
    """|A_U| and |A_P| of [P22, Eq. 14]."""
    au = 1
    ap = 1
    for s in S:
        au *= 1 << (s * (s - 1) // 2)
        for k in range(2, s + 1):
            ap *= k
    return au, ap


# ------------------------------------------------------------- matrix helpers

def inv_mod2(A):
    n = A.shape[0]
    M = np.concatenate([A.copy() % 2, np.eye(n, dtype=np.uint8)], axis=1)
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, n):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            raise ValueError("singular matrix")
        M[[r, piv]] = M[[piv, r]]
        for i in range(n):
            if i != r and M[i, c]:
                M[i, :] ^= M[r, :]
        r += 1
    return M[:, n:]


def in_blta(A, S):
    """Is A block-lower-triangular for profile S? (entries above blocks zero)"""
    n = A.shape[0]
    block_of = []
    for b, s in enumerate(S):
        block_of += [b] * s
    for i in range(n):
        for j in range(n):
            if A[i, j] and block_of[j] > block_of[i]:
                return False
    return True


def same_ec(A1, A2, n, S1):
    """[P22, Lem. 6]: pi_1 ~ pi_2  <=>  A1 . A2^{-1} in [1] = BLTA(S1)."""
    return in_blta((A1 @ inv_mod2(A2)) % 2, S1)


# ---------------------------------------------------------- report / CLI

def report(n, I):
    adm = admissible_matrix(n, I)
    dec = is_decreasing(n, I, adm)
    S = block_structure(n, adm)
    S1 = absorption_profile(n, S)
    exact = blta_is_exact(n, adm, S)
    ut = [(i, j) for i in range(n) for j in range(i + 1, n) if adm[i][j]]
    au, ap = au_ap_sizes(S)

    print("  n = %d, N = %d, K = %d" % (n, 1 << n, len(I)))
    print("  decreasing / UPO-compliant : %s%s"
          % (dec, "" if dec else "   <-- LTA is NOT contained in Aut!"))
    print("  admissible UT positions    : %d  %s" % (len(ut), ut))
    print("  block structure S          : %s" % (tuple(S),))
    print("  Aut_aff == BLTA(S) exactly : %s" % exact)
    print("  |BLTA(S)|                  : %d" % blta_size(n, S))
    print("  [1] profile (SC-absorbed)  : %s   |[1]| = %d"
          % (tuple(S1), blta_size(n, S1)))
    print("  equivalence classes E      : %d" % num_ec(n, S))
    print("  |A_U| = %d, |A_P| = %d, pure-UTL elements = %d"
          % (au, ap, 1 << len(ut)))
    for i in range(n):
        print("    " + "".join("D" if i == j else
                               ("1" if adm[i][j] else ".") for j in range(n)))
    return dict(adm=adm, S=S, S1=S1, ut=ut, decreasing=dec, exact=exact)


def load_I(path):
    with open(path) as f:
        return sorted(int(x) for x in f.read().split())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I_file", required=True)
    ap.add_argument("--n", type=int, required=True)
    args = ap.parse_args()
    report(args.n, load_I(args.I_file))


if __name__ == "__main__":
    main()
