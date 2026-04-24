# import json
# import numpy as np
# import os
# from itertools import permutations as it_permutations

# import galois
# import numpy as np


# # ================================================================
# #  Linear algebra utilities (nullspace, rank) over GF(2)
# # ================================================================

# def nullspace_mod2(G):
#     """
#     Basis of right nullspace of G over GF(2): { x in F2^n : G x^T = 0 }.
#     Returns H whose rows form a basis.
#     """
#     G = (G.astype(np.uint8) & 1)
#     k, n = G.shape

#     A = G.copy()
#     m, n = A.shape
#     r = 0
#     pivots = []
#     for c in range(n):
#         if r == m:
#             break
#         piv = None
#         for i in range(r, m):
#             if A[i, c]:
#                 piv = i
#                 break
#         if piv is None:
#             continue
#         if piv != r:
#             A[[r, piv]] = A[[piv, r]]
#         for i in range(m):
#             if i != r and A[i, c]:
#                 A[i, :] ^= A[r, :]
#         pivots.append(c)
#         r += 1

#     pivot_set = set(pivots)
#     free_cols = [j for j in range(n) if j not in pivot_set]

#     if not free_cols:
#         return np.zeros((0, n), dtype=np.uint8)

#     pivot_col_to_row = {c: i for i, c in enumerate(pivots)}
#     H_rows = []
#     for free in free_cols:
#         v = np.zeros(n, dtype=np.uint8)
#         v[free] = 1
#         for p in pivots:
#             row = pivot_col_to_row[p]
#             v[p] = A[row, free]
#         H_rows.append(v)

#     return np.vstack(H_rows)


# def rank_mod2(A):
#     A = (A.astype(np.uint8) & 1).copy()
#     m, n = A.shape
#     r = 0
#     for c in range(n):
#         if r == m:
#             break
#         piv = None
#         for i in range(r, m):
#             if A[i, c]:
#                 piv = i
#                 break
#         if piv is None:
#             continue
#         if piv != r:
#             A[[r, piv]] = A[[piv, r]]
#         for i in range(r + 1, m):
#             if A[i, c]:
#                 A[i, :] ^= A[r, :]
#         r += 1
#     return r


# # ================================================================
# #  Canonical eBCH generation with GF-indexed coordinates
# # ================================================================

# def generate_eBCH(m, t):
#     """
#     Canonical construction of extended primitive narrow-sense BCH over GF(2):

#       - length q = 2^m
#       - designed distance δ >= 2t+1
#       - coordinates = β_j, where:
#             β_0 = 0, β_1 = 1, β_2 = α, ..., β_{q-1} = α^{q-2}

#     Construction steps:
#       1) Build GF(2^m)-valued parity-check H_q using exponent set S = {1,...,2t}
#       2) Expand H_q to binary H_bin using the polynomial basis (default representation)
#       3) Add one all-ones parity row for the extension
#       4) Generator G_bin = nullspace(H_bin)

#     Returns:
#         G_bin, H_bin, elements, GF
#     """
#     GF = galois.GF(2**m)
#     alpha = GF.primitive_element

#     q = 2**m                      # length of extended code
#     elements = [GF(0)] + [alpha**(j-1) for j in range(1, q)]  # β_j
#     S = list(range(1, 2*t + 1))   # BCH exponent set = {1,...,2t}

#     # ---------------------------- Step 1: H_q over GF(2^m) ----------------------------
#     # H_q = GF.Zeros(len(S), q)
#     H_q = GF.Zeros((len(S), q))
#     for row, ell in enumerate(S):
#         H_q[row, 0] = GF(0)
#         for j in range(1, q):
#             H_q[row, j] = elements[j] ** ell

#     # ---------------------------- Step 2: expand to binary H_bin ----------------------
#     H_rows_bin = []
#     for row in range(H_q.shape[0]):
#         for bit in range(m):
#             v = np.zeros(q, dtype=np.uint8)
#             for j in range(q):
#                 coeffs = int(H_q[row, j])   # GF element → integer → bit decomposition
#                 v[j] = (coeffs >> bit) & 1
#             H_rows_bin.append(v)

#     H_bin = np.stack(H_rows_bin, axis=0)

#     # ---------------------------- Step 3: add extension parity row --------------------
#     parity_row = np.ones(q, dtype=np.uint8)
#     H_bin = np.vstack([parity_row, H_bin])

#     # ---------------------------- Step 4: generator = nullspace -----------------------
#     G_bin = nullspace_mod2(H_bin)

#     return G_bin, H_bin, elements, GF


# # ================================================================
# #  Automorphism check using parity-check H
# # ================================================================

# def is_automorphism(H, perm):
#     """
#     Check if a permutation 'perm' is an automorphism of the code with parity-check H.
#     Condition:
#         RowSpan(H) = RowSpan(H[:, perm])
#     ⇔      rank([H ; H_perm]) = rank(H)
#     """
#     H = (H.astype(np.uint8) & 1)
#     H_perm = H[:, perm]

#     rH = rank_mod2(H)
#     stacked = np.vstack([H, H_perm])
#     r_stacked = rank_mod2(stacked)

#     return r_stacked == rH


# def test_affine_semilinear(H, elements, GF):
#     """
#     Full correctness test:
#         - test all translations b in GF(2^m)
#         - primitive multiplier αx
#         - Frobenius x^2
#     """
#     n = len(elements)
#     idx = {int(elements[i]): i for i in range(n)}
#     failures = []

#     # translations
#     for b in elements:
#         perm = [idx[int(a + b)] for a in elements]
#         if not is_automorphism(H, perm):
#             failures.append(("translation", int(b)))

#     # primitive multiplier
#     alpha = GF.primitive_element
#     perm_mult = [idx[int(alpha * x)] for x in elements]
#     if not is_automorphism(H, perm_mult):
#         failures.append(("multiplier", int(alpha)))

#     # Frobenius
#     m = GF.degree
#     if m > 1:
#         perm_frob = [idx[int(x**2)] for x in elements]
#         if not is_automorphism(H, perm_frob):
#             failures.append(("frobenius", 1))

#     return failures


# # ================================================================
# #  Automorphism-generators: multiplier + Frobenius + translation basis
# # ================================================================

# def aut_generators_multipliers_frobenius(elements, GF, H=None, verify=True):
#     """
#     Construct a generating set of the automorphism group:
#       - identity
#       - primitive multiplier x -> α x
#       - Frobenius x -> x^2
#       - translations by a basis of GF(2^m) over GF(2)

#     If verify=True, each candidate is tested with is_automorphism(H, perm).
#     """
#     n = len(elements)
#     idx = {int(elements[i]): i for i in range(n)}
#     gens = []

#     # identity
#     id_perm = list(range(n))
#     if (not verify) or (H is None) or is_automorphism(H, id_perm):
#         gens.append(("identity", None, id_perm))

#     # primitive multiplier
#     alpha = GF.primitive_element
#     perm_mult = [idx[int(alpha * x)] for x in elements]
#     if (not verify) or (H is None) or is_automorphism(H, perm_mult):
#         gens.append(("multiplier", int(alpha), perm_mult))

#     # Frobenius
#     m = GF.degree
#     if m > 1:
#         perm_frob = [idx[int(x**2)] for x in elements]
#         if (not verify) or (H is None) or is_automorphism(H, perm_frob):
#             gens.append(("frobenius", 1, perm_frob))

#     # translation basis = {1, α, α^2, ..., α^{m-1}}
#     basis = [GF(1)] + [alpha**i for i in range(1, m)]
#     for b in basis:
#         perm_trans = [idx[int(a + b)] for a in elements]
#         if (not verify) or (H is None) or is_automorphism(H, perm_trans):
#             gens.append(("translation", int(b), perm_trans))

#     return gens


# # ================================================================
# #  JSON utility
# # ================================================================

# def save_checkpoint_json(path, data):
#     tmp = path + ".tmp"
#     with open(tmp, "w") as f:
#         json.dump(data, f, indent=2)
#     os.replace(tmp, path)
#     print(f"[JSON saved → {path}]")

# def save_matrix_txt(path, M):
#     """
#     Save a binary matrix M (numpy array) as:
#       rows cols
#       <row 0 entries separated by spaces>
#       ...
#     """
#     M = (M.astype(np.uint8) & 1)
#     r, c = M.shape
#     with open(path, "w") as f:
#         f.write(f"{r} {c}\n")
#         for i in range(r):
#             f.write(" ".join(str(int(x)) for x in M[i, :]) + "\n")


# # ================================================================
# #  MAIN
# # ================================================================

# if __name__ == "__main__":
#     m = 4 
#     t = 2 

#     print(f"[info] Building eBCH(m={m}, t={t}) with canonical GF-indexing...")
#     G, H, elements, GF = generate_eBCH(m, t)

#     # Export G in the desired .matrix format for AdjustPolarDecoder
#     codetype = f"eBCH_m{m}_t{t}"
#     G_path = f"{codetype}.matrix"
#     save_matrix_txt(G_path, G)
#     print(f"[info] Saved G matrix -> {G_path}")

#     print("G shape:", G.shape)
#     print("H shape:", H.shape)

#     # build automorphism generator set
#     aut_gens = aut_generators_multipliers_frobenius(elements, GF, H=H, verify=True)

#     # sanity-check full affine/semilinear family
#     failures = test_affine_semilinear(H, elements, GF)
#     print("Failures:", failures)

#     # export for Schreier–Sims pipeline
#     H_perms = [perm for (_, _, perm) in aut_gens]

#     data = {
#         "n": G.shape[1],
#         "k": G.shape[0],
#         "|H|": len(H_perms),
#         "H_perms": H_perms,
#         "generator_records": [
#             {"type": typ, "meta": meta}
#             for (typ, meta, _) in aut_gens
#         ],
#     }

#     out_path = f"generator_set_eBCH_m{m}_t{t}.json"
#     save_checkpoint_json(out_path, data)
#     print(f"[info] Saved automorphism generator set -> {out_path}")

import json
import numpy as np
import os
from itertools import permutations as it_permutations

import galois
import numpy as np


# ================================================================
#  Linear algebra utilities (nullspace, rank) over GF(2)
# ================================================================

def nullspace_mod2(G):
    """
    Basis of right nullspace of G over GF(2): { x in F2^n : G x^T = 0 }.
    Returns H whose rows form a basis.
    """
    G = (G.astype(np.uint8) & 1)
    k, n = G.shape

    A = G.copy()
    m, n = A.shape
    r = 0
    pivots = []
    for c in range(n):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        pivots.append(c)
        r += 1

    pivot_set = set(pivots)
    free_cols = [j for j in range(n) if j not in pivot_set]

    if not free_cols:
        return np.zeros((0, n), dtype=np.uint8)

    pivot_col_to_row = {c: i for i, c in enumerate(pivots)}
    H_rows = []
    for free in free_cols:
        v = np.zeros(n, dtype=np.uint8)
        v[free] = 1
        for p in pivots:
            row = pivot_col_to_row[p]
            v[p] = A[row, free]
        H_rows.append(v)

    return np.vstack(H_rows)


def rank_mod2(A):
    A = (A.astype(np.uint8) & 1).copy()
    m, n = A.shape
    r = 0
    for c in range(n):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        for i in range(r + 1, m):
            if A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
    return r


def reduce_mod2(H):
    """
    Row-reduce H over GF(2) and return only the non-zero (linearly independent) rows.
    This gives a reduced parity-check matrix with exactly rank(H) rows.
    """
    A = (H.astype(np.uint8) & 1).copy()
    m, n = A.shape
    r = 0
    for c in range(n):
        if r == m:
            break
        piv = None
        for i in range(r, m):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i, :] ^= A[r, :]
        r += 1
    return A[:r, :]


# ================================================================
#  Canonical eBCH generation with GF-indexed coordinates
# ================================================================

def generate_eBCH(m, t):
    """
    Canonical construction of extended primitive narrow-sense BCH over GF(2):

      - length q = 2^m
      - designed distance δ >= 2t+1
      - coordinates = β_j, where:
            β_0 = 0, β_1 = 1, β_2 = α, ..., β_{q-1} = α^{q-2}

    Construction steps:
      1) Build GF(2^m)-valued parity-check H_q using exponent set S = {1,...,2t}
      2) Expand H_q to binary H_bin using the polynomial basis (default representation)
      3) Add one all-ones parity row for the extension
      4) Generator G_bin = nullspace(H_bin)

    Returns:
        G_bin, H_bin, elements, GF
    """
    GF = galois.GF(2**m)
    alpha = GF.primitive_element

    q = 2**m                      # length of extended code
    elements = [GF(0)] + [alpha**(j-1) for j in range(1, q)]  # β_j
    S = list(range(1, 2*t + 1))   # BCH exponent set = {1,...,2t}

    # ---------------------------- Step 1: H_q over GF(2^m) ----------------------------
    H_q = GF.Zeros((len(S), q))
    for row, ell in enumerate(S):
        H_q[row, 0] = GF(0)
        for j in range(1, q):
            H_q[row, j] = elements[j] ** ell

    # ---------------------------- Step 2: expand to binary H_bin ----------------------
    H_rows_bin = []
    for row in range(H_q.shape[0]):
        for bit in range(m):
            v = np.zeros(q, dtype=np.uint8)
            for j in range(q):
                coeffs = int(H_q[row, j])   # GF element → integer → bit decomposition
                v[j] = (coeffs >> bit) & 1
            H_rows_bin.append(v)

    H_bin = np.stack(H_rows_bin, axis=0)

    # ---------------------------- Step 3: add extension parity row --------------------
    parity_row = np.ones(q, dtype=np.uint8)
    H_bin = np.vstack([parity_row, H_bin])

    # ---------------------------- Step 4: reduce H_bin --------------------------------
    H_bin = reduce_mod2(H_bin)

    # ---------------------------- Step 5: generator = nullspace -----------------------
    G_bin = nullspace_mod2(H_bin)

    return G_bin, H_bin, elements, GF


# ================================================================
#  Automorphism check using parity-check H
# ================================================================

def is_automorphism(H, perm):
    """
    Check if a permutation 'perm' is an automorphism of the code with parity-check H.
    Condition:
        RowSpan(H) = RowSpan(H[:, perm])
    ⇔      rank([H ; H_perm]) = rank(H)
    """
    H = (H.astype(np.uint8) & 1)
    H_perm = H[:, perm]

    rH = rank_mod2(H)
    stacked = np.vstack([H, H_perm])
    r_stacked = rank_mod2(stacked)

    return r_stacked == rH


def test_affine_semilinear(H, elements, GF):
    """
    Full correctness test:
        - test all translations b in GF(2^m)
        - primitive multiplier αx
        - Frobenius x^2
    """
    n = len(elements)
    idx = {int(elements[i]): i for i in range(n)}
    failures = []

    # translations
    for b in elements:
        perm = [idx[int(a + b)] for a in elements]
        if not is_automorphism(H, perm):
            failures.append(("translation", int(b)))

    # primitive multiplier
    alpha = GF.primitive_element
    perm_mult = [idx[int(alpha * x)] for x in elements]
    if not is_automorphism(H, perm_mult):
        failures.append(("multiplier", int(alpha)))

    # Frobenius
    m = GF.degree
    if m > 1:
        perm_frob = [idx[int(x**2)] for x in elements]
        if not is_automorphism(H, perm_frob):
            failures.append(("frobenius", 1))

    return failures


# ================================================================
#  Automorphism-generators: multiplier + Frobenius + translation basis
# ================================================================

def aut_generators_multipliers_frobenius(elements, GF, H=None, verify=True):
    """
    Construct a generating set of the automorphism group:
      - identity
      - primitive multiplier x -> α x
      - Frobenius x -> x^2
      - translations by a basis of GF(2^m) over GF(2)

    If verify=True, each candidate is tested with is_automorphism(H, perm).
    """
    n = len(elements)
    idx = {int(elements[i]): i for i in range(n)}
    gens = []

    # identity
    id_perm = list(range(n))
    if (not verify) or (H is None) or is_automorphism(H, id_perm):
        gens.append(("identity", None, id_perm))

    # primitive multiplier
    alpha = GF.primitive_element
    perm_mult = [idx[int(alpha * x)] for x in elements]
    if (not verify) or (H is None) or is_automorphism(H, perm_mult):
        gens.append(("multiplier", int(alpha), perm_mult))

    # Frobenius
    m = GF.degree
    if m > 1:
        perm_frob = [idx[int(x**2)] for x in elements]
        if (not verify) or (H is None) or is_automorphism(H, perm_frob):
            gens.append(("frobenius", 1, perm_frob))

    # translation basis = {1, α, α^2, ..., α^{m-1}}
    basis = [GF(1)] + [alpha**i for i in range(1, m)]
    for b in basis:
        perm_trans = [idx[int(a + b)] for a in elements]
        if (not verify) or (H is None) or is_automorphism(H, perm_trans):
            gens.append(("translation", int(b), perm_trans))

    return gens


# ================================================================
#  JSON utility
# ================================================================

def save_checkpoint_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, path)
    print(f"[JSON saved → {path}]")

def save_matrix_txt(path, M):
    """
    Save a binary matrix M (numpy array) as:
      rows cols
      <row 0 entries separated by spaces>
      ...
    """
    M = (M.astype(np.uint8) & 1)
    r, c = M.shape
    with open(path, "w") as f:
        f.write(f"{r} {c}\n")
        for i in range(r):
            f.write(" ".join(str(int(x)) for x in M[i, :]) + "\n")


# ================================================================
#  MAIN
# ================================================================

if __name__ == "__main__":
    m = 7 
    t = 10 

    print(f"[info] Building eBCH(m={m}, t={t}) with canonical GF-indexing...")
    G, H, elements, GF = generate_eBCH(m, t)

    codetype = f"eBCH_m{m}_t{t}"

    # Export G
    G_path = f"{codetype}.matrix"
    save_matrix_txt(G_path, G)
    print(f"[info] Saved G matrix -> {G_path}  (shape {G.shape})")

    # Export reduced H
    H_path = f"{codetype}_H.matrix"
    save_matrix_txt(H_path, H.T)
    print(f"[info] Saved H matrix -> {H_path}  (shape {H.shape})")

    print("G shape:", G.shape)
    print("H shape:", H.shape)

    # Verify G H^T = 0
    check = (G @ H.T) % 2
    print(f"[check] G · H^T = 0 ? {np.all(check == 0)}")

    # build automorphism generator set
    aut_gens = aut_generators_multipliers_frobenius(elements, GF, H=H, verify=True)

    # sanity-check full affine/semilinear family
    failures = test_affine_semilinear(H, elements, GF)
    print("Failures:", failures)

    # export for Schreier–Sims pipeline
    H_perms = [perm for (_, _, perm) in aut_gens]

    data = {
        "n": G.shape[1],
        "k": G.shape[0],
        "|H|": len(H_perms),
        "H_perms": H_perms,
        "generator_records": [
            {"type": typ, "meta": meta}
            for (typ, meta, _) in aut_gens
        ],
    }

    out_path = f"generator_set_eBCH_m{m}_t{t}.json"
    save_checkpoint_json(out_path, data)
    print(f"[info] Saved automorphism generator set -> {out_path}")