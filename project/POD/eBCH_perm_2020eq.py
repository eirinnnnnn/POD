import numpy as np
import argparse
import json
from generator_set_eBCH import generate_eBCH, save_matrix_txt

def build_P_proposed(m: int, t: int):
    """
    Construct the permutation matrix P exactly as in
    'Soft Polar Decoding of BCH Codes, Golay Codes and LDPC Codes' (2020).

    Paper convention:
        coordinate order = [α^0, α^1, ..., α^(n-2), 0]
        P_{i,j} = 1 if
            i = n-1 - α^j,  0 <= j <= n-2
            i = j = n-1
        P_{i,j} = 0 otherwise

    The returned permutation matrix P is expressed in primitive-sense canonical
    coordinate indexing (i.e., compatible with G from generate_eBCH).
    """

    # ------------------------------------------------------------------
    # Step 0. Generate eBCH and basic objects
    # ------------------------------------------------------------------
    G, H, elements, GF = generate_eBCH(m, t)
    n = 2 ** m
    alpha = GF.primitive_element

    # Map: integer label of field element -> your coordinate index
    idx_of_int = {int(elements[i]): i for i in range(n)}

    # ------------------------------------------------------------------
    # Step 1. Paper coordinate order
    #   index 0..n-2 : α^0, α^1, ..., α^(n-2)
    #   index n-1    : 0
    # ------------------------------------------------------------------
    paper_coords = [alpha ** j for j in range(n - 1)] + [GF(0)]

    # ------------------------------------------------------------------
    # Step 2. Build permutation in PAPER index space
    # ------------------------------------------------------------------
    perm_paper = [None] * n

    for j in range(n - 1):
        aj_int = int(alpha ** j)   # integer label in {1,...,n-1}
        i = (n - 1) - aj_int       # NATURAL-NUMBER subtraction
        perm_paper[j] = i

    # Special zero coordinate
    perm_paper[n - 1] = n - 1

    # Sanity check: must be bijection
    if len(set(perm_paper)) != n:
        raise RuntimeError(
            f"paper permutation not bijective: {perm_paper}"
        )

    # ------------------------------------------------------------------
    # Step 3. Convert paper permutation to YOUR coordinate indexing
    # ------------------------------------------------------------------
    perm = [None] * n

    for j_paper, i_paper in enumerate(perm_paper):
        j_elem_int = int(paper_coords[j_paper])
        i_elem_int = int(paper_coords[i_paper])

        j_idx = idx_of_int[j_elem_int]
        i_idx = idx_of_int[i_elem_int]

        perm[j_idx] = i_idx

    if any(v is None for v in perm):
        raise RuntimeError("Permutation incomplete after index conversion.")

    # ------------------------------------------------------------------
    # Step 4. Build permutation matrix P in your indexing
    # ------------------------------------------------------------------
    P = np.zeros((n, n), dtype=np.uint8)
    for j in range(n):
        P[perm[j], j] = 1

    # Final sanity check
    if not (np.all(P.sum(axis=0) == 1) and np.all(P.sum(axis=1) == 1)):
        raise RuntimeError("Constructed P is not a valid permutation matrix.")

    return perm, P

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--t", type=int, required=True)
    ap.add_argument("--out_prefix", type=str, default=None,
                    help="prefix for outputs (default: eBCH_m{m}_t{t}_P)")
    args = ap.parse_args()

    m, t = args.m, args.t
    n = 2**m
    out_prefix = args.out_prefix or f"eBCH_m{m}_t{t}_P"

    perm, P = build_P_proposed(m, t)

    # Save permutation vector
    with open(f"{out_prefix}_perm.json", "w") as f:
        json.dump(
            {
                "m": m,
                "t": t,
                "n": n,
                "definition": "paper-style: perm[j]=i where P_{i,j}=1, mapped to canonical GF-indexed coordinates",
                "perm": perm,
            },
            f,
            indent=2,
        )

    # Save matrix in .matrix format
    save_matrix_txt(f"{out_prefix}.matrix", P)

    print(f"[ok] n={n}, saved:")
    print(f"  - {out_prefix}_perm.json")
    print(f"  - {out_prefix}.matrix")

if __name__ == "__main__":
    main()