#!/usr/bin/env python3
"""
generator_set_eGolay24.py

What this script does
---------------------
1) Embeds YOUR [24,12,8] extended Golay generator matrix G (12x24, systematic).
2) Builds the Steiner system S(5,8,24) "octads" from YOUR code by enumerating all 2^12 codewords,
   extracting supports of weight 8.
3) Builds the ATLAS-labeled Steiner octads by orbiting a seed octad under the ATLAS standard
   generators b11, b21 of M24 (in ATLAS coordinate labeling).
4) Finds a relabeling permutation s (length 24, 0-based) such that:
      YOUR_octads = { s(A) : A in ATLAS_octads }
   i.e. the Steiner systems match under coordinate relabeling.
5) Outputs conjugated generators for YOUR labeling:
      g' = s^{-1} g s
   which should pass the automorphism check on YOUR embedded G.
6) Stores JSON in the same style as generator_set_eBCH.py: {n,k,|H|,H_perms,...}
   and includes an automorphism-check section.

Notes
-----
- ATLAS generators b11, b21 generate M24 in a specific coordinate labeling (Steiner system labeling).
  Your G uses some (unknown) labeling, hence we compute the conjugating relabeling s.
- All permutations are stored as 0-based "one-line" lists perm[i] = image of i.
  Apply to columns as: G[:, perm].
"""

import json
import itertools
from collections import deque
import numpy as np


# ================================================================
#  GF(2) linear algebra utilities
# ================================================================

def gf2_rank(A: np.ndarray) -> int:
    """Rank over GF(2) for a uint8 matrix."""
    A = (A.copy().astype(np.uint8)) & 1
    m, n = A.shape
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if A[i, c] == 1:
                piv = i
                break
        if piv is None:
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        for i in range(m):
            if i != r and A[i, c] == 1:
                A[i, :] ^= A[r, :]
        r += 1
        if r == m:
            break
    return r


def apply_col_perm(G: np.ndarray, perm0: list[int]) -> np.ndarray:
    """
    Apply coordinate permutation i -> perm0[i] to the *codeword positions*.

    Convention:
      - perm0[i] = image of coordinate i
      - action on a codeword support set S: S -> {perm0[i] : i in S}

    For generator matrix columns, that means:
      new[:, perm0[i]] = old[:, i]
    hence:
      new[:, j] = old[:, inv_perm[j]]
    """
    inv = perm_inv(perm0)
    return G[:, inv]


def is_code_automorphism(G: np.ndarray, perm0: list[int], verbose: bool = False) -> bool:
    """
    Check RowSpan(G) == RowSpan(GP) using rank([G; GP]) == rank(G).
    """
    Gp = apply_col_perm(G, perm0)
    rG = gf2_rank(G)
    rStack = gf2_rank(np.vstack([G, Gp]))
    ok = (rStack == rG)
    if verbose:
        print(f"[check] rank(G)={rG}, rank([G;Gperm])={rStack} -> {'OK' if ok else 'FAIL'}")
    return ok


# ================================================================
#  Permutation utilities
# ================================================================

def cycles_to_perm(n: int, cycles_1based: list[tuple[int, ...]]) -> list[int]:
    """1-based disjoint cycles -> 0-based one-line permutation perm[i]=image(i)."""
    img = {i: i for i in range(1, n + 1)}
    for cyc in cycles_1based:
        k = len(cyc)
        for t in range(k):
            img[cyc[t]] = cyc[(t + 1) % k]
    return [img[i + 1] - 1 for i in range(n)]


def perm_inv(p: list[int]) -> list[int]:
    q = [0] * len(p)
    for i, j in enumerate(p):
        q[j] = i
    return q


def perm_comp(p: list[int], q: list[int]) -> list[int]:
    """Composition p∘q (apply q then p): (p∘q)[i] = p[q[i]]."""
    return [p[q[i]] for i in range(len(p))]


def conj_s_inv_g_s(g, s):
    """s^{-1} g s"""
    return perm_comp(perm_inv(s), perm_comp(g, s))

def conj_s_g_s_inv(g, s):
    """s g s^{-1}"""
    return perm_comp(s, perm_comp(g, perm_inv(s)))

def perm_conj(g: list[int], s: list[int]) -> list[int]:
    """Conjugate g by s: s^{-1} g s."""
    s_inv = perm_inv(s)
    return perm_comp(s_inv, perm_comp(g, s))


def perm_apply_to_set(p: list[int], S: frozenset[int]) -> frozenset[int]:
    return frozenset(p[i] for i in S)


def perm_order(perm: list[int]) -> int:
    """Permutation order = lcm of disjoint cycle lengths."""
    n = len(perm)
    seen = [False] * n
    from math import gcd
    def lcm(a, b): return a // gcd(a, b) * b
    ordv = 1
    for i in range(n):
        if seen[i]:
            continue
        j = i
        cnt = 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            cnt += 1
        if cnt > 1:
            ordv = lcm(ordv, cnt)
    return ordv


# ================================================================
#  Your eGolay [24,12,8] generator matrix (12x24)
# ================================================================

def egolay_G24_systematic_yours() -> np.ndarray:
    # rows = [
    #    '110001110101000000000001',
    #    '011000111010100000000001',
    #    '001100011101010000000001',
    #    '000110001110101000000001',
    #    '000011000111010100000001',
    #    '000001100011101010000001',
    #    '000000110001110101000001',
    #    '000000011000111010100001',
    #    '000000001100011101010001',
    #    '000000000110001110101001',
    #    '000000000011000111010101',
    #    '000000000001100011101011',
    # ]
    rows = [
        "100000000000110001110101",
        "010000000000011000111011",
        "001000000000111101101000",
        "000100000000011110110100",
        "000010000000001111011010",
        "000001000000110110011001",
        "000000100000011011001101",
        "000000010000001101100111",
        "000000001000110111000110",
        "000000000100101010010111",
        "000000000010100100111110",
        "000000000001100011101011",
    ]
    G = np.array([[int(ch) for ch in r.strip()] for r in rows], dtype=np.uint8)
    assert G.shape == (12, 24)
    return G


# ================================================================
#  Enumerate codewords and octads (weight-8 supports)
# ================================================================

def all_codewords(G: np.ndarray) -> np.ndarray:
    """Enumerate all 2^k codewords for k=12. Returns array shape (4096, 24) uint8."""
    k, n = G.shape
    assert k == 12 and n == 24
    C = np.zeros((1 << k, n), dtype=np.uint8)
    for u in range(1 << k):
        v = np.zeros(n, dtype=np.uint8)
        x = u
        r = 0
        while r < k:
            if x & 1:
                v ^= G[r, :]
            x >>= 1
            r += 1
        C[u, :] = v
    return C


def octads_from_code(G: np.ndarray) -> set[frozenset[int]]:
    """Set of octads (supports of weight-8 codewords), indices 0..23."""
    C = all_codewords(G)
    octads: set[frozenset[int]] = set()
    for v in C:
        if int(v.sum()) == 8:
            S = frozenset(np.nonzero(v)[0].tolist())
            octads.add(S)
    return octads


def build_5_to_octad_map(octads: set[frozenset[int]]) -> dict[frozenset[int], frozenset[int]]:
    """
    For Steiner S(5,8,24): every 5-subset sits in a UNIQUE octad.
    Build lookup: key=5-subset -> value=octad.
    """
    mp: dict[frozenset[int], frozenset[int]] = {}
    for O in octads:
        pts = sorted(O)
        for comb in itertools.combinations(pts, 5):
            key = frozenset(comb)
            if key in mp:
                # In a true S(5,8,24) this should never happen (uniqueness).
                # But keep a defensive check.
                if mp[key] != O:
                    raise RuntimeError("Not a Steiner S(5,8,24): a 5-set belongs to two different octads.")
            else:
                mp[key] = O
    # sanity: number of 5-subsets = C(24,5) = 42504
    if len(mp) != 42504:
        raise RuntimeError(f"5->octad map size {len(mp)} != C(24,5)=42504. Something is inconsistent.")
    return mp


# ================================================================
#  ATLAS standard generators b11, b21 for M24 (ATLAS labeling)
# ================================================================

def m24_atlas_generators() -> list[tuple[str, dict, list[int]]]:
    # b11 := (1,4)(2,7)(3,17)(5,13)(6,9)(8,15)(10,19)(11,18)(12,21)(14,16)(20,24)(22,23)
    b11_cycles = [
        (1, 4), (2, 7), (3, 17), (5, 13), (6, 9), (8, 15),
        (10, 19), (11, 18), (12, 21), (14, 16), (20, 24), (22, 23),
    ]
    b11 = cycles_to_perm(24, b11_cycles)

    # b21 := (1,4,6)(2,21,14)(3,9,15)(5,18,10)(13,17,16)(19,24,23)
    b21_cycles = [
        (1, 4, 6),
        (2, 21, 14),
        (3, 9, 15),
        (5, 18, 10),
        (13, 17, 16),
        (19, 24, 23),
    ]
    b21 = cycles_to_perm(24, b21_cycles)

    gens = [
        ("atlas_standard", {"name": "b11", "order": perm_order(b11)}, b11),
        ("atlas_standard", {"name": "b21", "order": perm_order(b21)}, b21),
    ]
    return gens


def orbit_octads_from_seed(seed_octad_1based: list[int], gens0: list[list[int]]) -> set[frozenset[int]]:
    """
    Build the set of octads in ATLAS labeling by orbit of a seed octad under the group generated by gens0.
    """
    seed = frozenset([x - 1 for x in seed_octad_1based])  # 0-based
    seen = set([seed])
    q = deque([seed])
    while q:
        O = q.popleft()
        for g in gens0:
            O2 = perm_apply_to_set(g, O)
            if O2 not in seen:
                seen.add(O2)
                q.append(O2)
    return seen


# ================================================================
#  Find relabeling s such that YOUR_octads = { s(A) : A in ATLAS_octads }
# ================================================================

def find_label_isomorphism(
    atlas_octads: set[frozenset[int]],
    yours_octads: set[frozenset[int]],
    max_seed_octads_to_try: int = 40,
) -> list[int]:
    """
    Find a permutation s (0-based list length 24) so that applying s to every ATLAS octad
    yields exactly your octads.

    Strategy:
    - Use the Steiner property via 5->octad maps for strong constraints.
    - Backtracking with pruning.
    - To reduce symmetry, we first try mapping ATLAS seed octad to some candidate octad in your system.
      We don't try all 759 by default; max_seed_octads_to_try controls how many candidates to attempt.

    Returns:
      s as list length 24 with s[a] = image of ATLAS point a in YOUR labeling.
    """
    atlas_5 = build_5_to_octad_map(atlas_octads)
    yours_5 = build_5_to_octad_map(yours_octads)

    atlas_octads_list = list(atlas_octads)
    # pick a "seed" atlas octad deterministically (first in sorted order for reproducibility)
    atlas_seed = min(atlas_octads_list, key=lambda O: tuple(sorted(O)))
    yours_octads_list = sorted(list(yours_octads), key=lambda O: tuple(sorted(O)))

    # helper: consistency check using 5-sets fully mapped
    def consistent_partial(f: dict[int, int], finv: dict[int, int]) -> bool:
        dom = sorted(f.keys())
        if len(dom) < 5:
            return True
        # check all 5-subsets of currently assigned atlas points
        for comb in itertools.combinations(dom, 5):
            A5 = frozenset(comb)
            A8 = atlas_5[A5]
            Y5 = frozenset(f[x] for x in A5)
            Y8 = yours_5.get(Y5)
            if Y8 is None:
                return False
            # already-mapped points inside A8 must land inside Y8
            for x in A8:
                if x in f and f[x] not in Y8:
                    return False
            # also, images must be injective
            if len(set(f.values())) != len(f):
                return False
        return True

    # propagation: if 5-set mapped, constrain the images of remaining points in its octad
    def allowed_images_for_unmapped_points(f: dict[int, int], finv: dict[int, int]) -> dict[int, set[int]]:
        """
        Build a domain restriction map for unmapped atlas points using current assignments.
        For each constraint induced by a mapped 5-set, restrict the 3 remaining points in that atlas octad
        to map into the remaining points of the target your octad.
        """
        dom_restrict: dict[int, set[int]] = {}
        dom = sorted(f.keys())
        if len(dom) < 5:
            return dom_restrict

        used_imgs = set(f.values())

        for comb in itertools.combinations(dom, 5):
            A5 = frozenset(comb)
            A8 = atlas_5[A5]
            Y5 = frozenset(f[x] for x in A5)
            Y8 = yours_5[Y5]
            # remaining available target points in that octad:
            # (excluding already used images, but keeping images of already-mapped points inside A8)
            already_inside = set()
            for x in A8:
                if x in f:
                    already_inside.add(f[x])
            if not already_inside.issubset(Y8):
                return {"__FAIL__": set()}  # signal failure
            avail = set(Y8) - (used_imgs - already_inside)

            # for each unmapped x in this A8, intersect its allowed set with 'avail'
            for x in A8:
                if x in f:
                    continue
                if x not in dom_restrict:
                    dom_restrict[x] = set(avail)
                else:
                    dom_restrict[x] &= avail
                if len(dom_restrict[x]) == 0:
                    return {"__FAIL__": set()}
        return dom_restrict

    # order of choosing next variable: prefer those with smallest inferred domain
    def choose_next_point(f: dict[int, int], dom_restrict: dict[int, set[int]]) -> int:
        unmapped = [i for i in range(24) if i not in f]
        # If we have inferred domain restrictions, pick the tightest
        best = None
        best_sz = 10**9
        for x in unmapped:
            if x in dom_restrict:
                sz = len(dom_restrict[x])
                if sz < best_sz:
                    best_sz = sz
                    best = x
        if best is not None:
            return best
        # else just pick smallest unmapped
        return unmapped[0]

    # try mapping the atlas_seed octad to a handful of your octads
    for idx_try, Yseed in enumerate(yours_octads_list[:max_seed_octads_to_try]):
        # Start with a partial mapping by mapping the entire seed octad as a set,
        # but we still need point-to-point mapping within it.
        # We'll backtrack; first fix one point of atlas_seed to one point of Yseed to break symmetry.
        atlas_seed_pts = sorted(atlas_seed)
        yseed_pts = sorted(Yseed)

        # Fix atlas_seed_pts[0] -> yseed_pts[0] to reduce branching.
        f0 = {atlas_seed_pts[0]: yseed_pts[0]}
        finv0 = {yseed_pts[0]: atlas_seed_pts[0]}

        # Also require that the whole atlas_seed octad maps into Yseed as a set:
        # i.e., images of its 8 points must be within Yseed.
        seed_allowed = {x: set(yseed_pts) for x in atlas_seed_pts}

        # recursive search
        def dfs(f: dict[int, int], finv: dict[int, int]) -> list[int] | None:
            # basic checks
            if not consistent_partial(f, finv):
                return None

            # infer restrictions from all mapped 5-sets
            dom_restrict = allowed_images_for_unmapped_points(f, finv)
            if "__FAIL__" in dom_restrict:
                return None

            # include seed-octad "must land in Yseed" restriction
            for x in atlas_seed_pts:
                if x in f:
                    if f[x] not in Yseed:
                        return None
                else:
                    dom_restrict[x] = dom_restrict.get(x, set(range(24))) & set(Yseed)
                    if len(dom_restrict[x]) == 0:
                        return None

            # done?
            if len(f) == 24:
                # Verify octad-set equality
                s = [None] * 24
                for a, y in f.items():
                    s[a] = y
                if any(v is None for v in s):
                    return None
                s = [int(v) for v in s]
                mapped = {perm_apply_to_set(s, O) for O in atlas_octads}
                if mapped == yours_octads:
                    return s
                return None

            # choose next point
            x = choose_next_point(f, dom_restrict)
            # candidate images
            if x in dom_restrict:
                cand = sorted(dom_restrict[x])
            else:
                cand = [y for y in range(24) if y not in finv]

            for y in cand:
                if y in finv:
                    continue
                # if x is in atlas_seed, force into Yseed (already enforced in dom_restrict, but keep safe)
                if x in atlas_seed and y not in Yseed:
                    continue
                f[x] = y
                finv[y] = x
                out = dfs(f, finv)
                if out is not None:
                    return out
                del finv[y]
                del f[x]
            return None

        s = dfs(dict(f0), dict(finv0))
        if s is not None:
            print(f"[info] Found label isomorphism using Yseed try #{idx_try+1}")
            return s

    raise RuntimeError("Failed to find a relabeling permutation s. Increase max_seed_octads_to_try.")


# ================================================================
#  JSON saving (atomic replace)
# ================================================================

def save_checkpoint_json(path: str, data: dict) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    import os
    os.replace(tmp, path)


# ================================================================
#  Main
# ================================================================

def main():
    # --- your code and its octads ---
    G = egolay_G24_systematic_yours()
    print("[info] Building octads from your embedded G (enumerate 2^12 codewords)...")
    yours_octads = octads_from_code(G)
    print(f"[info] #your_octads = {len(yours_octads)} (expected 759)")

    if len(yours_octads) != 759:
        raise RuntimeError("Your embedded G does not yield 759 octads; not a valid eGolay labeling?")

    # --- ATLAS generators and ATLAS octads ---
    gens = m24_atlas_generators()
    gens0 = [perm for (_, _, perm) in gens]

    # ATLAS seed octad (1-based): {1,2,3,4,5,11,17,24}
    atlas_seed_octad = [1, 2, 3, 4, 5, 11, 17, 24]

    print("[info] Building ATLAS octads by orbiting the ATLAS seed octad under <b11,b21>...")
    atlas_octads = orbit_octads_from_seed(atlas_seed_octad, gens0)
    print(f"[info] #atlas_octads = {len(atlas_octads)} (expected 759)")

    if len(atlas_octads) != 759:
        raise RuntimeError("ATLAS octad orbit did not produce 759 octads. Generator definitions wrong?")

    # --- find relabeling s ---
    print("[info] Finding relabeling permutation s to match ATLAS octads to your octads...")
    s = find_label_isomorphism(atlas_octads, yours_octads, max_seed_octads_to_try=60)
    print("[info] Found s (0-based):", s)

    # --- conjugate generators into YOUR labeling ---
    gens_conj = []
    for typ, meta, g in gens:
        # g2 = perm_conj(g, s)  # s^{-1} g s
        g2 = conj_s_g_s_inv(g, s)
        
        # # ========= debug section ============
        # # pick any atlas octad O, map it to your labels, then apply g2
        # O = next(iter(atlas_octads))
        # lhs = perm_apply_to_set(g2, perm_apply_to_set(s, O))

        # # apply g in atlas first, then map by s
        # rhs = perm_apply_to_set(s, perm_apply_to_set(g, O))

        # print("[debug] octad equiv:", lhs == rhs)
        # # ========= debug section ============

        meta2 = dict(meta)
        meta2["name"] = meta2["name"] + "_conj_to_yours"
        meta2["order"] = perm_order(g2)
        gens_conj.append(("m24_generator", meta2, g2))

    # --- automorphism check section ---
    print("[info] Checking conjugated generators against your embedded G...")
    ok_all = True
    for typ, meta, perm in gens_conj:
        ok = is_code_automorphism(G, perm, verbose=True)
        ok_all &= ok
        print(f"        - {meta['name']}: order={meta['order']} -> {'OK' if ok else 'FAIL'}")
    if not ok_all:
        raise RuntimeError("Conjugated generators still failed automorphism check. Something inconsistent.")

    # --- store JSON (style similar to generator_set_eBCH.py) ---
    H_perms = [perm for (_, _, perm) in gens_conj]
    data = {
        "code": "eGolay24",
        "n": int(G.shape[1]),
        "k": int(G.shape[0]),
        "|H|": len(H_perms),
        "H_perms": H_perms,
        "generator_records": [{"type": typ, "meta": meta} for (typ, meta, _) in gens_conj],
        "relabeling_s_atlas_to_yours": s,
        "notes": [
            "All perms are 0-based one-line lists: perm[i] = image(i).",
            "Apply to columns as G[:, perm].",
            "These generators are ATLAS (b11,b21) conjugated by s to match YOUR coordinate labeling.",
        ],
        "G_rows_bits": ["".join(str(int(b)) for b in row.tolist()) for row in G],
    }

    out_path = "generator_set_eGolay24.json"
    save_checkpoint_json(out_path, data)
    print(f"[info] Saved -> {out_path}")


if __name__ == "__main__":
    main()
