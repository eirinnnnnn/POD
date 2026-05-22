import argparse
import json
import os


# ================================================================
#  JSON utility
# ================================================================

def save_checkpoint_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, path)
    print(f"[JSON saved → {path}]")


# ================================================================
#  LTA(m, 2) permutations on coordinates {0, ..., 2^m - 1}
# ================================================================


def int_to_bits(x, m):
    """
    LSB-first binary vector:
        x = sum_{i=0}^{m-1} bits[i] 2^i.
    """
    return [(x >> i) & 1 for i in range(m)]


def bits_to_int(bits):
    """Inverse of int_to_bits under the same LSB-first convention."""
    y = 0
    for i, b in enumerate(bits):
        if b & 1:
            y |= (1 << i)
    return y


def translation_perm(m, i):
    """
    Translation by the i-th standard basis vector:
        x_i -> x_i + 1.

    As a permutation list perm, coordinate j is sent to perm[j].
    """
    n = 1 << m
    return [x ^ (1 << i) for x in range(n)]


def lower_shear_perm(m, i, j):
    """
    Elementary lower-triangular shear over GF(2):
        x_i -> x_i + x_j,       for 0 <= j < i < m,
        x_l -> x_l              for l != i.

    With LSB-first coordinates, this generates the unit lower-triangular
    linear part together with all such pairs (i, j).

    As a permutation list perm, coordinate x is sent to perm[x].
    """
    if not (0 <= j < i < m):
        raise ValueError(f"Require 0 <= j < i < m, got i={i}, j={j}, m={m}")

    n = 1 << m
    perm = []
    for x in range(n):
        bits = int_to_bits(x, m)
        bits[i] ^= bits[j]
        perm.append(bits_to_int(bits))
    return perm


def is_permutation(perm):
    return sorted(perm) == list(range(len(perm)))


def compose_perm(p, q):
    """
    Composition p ∘ q under image-list convention:
        (p ∘ q)[x] = p[q[x]].
    """
    if len(p) != len(q):
        raise ValueError("Cannot compose permutations of different lengths")
    return [p[q[x]] for x in range(len(p))]


def generate_LTA_generators(m, include_identity=True, verify=True):
    """
    Generate a standard generator set for LTA(m, 2):

        x -> A x + b,

    where A is an m x m unit lower-triangular matrix over GF(2), and
    b in GF(2)^m.

    Generators:
      - identity, optional, for compatibility with existing JSON files
      - m translations x_i -> x_i + 1
      - m(m-1)/2 elementary lower shears x_i -> x_i + x_j, j < i

    The coordinate index convention is LSB-first:
        index(x_0, ..., x_{m-1}) = sum_i x_i 2^i.
    """
    if m < 1:
        raise ValueError("m must be >= 1")

    n = 1 << m
    gens = []

    if include_identity:
        gens.append(("identity", None, list(range(n))))

    for i in range(m):
        gens.append(("translation", {"bit": i}, translation_perm(m, i)))

    for i in range(1, m):
        for j in range(i):
            gens.append(("lower_shear", {"target_bit": i, "source_bit": j}, lower_shear_perm(m, i, j)))

    if verify:
        failures = []
        for typ, meta, perm in gens:
            if not is_permutation(perm):
                failures.append((typ, meta))
        if failures:
            raise RuntimeError(f"Non-permutation generators detected: {failures}")

    return gens


def expected_group_order_LTA(m):
    """
    |LTA(m,2)| = 2^m * 2^{m(m-1)/2}, since the translation part has
    2^m choices and the unit lower-triangular part has one free GF(2)
    entry below each diagonal position.
    """
    return 1 << (m + (m * (m - 1)) // 2)


# ================================================================
#  MAIN
# ================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Generate a Schreier-Sims-compatible generator set for LTA(m,2)."
    )
    parser.add_argument("--m", type=int, default=5, help="number of binary variables; n = 2^m")
    parser.add_argument(
        "--k",
        type=int,
        default=None,
        help="optional code dimension field kept only for compatibility with the existing JSON schema",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="output JSON path; default: generator_set_LTA_m{m}.json",
    )
    parser.add_argument(
        "--no-identity",
        action="store_true",
        help="do not include identity in H_perms; default includes it to match generator_set_eBCH.py style",
    )
    args = parser.parse_args()

    m = args.m
    n = 1 << m
    k = args.k if args.k is not None else n

    print(f"[info] Building LTA(m={m}, 2) generator set with n={n} coordinates...")
    gens = generate_LTA_generators(m, include_identity=(not args.no_identity), verify=True)
    H_perms = [perm for (_, _, perm) in gens]

    data = {
        "n": n,
        "k": k,
        "|H|": len(H_perms),
        "H_perms": H_perms,
        "generator_records": [
            {"type": typ, "meta": meta}
            for (typ, meta, _) in gens
        ],
        "group": "LTA(m,2)",
        "m": m,
        "coordinate_order": "LSB-first: index = sum_i x_i 2^i",
        "permutation_convention": "image list: perm[x] = transformed coordinate of x",
        "expected_group_order": expected_group_order_LTA(m),
    }

    out_path = args.output or f"generator_set_LTA_m{m}.json"
    save_checkpoint_json(out_path, data)

    print(f"[info] Saved LTA generator set -> {out_path}")
    print(f"[info] n={n}, k={k}, |H|={len(H_perms)}")
    print(f"[info] expected |LTA(m,2)|={expected_group_order_LTA(m)}")


if __name__ == "__main__":
    main()
