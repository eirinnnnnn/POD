#!/usr/bin/env python3
"""
convert_bsgs_to_perm.py  (BSGS sampling; no group enumeration)

Input JSON (from patched schreier_sims.py) must contain:
  - n: int
  - base: list[int]
  - levels: list of dicts, each with:
      * base_point: int
      * orbit: list[int]
      * trans: list[list[int]]   # aligned with orbit: trans[j](base_point) = orbit[j]

Output (AED loader format):
  Header: "<num_perms> <n>"
  Rows:   space-separated permutation of 0..n-1, one per line

Conventions:
  - A permutation p is a list[int] length n with action: i -> p[i]
  - Composition compose(a,b) = a ∘ b (apply b then a): (a∘b)[i] = a[b[i]]
"""

import argparse
import json
import random
from typing import List, Tuple, Iterable


Perm = Tuple[int, ...]


def compose(a: Perm, b: Perm) -> Perm:
    """a ∘ b (apply b then a)."""
    return tuple(a[b[i]] for i in range(len(a)))


def identity(n: int) -> Perm:
    return tuple(range(n))


def sample_from_levels(levels: List[dict], n: int, rng: random.Random) -> Perm:
    """
    Sample a random-ish element from the stabilizer chain:
      g = t_1,x1 ∘ t_2,x2 ∘ ... ∘ t_r,xr
    where x_i is chosen uniformly from orbit_i and t_i,xi is the stored transversal.
    """
    g = identity(n)
    for lvl in levels:
        orbit: List[int] = lvl["orbit"]
        trans_list: List[List[int]] = lvl["trans"]  # aligned with orbit
        j = rng.randrange(len(orbit))
        t = tuple(trans_list[j])
        g = compose(t, g)  # left-multiply
    return g


def write_perms(path: str, perms: Iterable[Perm], n: int) -> None:
    perms = list(perms)
    with open(path, "w", encoding="ascii") as f:
        f.write(f"{len(perms)} {n}\n")
        for p in perms:
            f.write(" ".join(str(x) for x in p))
            f.write("\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Schreier–Sims BSGS JSON to permutation rows for AED (sampling, no enumeration).")
    ap.add_argument("json_in", help="Input *_schreier_sims.json (must include levels[*].orbit and levels[*].trans)")
    ap.add_argument("out_path", help="Output path for permutation rows")
    ap.add_argument("--L", type=int, default=16, help="Number of permutations to output (default: 16)")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed (default: None)")
    ap.add_argument("--include_identity", action="store_true", help="Include identity as the first permutation")
    ap.add_argument("--max_tries", type=int, default=100000, help="Max attempts to collect L unique perms")

    args = ap.parse_args()

    with open(args.json_in, "r", encoding="utf-8") as f:
        data = json.load(f)

    if "levels" not in data:
        raise KeyError("Missing key 'levels'. Did you patch schreier_sims.py to store transversals?")
    n = int(data["n"])
    levels = data["levels"]
    print(levels[0].keys())

    # Validate that trans is present and aligned
    for i, lvl in enumerate(levels, start=1):
        if "orbit" not in lvl or "trans" not in lvl:
            raise KeyError(f"Level {i} missing 'orbit' or 'trans'. Re-run patched schreier_sims.py.")
        if len(lvl["orbit"]) != len(lvl["trans"]):
            raise ValueError(f"Level {i}: len(orbit) != len(trans). Expected aligned storage.")

    rng = random.Random(args.seed)

    perms: List[Perm] = []
    seen = set()

    if args.include_identity:
        e = identity(n)
        perms.append(e)
        seen.add(e)

    tries = 0
    while len(perms) < args.L and tries < args.max_tries:
        g = sample_from_levels(levels, n, rng)
        if g not in seen:
            seen.add(g)
            perms.append(g)
        tries += 1

    if len(perms) < args.L:
        raise RuntimeError(f"Only collected {len(perms)} unique perms after {tries} tries. Increase --max_tries or reduce --L.")

    write_perms(args.out_path, perms, n)
    print(f"Wrote {len(perms)} permutations (n={n}) to {args.out_path}")


if __name__ == "__main__":
    main()
