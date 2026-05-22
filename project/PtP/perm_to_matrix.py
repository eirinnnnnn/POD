#!/usr/bin/env python3
import argparse
import ast


def write_perm_matrix(perm, output_name):
    n = len(perm)

    # Basic validation
    if sorted(perm) != list(range(n)):
        raise ValueError(
            f"Input list must be a permutation of 0,...,{n-1}."
        )

    if not output_name.endswith(".matrix"):
        output_name += ".matrix"

    with open(output_name, "w") as f:
        f.write(f"{n} {n}\n")

        for i in range(n):
            row = ["0"] * n
            row[perm[i]] = "1"
            f.write(" ".join(row) + "\n")

    print(f"Wrote permutation matrix to {output_name}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a permutation list into .matrix format."
    )

    parser.add_argument(
        "--perm",
        type=str,
        required=True,
        help='Permutation list, e.g. "[2,0,1,3]"',
    )

    parser.add_argument(
        "--name",
        type=str,
        required=True,
        help='Output filename, e.g. "pi.matrix" or "pi"',
    )

    args = parser.parse_args()

    perm = ast.literal_eval(args.perm)

    if not isinstance(perm, list):
        raise ValueError("--perm must be a Python-style list.")

    write_perm_matrix(perm, args.name)


if __name__ == "__main__":
    main()