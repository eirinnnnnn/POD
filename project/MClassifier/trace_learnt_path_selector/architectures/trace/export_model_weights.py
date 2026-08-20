#!/usr/bin/env python3
"""Export a trained trace-ranker checkpoint (.pt) to a plain-text weight
file a streaming C++ simulator can load without linking torch. The MLP is
just Linear->ReLU->Linear->ReLU->Linear, small enough to reimplement
directly (input_dim<=~100, hidden=64) -- this lets the BLER simulator
score branches inline during a single decode pass instead of writing
per-sample trace features to disk for a separate Python prediction step."""
import argparse
from pathlib import Path

import numpy as np
import torch


def write_vec(f, name, arr):
    arr = np.asarray(arr, dtype=np.float64).reshape(-1)
    f.write(f"{name} {len(arr)}\n")
    f.write(" ".join(repr(float(x)) for x in arr) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ckpt = torch.load(args.model, map_location="cpu", weights_only=False)
    sd = ckpt["state_dict"]

    out = Path(args.out)
    with out.open("w", encoding="utf-8") as f:
        f.write(f"input_dim {int(ckpt['input_dim'])}\n")
        f.write(f"hidden {int(ckpt['hidden'])}\n")
        f.write(f"M {int(ckpt['M'])}\n")
        f.write(f"branch_onehot {1 if ckpt['branch_onehot'] else 0}\n")
        f.write(f"direction_largest {1 if ckpt['direction'] == 'largest' else 0}\n")
        write_vec(f, "mean", ckpt["mean"])
        write_vec(f, "std", ckpt["std"])
        write_vec(f, "W1", sd["0.weight"].numpy())   # [hidden, input_dim], row-major
        write_vec(f, "b1", sd["0.bias"].numpy())
        write_vec(f, "W2", sd["2.weight"].numpy())   # [hidden, hidden]
        write_vec(f, "b2", sd["2.bias"].numpy())
        write_vec(f, "W3", sd["4.weight"].numpy())   # [1, hidden]
        write_vec(f, "b3", sd["4.bias"].numpy())
        f.write("feature_cols " + str(len(ckpt["feature_cols"])) + "\n")
        f.write(" ".join(ckpt["feature_cols"]) + "\n")

    print(f"wrote {out}")


if __name__ == "__main__":
    main()
