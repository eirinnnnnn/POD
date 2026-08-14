#!/usr/bin/env python3
"""Model-free branch ranking: at checkpoint t, rank the M branches by their
own raw partial path metric (t{t}_pm_min) and keep the top-k -- no training,
no learned scorer. This is the "just keep the currently-best k parallel
paths" baseline against which the trained trace-checkpoint classifiers
should be measured."""
import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", required=True, help="trace-features CSV (mclass_trace_dataset output)")
    parser.add_argument("--step", type=int, required=True, help="checkpoint t, e.g. 8/16/32/64/128")
    parser.add_argument("--out", required=True)
    parser.add_argument("--field", default="pm_min", choices=["pm_min", "pm_gap", "pm_mean", "pm_max"])
    args = parser.parse_args()

    col_name = f"t{args.step}_{args.field}"
    with open(args.trace, newline="", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    if col_name not in header:
        raise SystemExit(f"column {col_name} not found in {args.trace}")
    col = {name: i for i, name in enumerate(header)}

    # Stream only the columns we need straight into an array -- csv.DictReader
    # materializing every row as a dict of strings first is ruinous on
    # multi-GB trace files (many GB of Python object overhead). float64 to
    # match the original float()-parsed precision exactly (avoids flipping
    # stable-sort tie-breaks between near-equal metrics).
    usecols = (col["sample"], col["branch"], col[col_name])
    data = np.loadtxt(args.trace, delimiter=",", skiprows=1, dtype=np.float64, usecols=usecols)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    samples = data[:, 0].astype(np.int64)
    branches = data[:, 1].astype(np.int64)
    metric = data[:, 2]

    n_samples = int(samples.max()) + 1
    M = int(branches.max()) + 1

    vals = np.full((n_samples, M), np.inf, dtype=np.float64)
    vals[samples, branches] = metric

    # Lower metric = better throughout this codebase (metric_argmin convention).
    ranked = np.argsort(vals, axis=1, kind="stable")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("sample")
        for i in range(M):
            f.write(f",idx{i}")
        f.write("\n")
        for s in range(n_samples):
            f.write(str(s))
            for idx in ranked[s]:
                f.write(f",{int(idx)}")
            f.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
