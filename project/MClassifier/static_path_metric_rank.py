#!/usr/bin/env python3
"""Model-free branch ranking: at checkpoint t, rank the M branches by their
own raw partial path metric (t{t}_pm_min) and keep the top-k -- no training,
no learned scorer. This is the "just keep the currently-best k parallel
paths" baseline against which the trained trace-checkpoint classifiers
should be measured."""
import argparse
import csv
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
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise SystemExit("empty trace csv")
    if col_name not in rows[0]:
        raise SystemExit(f"column {col_name} not found in {args.trace}")

    n_samples = max(int(r["sample"]) for r in rows) + 1
    M = max(int(r["branch"]) for r in rows) + 1

    vals = np.full((n_samples, M), np.inf, dtype=np.float64)
    for r in rows:
        vals[int(r["sample"]), int(r["branch"])] = float(r[col_name])

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
