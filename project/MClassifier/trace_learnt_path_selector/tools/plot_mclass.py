#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt


def read_eval(path):
    vals = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        key, val = line.split(",", 1)
        vals[key] = float(val)
    return vals


def main():
    parser = argparse.ArgumentParser(description="Plot M-classifier eval summaries.")
    parser.add_argument("--eval", nargs="+", required=True, help="mclass_eval CSV summaries.")
    parser.add_argument("--labels", nargs="*", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    labels = args.labels or [Path(p).stem for p in args.eval]
    vals = [read_eval(p) for p in args.eval]

    x = range(len(vals))
    plt.figure(figsize=(7, 4))
    plt.semilogy(x, [v["teacher_bler"] for v in vals], marker="o", label="teacher")
    plt.semilogy(x, [v["top1_bler"] for v in vals], marker="s", label="top1")
    plt.semilogy(x, [v["topk_bler"] for v in vals], marker="^", label="topk")
    plt.xticks(list(x), labels, rotation=30, ha="right")
    plt.ylabel("BLER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
