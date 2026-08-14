#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_eval(path):
    vals = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            a, b = line.strip().split(",")
            vals[a] = float(b)
    return vals


def main():
    parser = argparse.ArgumentParser(description="Evaluate best-basin hit rate and BLER.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--pred", action="append", nargs=2, metavar=("NAME", "CSV"), required=True)
    args = parser.parse_args()

    ks = [1, 2, 4, 8, 16, 32, 64]
    rows = list(csv.DictReader(open(args.dataset, encoding="utf-8")))
    M = len([c for c in rows[0] if c.startswith("metric_") and c[7:].isdigit()])
    metrics = np.array([[float(r[f"metric_{i}"]) for i in range(M)] for r in rows])
    best = metrics.min(axis=1, keepdims=True)
    basin = np.abs(metrics - best) <= args.delta
    labels = np.array([int(r["metric_argmin"]) for r in rows])

    summary = {}
    print("model,k,bler,basin_hit,argmin_hit")
    for name, pred_path in args.pred:
        pred_path = Path(pred_path)
        for k in ks:
            out = subprocess.check_output(
                ["./build/mclass_eval", "--dataset", args.dataset, "--pred", str(pred_path), "--topk", str(k)],
                text=True,
            )
            Path(f"build/eval_{args.prefix}_{name}_top{k}.csv").write_text(out, encoding="utf-8")

        pred_rows = list(csv.DictReader(open(pred_path, encoding="utf-8")))
        pred_order = np.array([[int(r[f"idx{i}"]) for i in range(M)] for r in pred_rows])
        blers = []
        hits = []
        arg_hits = []
        for k in ks:
            vals = read_eval(f"build/eval_{args.prefix}_{name}_top{k}.csv")
            top = pred_order[:, :k]
            hit = basin[np.arange(len(rows))[:, None], top].any(axis=1).mean()
            arg_hit = (top == labels[:, None]).any(axis=1).mean()
            bler = vals["topk_bler"]
            blers.append(bler)
            hits.append(hit)
            arg_hits.append(arg_hit)
            print(f"{name},{k},{bler:.6f},{hit:.6f},{arg_hit:.6f}")
        summary[name] = (blers, hits, arg_hits)

    plt.figure(figsize=(7, 4))
    for name, (_, hits, _) in summary.items():
        plt.plot(ks, hits, marker="o", label=name)
    plt.xscale("log", base=2)
    plt.xticks(ks, [str(k) for k in ks])
    plt.xlabel("predicted top-k branches")
    plt.ylabel("best-basin hit rate")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    hit_plot = f"build/{args.prefix}_basin_hit_compare.png"
    plt.savefig(hit_plot, dpi=200)
    print(f"wrote {hit_plot}")

    teacher = []
    first_name = args.pred[0][0]
    for k in ks:
        teacher.append(read_eval(f"build/eval_{args.prefix}_{first_name}_top{k}.csv")["teacher_bler"])

    plt.figure(figsize=(7, 4))
    plt.semilogy(ks, teacher, marker="o", label="full PED teacher")
    for name, (blers, _, _) in summary.items():
        plt.semilogy(ks, blers, marker="o", label=name)
    plt.xscale("log", base=2)
    plt.xticks(ks, [str(k) for k in ks])
    plt.xlabel("decoded AED branches k out of M=64")
    plt.ylabel("BLER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    bler_plot = f"build/{args.prefix}_basin_bler_compare.png"
    plt.savefig(bler_plot, dpi=200)
    print(f"wrote {bler_plot}")


if __name__ == "__main__":
    main()
