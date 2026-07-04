#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Compare top-k BLER and best-basin hit rate.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--ks", default="1,2,4,8,16,32,64")
    parser.add_argument("--pred", action="append", nargs=2, metavar=("NAME", "CSV"), required=True)
    args = parser.parse_args()

    ks = [int(x) for x in args.ks.split(",") if x]
    rows = list(csv.DictReader(open(args.dataset, encoding="utf-8")))
    M = len([c for c in rows[0] if c.startswith("metric_") and c[7:].isdigit()])
    metrics = np.array([[float(r[f"metric_{i}"]) for i in range(M)] for r in rows])
    basin = np.abs(metrics - metrics.min(axis=1, keepdims=True)) <= args.delta
    labels = np.array([int(r["metric_argmin"]) for r in rows])

    summary = {}
    teacher = None
    print("model,k,bler,basin_hit,argmin_hit")
    for name, pred_path in args.pred:
        pred_rows = list(csv.DictReader(open(pred_path, encoding="utf-8")))
        pred_order = np.array([[int(r[f"idx{i}"]) for i in range(M)] for r in pred_rows])
        blers = []
        hits = []
        ahits = []
        local_teacher = []
        for k in ks:
            out = subprocess.check_output(
                ["./build/mclass_eval", "--dataset", args.dataset, "--pred", pred_path, "--topk", str(k)],
                text=True,
            )
            vals = {}
            for line in out.strip().splitlines():
                a, b = line.split(",")
                vals[a] = float(b)
            top = pred_order[:, :k]
            hit = basin[np.arange(len(rows))[:, None], top].any(axis=1).mean()
            ahit = (top == labels[:, None]).any(axis=1).mean()
            bler = vals["topk_bler"]
            local_teacher.append(vals["teacher_bler"])
            blers.append(bler)
            hits.append(hit)
            ahits.append(ahit)
            print(f"{name},{k},{bler:.6f},{hit:.6f},{ahit:.6f}")
        if teacher is None:
            teacher = local_teacher
        summary[name] = (blers, hits, ahits)

    plt.figure(figsize=(7, 4))
    plt.semilogy(ks, teacher, marker="o", label="full PED teacher")
    for name, (blers, _, _) in summary.items():
        plt.semilogy(ks, blers, marker="o", label=name)
    plt.xscale("log", base=2)
    plt.xticks(ks, [str(k) for k in ks])
    plt.xlabel(f"decoded AED branches k out of M={M}")
    plt.ylabel("BLER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    bler_plot = f"build/{args.prefix}_bler.png"
    plt.savefig(bler_plot, dpi=200)
    print(f"wrote {bler_plot}")

    plt.figure(figsize=(7, 4))
    for name, (_, hits, _) in summary.items():
        plt.plot(ks, hits, marker="o", label=name)
    plt.xscale("log", base=2)
    plt.xticks(ks, [str(k) for k in ks])
    plt.xlabel(f"predicted top-k branches out of M={M}")
    plt.ylabel("best-basin hit rate")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    hit_plot = f"build/{args.prefix}_basin_hit.png"
    plt.savefig(hit_plot, dpi=200)
    print(f"wrote {hit_plot}")


if __name__ == "__main__":
    main()
