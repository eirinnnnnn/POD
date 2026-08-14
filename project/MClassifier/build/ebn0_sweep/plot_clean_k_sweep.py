#!/usr/bin/env python3
"""BLER vs k, one subplot per checkpoint t, clean-automorphism-file rerun.
Compares static-pm-rank (model-free) vs trace-learned (retrained on clean
data) against the full-PED (M=64) reference line."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "clean_k_sweep_results.txt"

STYLE = {
    "static_pm_rank": dict(color="#e34948", marker="^", linestyle="-", label="static pm-rank"),
    "trace_learned": dict(color="#008300", marker="s", linestyle="--", label="trace-learned (clean retrain)"),
}


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    ts = sorted({int(r["t"]) for r in rows})
    series = defaultdict(lambda: defaultdict(list))
    full_ped = {}
    for r in rows:
        t = int(r["t"])
        if r["method"] == "full_ped":
            full_ped[t] = float(r["bler"])
        else:
            series[t][r["method"]].append((int(r["k"]), float(r["bler"])))

    fig, axes = plt.subplots(1, len(ts), figsize=(6.2 * len(ts), 5.2), sharey=True)
    for ax, t in zip(axes, ts):
        ax.axhline(full_ped[t], color="#2a78d6", linestyle="-", linewidth=1.6, label="full PED (M=64)")
        for method, style in STYLE.items():
            pts = sorted(series[t][method])
            ax.plot([k for k, _ in pts], [b for _, b in pts], **style)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("k")
        ax.set_title(f"t = {t}")
        ax.grid(True, which="both", alpha=0.3)

    axes[0].set_ylabel("BLER")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("m7t10_m64, SNR=2.0dB: static-pm-rank vs trace-learned vs full-PED (clean automorphism file)")
    fig.tight_layout()

    out = HERE / "clean_k_sweep_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
