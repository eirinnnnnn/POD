#!/usr/bin/env python3
"""BLER-vs-Eb/N0 plot: full-PED (M=16) vs static-pm-rank vs trace-learned,
checkpoint t=8, list_size=2, m in {1,2,3,4,8}, over whatever SNR points
are currently done."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "sweep" / "golay_trace_vs_static_results.txt"

M_LIST = [1, 2, 3, 4, 8]


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


STATIC_COLORS = dict(zip(M_LIST, _ramp("Reds", len(M_LIST))))
TRACE_COLORS = dict(zip(M_LIST, _ramp("Greens", len(M_LIST))))


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    full_ped = []
    static = defaultdict(list)
    trace = defaultdict(list)
    for r in rows:
        pt = (float(r["snr"]), float(r["bler"]))
        if r["method"] == "full_ped":
            full_ped.append(pt)
        elif r["method"] == "static_pm_rank":
            static[int(r["m"])].append(pt)
        elif r["method"] == "trace_learned":
            trace[int(r["m"])].append(pt)

    fig, ax = plt.subplots(figsize=(10, 7.5))

    pts = sorted(full_ped)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#2a78d6",
                 linestyle="-", marker="o", markersize=6, linewidth=2,
                 label="full PED (M=16)", zorder=3)

    for m in M_LIST:
        pts = sorted(static.get(m, []))
        if pts:
            ax.plot([s for s, _ in pts], [b for _, b in pts], color=STATIC_COLORS[m],
                     linestyle="--", marker="^", markersize=5, linewidth=1.4,
                     label=f"static pm-rank (m={m})", zorder=2)
        pts = sorted(trace.get(m, []))
        if pts:
            ax.plot([s for s, _ in pts], [b for _, b in pts], color=TRACE_COLORS[m],
                     linestyle=":", marker="s", markersize=5, linewidth=1.4,
                     label=f"trace-learned (m={m})", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Golay(24,12) AED (M=16, list_size=2, t=8): full-PED vs static-pm-rank vs trace-learned")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "golay_trace_vs_static_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
