#!/usr/bin/env python3
"""BLER-vs-Eb/N0 plot for the Golay24 M=16/list_size=2/t=8 checkpoint-ranked
static-pm-rank sweep: full-PED (M=16) vs static-pm-rank at m=1,2,3,4,8,
over whatever SNR points are currently done."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "sweep" / "golay_m16_t8_L2_results.txt"

M_LIST = [1, 2, 3, 4, 8]


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


COLORS = dict(zip(M_LIST, _ramp("Reds", len(M_LIST))))


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    full_ped = defaultdict(list)
    static = defaultdict(list)
    for r in rows:
        pt = (float(r["snr"]), float(r["bler"]))
        if r["method"] == "full_ped":
            full_ped[16].append(pt)
        elif r["method"] == "static_pm_rank":
            static[int(r["m"])].append(pt)

    fig, ax = plt.subplots(figsize=(9, 7))

    pts = sorted(full_ped[16])
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#2a78d6",
                 linestyle="-", marker="o", markersize=6, linewidth=2,
                 label="full PED (M=16)", zorder=3)

    for m in M_LIST:
        pts = sorted(static.get(m, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS[m],
                 linestyle="--", marker="^", markersize=5, linewidth=1.4,
                 label=f"static pm-rank (t=8, m={m})", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Golay(24,12) AED (M=16, list_size=2, checkpoint t=8): full-PED vs static-pm-rank")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "golay_m16_t8_L2_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
