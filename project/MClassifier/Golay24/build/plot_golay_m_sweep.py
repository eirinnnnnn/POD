#!/usr/bin/env python3
"""BLER-vs-Eb/N0 plot for the Golay24 static-pm-rank sweep: full-PED (M=64)
vs static-pm-rank at m=1,4,8,16,32 vs pure-SCL(L=8), over whatever SNR
points are currently done."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "sweep" / "golay_m_sweep_results.txt"

M_LIST = [1, 4, 8, 16, 32]


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


COLORS = dict(zip(M_LIST, _ramp("Reds", len(M_LIST))))


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    full_ped = defaultdict(list)
    static = defaultdict(list)
    pure_scl8 = []
    for r in rows:
        pt = (float(r["snr"]), float(r["bler"]))
        if r["method"] == "full_ped":
            full_ped[64].append(pt)
        elif r["method"] == "static_pm_rank":
            static[int(r["m"])].append(pt)
        elif r["method"] == "pure_scl8":
            pure_scl8.append(pt)

    fig, ax = plt.subplots(figsize=(9, 7))

    pts = sorted(full_ped[64])
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#2a78d6",
                 linestyle="-", marker="o", markersize=6, linewidth=2,
                 label="full PED (M=64)", zorder=3)

    for m in M_LIST:
        pts = sorted(static.get(m, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS[m],
                 linestyle="--", marker="^", markersize=5, linewidth=1.4,
                 label=f"static pm-rank (m={m})", zorder=2)

    pts = sorted(pure_scl8)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#333333",
                 linestyle=":", marker="D", markersize=5, linewidth=1.4,
                 label="pure SCL (L=8)", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Golay(24,12) AED: full-PED vs static-pm-rank")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "golay_m_sweep_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
