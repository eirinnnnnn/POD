#!/usr/bin/env python3
"""Combined BLER-vs-Eb/N0 plot for the t_sweep_m8 trial (clean automorphism
file, m=8 pruning width fixed, t in {4,8,16,32,64}), over whatever SNR
points are currently done. full_ped (M=64) is the reference; static-pm-rank/
trace-learned at t=4,8,16,32,64, m=8. full_ped_m8 rows (leftover from before
that step was dropped) are excluded. Also overlays the t_sweep_m8_SCL_baseline
trial (pure SCL, no AED, L=65,68,76,96) and the PED_M-SCL4 baseline (genuine
M-branch ensemble, SCL list_size=4, M chosen so M*4~=L) as reference lines."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "t_sweep_m8_results.txt"
SCL_RESULTS = HERE / "t_sweep_m8_SCL_baseline_results.txt"
PED_RESULTS = HERE / "t_sweep_m8_PED_baseline_results.txt"

STATIC_T = [4, 8, 16, 32, 64]
TRACE_T = [4, 8, 16, 32, 64]
SCL_L = [65, 68, 76, 96]
PED_M_TO_L = {8: 32, 16: 65, 17: 68, 19: 76, 24: 96, 36: 144}
PED_M = list(PED_M_TO_L)


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


STATIC_COLORS = dict(zip(STATIC_T, _ramp("Reds", len(STATIC_T))))
TRACE_COLORS = dict(zip(TRACE_T, _ramp("Greens", len(TRACE_T))))
SCL_COLORS = dict(zip(SCL_L, _ramp("Blues", len(SCL_L))))
PED_COLORS = dict(zip(PED_M, _ramp("Purples", len(PED_M))))


def series_style(method, t):
    if method == "static_pm_rank":
        return dict(color=STATIC_COLORS[t], linestyle="-", marker="^", markersize=6,
                    linewidth=1.6, label=f"static pm-rank (m=8, t={t})", zorder=3)
    return dict(color=TRACE_COLORS[t], linestyle="--", marker="s", markersize=6,
                linewidth=1.6, label=f"trace-learned (m=8, t={t})", zorder=3)


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["method"] not in ("full_ped", "full_ped_m8")]

    series = defaultdict(list)
    for r in rows:
        t = r["t"] if r["t"] == "all" else int(r["t"])
        series[(r["method"], t)].append((float(r["snr"]), float(r["bler"])))

    scl_series = defaultdict(list)
    if SCL_RESULTS.exists():
        with open(SCL_RESULTS, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                scl_series[int(r["L"])].append((float(r["snr"]), float(r["bler"])))

    ped_series = defaultdict(list)
    if PED_RESULTS.exists():
        with open(PED_RESULTS, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ped_series[int(r["M"])].append((float(r["snr"]), float(r["bler"])))

    fig, ax = plt.subplots(figsize=(10, 7.5))
    order = []
    for t in STATIC_T:
        order.append(("static_pm_rank", t))
    for t in TRACE_T:
        order.append(("trace_learned", t))

    for method, t in order:
        pts = sorted(series.get((method, t), []))
        if not pts:
            continue
        style = series_style(method, t if t != "all" else None)
        ax.plot([s for s, _ in pts], [b for _, b in pts], **style)

    for L in SCL_L:
        pts = sorted(scl_series.get(L, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=SCL_COLORS[L],
                 linestyle=":", marker="o", markersize=5, linewidth=1.4,
                 label=f"pure SCL (L={L})", zorder=2)

    for M in PED_M:
        pts = sorted(ped_series.get(M, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=PED_COLORS[M],
                 linestyle="-.", marker="D", markersize=5, linewidth=1.4,
                 label=f"PED_{M}-SCL4 (L~{PED_M_TO_L[M]})", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("(128, 64) eBCH M=64, m=8. Static Path Metric Ranker vs Trace Learned Selector, over different t")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "t_sweep_m8_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
