#!/usr/bin/env python3
"""Combined BLER-vs-Eb/N0 plot for the m_sweep_t8 trial ('M choose m'
early-path-selecting-PED sweep): fixed checkpoint t=8, fixed SCL list_size
L=4, pruning width m in {1,4,16,32}, static-pm-rank vs trace-learned, over
whatever SNR points are currently done. Also overlays the
m_sweep_t8_SCL_baseline trial (pure SCL, no AED, L=45,55,95,149) and two
PED_x-SCL4 baselines (genuine x-branch ensemble, SCL list_size=4): a
"complexity match" set where M*4 equals a target total-path-complexity
budget L (L=200/M=50 is the complexity-matched baseline for m=48), and an
"m match" set where M exactly equals one of the m-sweep's own m values."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "m_sweep_t8_results.txt"
SCL_RESULTS = HERE / "m_sweep_t8_SCL_baseline_results.txt"
PED_RESULTS = HERE / "m_sweep_t8_PED_baseline_results.txt"

M_LIST = [1, 4, 16, 32, 48]
SCL_L = [45, 55, 95, 149]
# complexity target L -> PED baseline ensemble size M (M*4~=L)
COMPLEXITY_L_TO_M = {45: 11, 55: 14, 95: 24, 149: 37, 200: 50}
PED_L_MATCHED_M = list(COMPLEXITY_L_TO_M.values())
# m-sweep value -> PED baseline ensemble size M, exact match (M==m)
M_MATCH_MAP = {1: 1, 4: 4, 16: 16, 32: 32, 48: 48}


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


STATIC_COLORS = dict(zip(M_LIST, _ramp("Reds", len(M_LIST))))
TRACE_COLORS = dict(zip(M_LIST, _ramp("Greens", len(M_LIST))))
SCL_COLORS = dict(zip(SCL_L, _ramp("Blues", len(SCL_L))))
PED_L_COLORS = dict(zip(PED_L_MATCHED_M, _ramp("Purples", len(PED_L_MATCHED_M))))
PED_M_COLORS = dict(zip(M_MATCH_MAP.values(), _ramp("Oranges", len(M_MATCH_MAP))))


def series_style(method, m):
    if method == "static_pm_rank":
        return dict(color=STATIC_COLORS[m], linestyle="-", marker="^", markersize=6,
                    linewidth=1.6, label=f"static pm-rank (t=8, m={m})", zorder=3)
    return dict(color=TRACE_COLORS[m], linestyle="--", marker="s", markersize=6,
                linewidth=1.6, label=f"trace-learned (t=8, m={m})", zorder=3)


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    series = defaultdict(list)
    for r in rows:
        m = int(r["m"])
        series[(r["method"], m)].append((float(r["snr"]), float(r["bler"])))

    scl_series = defaultdict(list)
    if SCL_RESULTS.exists():
        with open(SCL_RESULTS, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                scl_series[int(r["L"])].append((float(r["snr"]), float(r["bler"])))

    ped_l_series = defaultdict(list)
    ped_m_series = defaultdict(list)
    if PED_RESULTS.exists():
        with open(PED_RESULTS, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                pt = (float(r["snr"]), float(r["bler"]))
                if r["kind"] == "L_matched":
                    ped_l_series[int(r["M"])].append(pt)
                else:
                    ped_m_series[int(r["M"])].append(pt)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    order = [("static_pm_rank", m) for m in M_LIST] + [("trace_learned", m) for m in M_LIST]

    for method, m in order:
        pts = sorted(series.get((method, m), []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], **series_style(method, m))

    for L in SCL_L:
        pts = sorted(scl_series.get(L, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=SCL_COLORS[L],
                 linestyle=":", marker="o", markersize=5, linewidth=1.4,
                 label=f"pure SCL (L={L})", zorder=2)

    M_TO_L = {M: L for L, M in COMPLEXITY_L_TO_M.items()}
    for M in PED_L_MATCHED_M:
        pts = sorted(ped_l_series.get(M, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=PED_L_COLORS[M],
                 linestyle="-.", marker="D", markersize=5, linewidth=1.4,
                 label=f"PED_{M}-SCL4 (complexity match, L~{M_TO_L[M]})", zorder=2)

    for m, M in M_MATCH_MAP.items():
        pts = sorted(ped_m_series.get(M, []))
        if not pts:
            continue
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=PED_M_COLORS[M],
                 linestyle="-.", marker="v", markersize=5, linewidth=1.4,
                 label=f"PED_{M}-SCL4 (m-matched)", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("(128, 64) eBCH M=64, t=8. Static Path Metric Ranker vs Trace Learned Selector, over different m")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "m_sweep_t8_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
