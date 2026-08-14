#!/usr/bin/env python3
"""BLER-vs-Eb/N0 plot for the Golay24 t-sweep trial (CUMULATIVE checkpoint
variant): fixed m=2, checkpoint t in {4,8,16,20}, static-pm-rank vs
trace-learned (models trained on cumulative upto-t features, SNR=2.0dB),
vs full-PED (M=16) and pure-SCL(L=8)/pure-SCL(L=2)/pure-SC(L=1) baselines
(all reused from the _single trial), over whatever SNR points are
currently done. Also includes a t=24 static-pm-rank sanity point (n=24 is
the full codeword length, so ranking by the checkpoint at the last decode
position collapses onto ranking by the final metric -- it should collide
with full-PED)."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "sweep" / "golay_t_sweep_m2_cumulative_results.txt"

T_LIST = [4, 8, 16, 20]
GAP_T_LIST = [21, 22, 23]  # static-only probe filling the t=20 -> t=24 tail gap


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


STATIC_COLORS = dict(zip(T_LIST, _ramp("Reds", len(T_LIST))))
TRACE_COLORS = dict(zip(T_LIST, _ramp("Greens", len(T_LIST))))
GAP_COLORS = dict(zip(GAP_T_LIST, _ramp("Oranges", len(GAP_T_LIST))))


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    full_ped, pure_scl8, pure_scl2, pure_sc, static_t24 = [], [], [], [], []
    static, trace = defaultdict(list), defaultdict(list)
    for r in rows:
        pt = (float(r["snr"]), float(r["bler"]))
        if r["method"] == "full_ped":
            full_ped.append(pt)
        elif r["method"] == "pure_scl8":
            pure_scl8.append(pt)
        elif r["method"] == "pure_scl2":
            pure_scl2.append(pt)
        elif r["method"] == "pure_sc":
            pure_sc.append(pt)
        elif r["method"] == "static_pm_rank" and int(r["t"]) == 24:
            static_t24.append(pt)
        elif r["method"] == "static_pm_rank":
            static[int(r["t"])].append(pt)
        elif r["method"] == "trace_learned":
            trace[int(r["t"])].append(pt)

    fig, ax = plt.subplots(figsize=(10, 7.5))

    pts = sorted(full_ped)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#2a78d6",
                 linestyle="-", marker="o", markersize=6, linewidth=2,
                 label="full PED (M=16)", zorder=3)

    pts = sorted(static_t24)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="black",
                 linestyle="none", marker="*", markersize=11, markerfacecolor="none",
                 markeredgewidth=1.6, label="static pm-rank (t=24, sanity check)", zorder=4)

    for t in T_LIST:
        pts = sorted(static.get(t, []))
        if pts:
            ax.plot([s for s, _ in pts], [b for _, b in pts], color=STATIC_COLORS[t],
                     linestyle="--", marker="^", markersize=5, linewidth=1.4,
                     label=f"static pm-rank (t={t}, m=2)", zorder=2)
        pts = sorted(trace.get(t, []))
        if pts:
            ax.plot([s for s, _ in pts], [b for _, b in pts], color=TRACE_COLORS[t],
                     linestyle=":", marker="s", markersize=5, linewidth=1.4,
                     label=f"trace-learned (t={t}, m=2)", zorder=2)

    for t in GAP_T_LIST:
        pts = sorted(static.get(t, []))
        if pts:
            ax.plot([s for s, _ in pts], [b for _, b in pts], color=GAP_COLORS[t],
                     linestyle="--", marker="D", markersize=5, linewidth=1.4,
                     label=f"static pm-rank (t={t}, m=2, gap probe)", zorder=2)

    pts = sorted(pure_scl8)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#333333",
                 linestyle="-.", marker="D", markersize=5, linewidth=1.4,
                 label="pure SCL (L=8)", zorder=2)
    pts = sorted(pure_scl2)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#6a3d9a",
                 linestyle="-.", marker="P", markersize=6, linewidth=1.4,
                 label="pure SCL (L=2)", zorder=2)
    pts = sorted(pure_sc)
    if pts:
        ax.plot([s for s, _ in pts], [b for _, b in pts], color="#8a8a8a",
                 linestyle="-.", marker="x", markersize=6, linewidth=1.4,
                 label="pure SC (L=1)", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Golay(24,12) AED (M=16, list_size=2, m=2, cumulative checkpoints): full-PED vs static-pm-rank vs trace-learned over t")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "golay_t_sweep_m2_cumulative_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
