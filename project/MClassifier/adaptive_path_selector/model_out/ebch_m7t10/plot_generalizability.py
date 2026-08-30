#!/usr/bin/env python3
"""ONE combined panel, pm and model_out overlaid directly, from the new
generalizability.txt trial (t in {8,32,64,96}, target-errors=500, streaming
only). Same color per checkpoint t, solid=model_out vs dashed=pm, so the
gap between a matched pair at the same t is the model_out-vs-pm gain,
read directly off the plot. Same style as plot_csv.py, different source
data (this trial's own generalizability.txt, not the older
cross_t_stream_results.txt)."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PM_RESULTS = HERE.parent.parent / "pm" / "ebch_m7t10" / "generalizability.txt"
MO_RESULTS = HERE / "generalizability.txt"


def load(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [(int(r["eval_t"]), float(r["snr"]), float(r["bler"])) for r in rows]


def main():
    pm_rows = load(PM_RESULTS)
    mo_rows = load(MO_RESULTS)
    print(f"pm: {len(pm_rows)} points, model_out: {len(mo_rows)} points")

    t_list = sorted(set(t for t, _, _ in pm_rows) | set(t for t, _, _ in mo_rows))
    t_colors = dict(zip(t_list, plt.get_cmap("viridis")(np.linspace(0.05, 0.85, len(t_list)))))

    fig, ax = plt.subplots(figsize=(11, 8))
    for t in t_list:
        mo_pts = sorted((s, b) for tt, s, b in mo_rows if tt == t)
        if mo_pts:
            snrs, blers = zip(*mo_pts)
            ax.plot(snrs, blers, linestyle="-", marker="o", markersize=6, linewidth=2,
                     color=t_colors[t], label=f"model_out t={t}", zorder=3)
        pm_pts = sorted((s, b) for tt, s, b in pm_rows if tt == t)
        if pm_pts:
            snrs, blers = zip(*pm_pts)
            ax.plot(snrs, blers, linestyle="--", marker="^", markersize=6, linewidth=1.6,
                     color=t_colors[t], label=f"pm t={t}", zorder=2)

    ax.set_yscale("log")
    ax.set_xlim(2.0, 4.5)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER (real pickBest-among-top-predicted-m)")
    ax.set_title("model_out vs pm generalization (target-errors=500 trial) -- solid=model_out, dashed=pm")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    fig.tight_layout()

    out = HERE / "generalizability.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
