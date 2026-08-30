#!/usr/bin/env python3
"""Compare matched-t APS (trained and evaluated at the same checkpoint,
matched_aps.txt) against the mismatched t=64-trained APS applied
out-of-distribution (generalizability.txt), for t in {8,32,96}, both
sides -- ALL in one combined figure. x=SNR, y=BLER. Encoding: color=t
(viridis ramp), marker o=model_out / ^=pm, linestyle solid=matched /
dashed=mismatched -- the gap between a matched pair (same color+marker,
solid vs dashed) is the accuracy lost purely from evaluating a selector
at a checkpoint it wasn't trained for."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"
MO_DIR = HERE

MARKERS = {
    ("pm", "matched"): "^",
    ("pm", "mismatched"): "v",
    ("model_out", "matched"): "o",
    ("model_out", "mismatched"): "s",
}
LINESTYLES = {
    ("model_out", "matched"): "-",
    ("model_out", "mismatched"): "--",
    ("pm", "matched"): "-.",
    ("pm", "mismatched"): ":",
}
T_LIST = [8, 32, 96]
T_COLORS = dict(zip(T_LIST, plt.get_cmap("viridis")(np.linspace(0.1, 0.85, len(T_LIST)))))


def load(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [(int(r["eval_t"]), float(r["snr"]), float(r["bler"])) for r in rows]


def main():
    mismatched = {"pm": load(PM_DIR / "generalizability.txt"), "model_out": load(MO_DIR / "generalizability.txt")}
    matched = {"pm": load(PM_DIR / "matched_aps.txt"), "model_out": load(MO_DIR / "matched_aps.txt")}
    for side in ("pm", "model_out"):
        print(f"{side}: matched={len(matched[side])} points, mismatched(t=64)={len(mismatched[side])} points")

    fig, ax = plt.subplots(figsize=(11, 8))
    for t in T_LIST:
        for side in ("pm", "model_out"):
            mm_pts = sorted((s, b) for tt, s, b in matched[side] if tt == t)
            if mm_pts:
                snrs, blers = zip(*mm_pts)
                ax.plot(snrs, blers, linestyle=LINESTYLES[(side, "matched")], marker=MARKERS[(side, "matched")],
                         markersize=6, linewidth=2, color=T_COLORS[t],
                         label=f"{side} (trained t={t}, evaluate t={t})", zorder=3)
            ms_pts = sorted((s, b) for tt, s, b in mismatched[side] if tt == t)
            if ms_pts:
                snrs, blers = zip(*ms_pts)
                ax.plot(snrs, blers, linestyle=LINESTYLES[(side, "mismatched")], marker=MARKERS[(side, "mismatched")],
                         markersize=6, linewidth=1.4, color=T_COLORS[t],
                         label=f"{side} (trained t=64, evaluate t={t})", zorder=2)

    # t=64 reference (trained@64, evaluate@64 -- trivially "matched" since
    # it's the same selector generalizability.txt already measures at its
    # own training checkpoint). pm only: pm and model_out are known to
    # collide closely, so a second near-identical line just adds clutter.
    ref_pts = sorted((s, b) for tt, s, b in mismatched["pm"] if tt == 64)
    if ref_pts:
        snrs, blers = zip(*ref_pts)
        ax.plot(snrs, blers, linestyle="-", marker="D", markersize=6, linewidth=2.2,
                 color="black", label="pm (trained t=64, evaluate t=64) -- reference", zorder=4)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER (real pickBest-among-top-predicted-m)")
    ax.set_title("APS matched (trained @ eval t) vs mismatched (trained @ t=64) -- color=t, "
                 "marker o/s/^/v = model_out matched/mismatched, pm matched/mismatched, "
                 "line solid/dashed/dashdot/dotted = same order")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="lower left")
    fig.tight_layout()
    out = HERE / "matched_vs_mismatched.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
