#!/usr/bin/env python3
"""BLER-vs-Eb/N0 for the two frozen-bit windows t in {14..18} (centered on
the span-4 dynamic-frozen bit at decode_idx=16) and t in {22..26} (centered
on the span-4 one at decode_idx=24). static_pm_rank only (no model), k=8.
Checkpoint t means positions 0..t-1 are decided, so idx=16's constraint
first lands at t=17 (dashed->solid within the 16-window) and idx=24's at
t=25 (dashed->solid within the 24-window)."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "frozen_window_static_results.txt"

WINDOW_16 = [14, 15, 16, 17, 18]
WINDOW_24 = [22, 23, 24, 25, 26]
SNR_MIN, SNR_MAX = 2.00, 3.50  # trial scope -- excludes a stray t=14/SNR=4.0
# row left over from before the sweep was capped at 3.5


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.3 + 0.6 * i / max(n - 1, 1)) for i in range(n)]


COLORS_16 = dict(zip(WINDOW_16, _ramp("Reds", len(WINDOW_16))))
COLORS_24 = dict(zip(WINDOW_24, _ramp("Blues", len(WINDOW_24))))


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    series = defaultdict(list)
    all_blers = []
    for r in rows:
        snr = float(r["snr"])
        if not (SNR_MIN <= snr <= SNR_MAX):
            continue
        bler = float(r["bler"])
        series[int(r["t"])].append((snr, bler))
        all_blers.append(bler)

    y_lo = min(all_blers) * 0.8
    y_hi = max(all_blers) * 1.25
    x_pad = 0.06

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)

    ax = axes[0]
    for t in WINDOW_16:
        pts = sorted(series.get(t, []))
        if not pts:
            continue
        after = t >= 17  # first checkpoint that has resolved decode_idx=16
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS_16[t],
                 linestyle="-" if after else "--",
                 marker="o" if t == 16 else ("s" if t == 17 else "^"),
                 markersize=6 if t in (16, 17) else 5,
                 linewidth=2.2 if t in (16, 17) else 1.3,
                 label=f"static pm-rank (m=8, t={t})" + ("  <- before idx16" if t == 16 else
                                                            "  <- after idx16 (span=4)" if t == 17 else ""),
                 zorder=5 if t in (16, 17) else 3)
    ax.set_yscale("log")
    ax.set_xlim(SNR_MIN - x_pad, SNR_MAX + x_pad)
    ax.set_xticks([2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50])
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Window centered on decode_idx=16 (frozen, span=4, [13,14,15,16])")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left", framealpha=0.9)

    ax = axes[1]
    for t in WINDOW_24:
        pts = sorted(series.get(t, []))
        if not pts:
            continue
        after = t >= 25  # first checkpoint that has resolved decode_idx=24
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS_24[t],
                 linestyle="-" if after else "--",
                 marker="o" if t == 24 else ("s" if t == 25 else "^" if t == 26 else "^"),
                 markersize=6 if t in (24, 25, 26) else 5,
                 linewidth=2.2 if t in (24, 25, 26) else 1.3,
                 label=f"static pm-rank (m=8, t={t})" + ("  <- before idx24" if t == 24 else
                                                            "  (after idx24, span=4)" if t == 25 else
                                                            "  (after idx25, span=5)" if t == 26 else ""),
                 zorder=5 if t in (24, 25, 26) else 3)
    ax.set_yscale("log")
    ax.set_xlim(SNR_MIN - x_pad, SNR_MAX + x_pad)
    ax.set_xticks([2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50])
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_title("Window centered on decode_idx=24 (frozen, span=4, [21,22,23,24])")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left", framealpha=0.9)

    fig.suptitle("(128,64) eBCH, static pm-rank only (k=8, no model) -- frozen-bit t±2 windows")
    fig.tight_layout()

    out = HERE / "frozen_window_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
