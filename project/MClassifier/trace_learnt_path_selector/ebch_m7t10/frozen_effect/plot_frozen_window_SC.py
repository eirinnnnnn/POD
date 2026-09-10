#!/usr/bin/env python3
"""BLER-vs-Eb/N0 for two windows, each extended to run from one info bit
before its frozen run to the first info bit after it: t in {14..20}
(decode_idx 13-15 info, 16/17/18 frozen [span 4/1/1], 19 info, 20 frozen
again) and t in {22..28} (decode_idx 21-23 info, 24/25/26 frozen [span
4/5/4], 27 info, 28 frozen again). static_pm_rank only (no model), k=8.
Checkpoint t means positions 0..t-1 are decided; every t whose
newly-consumed bit t-1 is frozen is marked "(frozen)" in the legend."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "frozen_window_static_results_SC.txt"

WINDOW_16 = [14, 15, 16, 17, 18, 19, 20]
WINDOW_24 = [22, 23, 24, 25, 26, 27, 28]
SNR_MIN, SNR_MAX = 2.00, 4.50  # extended range for the SC (list_size=1) trial

# checkpoint t consumes position t-1 -- from ebch_dynamic_frozen_matrix.txt,
# decode_idx 16/17/18 and 24/25/26 are each 3 consecutive frozen bits, each
# run bounded by a single info bit on either side (19 and 27 respectively).
# "(frozen)" marks every t whose newly-consumed bit (t-1) is frozen.
FROZEN_AT = {17: "(frozen)", 18: "(frozen)", 19: "(frozen)",
             25: "(frozen)", 26: "(frozen)", 27: "(frozen)"}


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
        after = t in FROZEN_AT  # newly-consumed bit (t-1) is frozen
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS_16[t],
                 linestyle="-" if after else "--",
                 marker="s" if after else "o",
                 markersize=6 if after else 5,
                 linewidth=2.2 if after else 1.3,
                 label=f"static pm-rank (m=8, t={t}) {FROZEN_AT.get(t, '')}".rstrip(),
                 zorder=5 if after else 3)
    ax.set_yscale("log")
    ax.set_xlim(SNR_MIN - x_pad, SNR_MAX + x_pad)
    ax.set_xticks([2.00, 2.50, 3.00, 3.50, 4.00, 4.50])
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("t=14..20 (idx16/17/18 frozen, idx19 info, idx20 frozen again)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left", framealpha=0.9)

    ax = axes[1]
    for t in WINDOW_24:
        pts = sorted(series.get(t, []))
        if not pts:
            continue
        after = t in FROZEN_AT  # newly-consumed bit (t-1) is frozen
        ax.plot([s for s, _ in pts], [b for _, b in pts], color=COLORS_24[t],
                 linestyle="-" if after else "--",
                 marker="s" if after else "o",
                 markersize=6 if after else 5,
                 linewidth=2.2 if after else 1.3,
                 label=f"static pm-rank (m=8, t={t}) {FROZEN_AT.get(t, '')}".rstrip(),
                 zorder=5 if after else 3)
    ax.set_yscale("log")
    ax.set_xlim(SNR_MIN - x_pad, SNR_MAX + x_pad)
    ax.set_xticks([2.00, 2.50, 3.00, 3.50, 4.00, 4.50])
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_title("t=22..28 (idx24/25/26 frozen, idx27 info, idx28 frozen again)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left", framealpha=0.9)

    fig.suptitle("(128,64) eBCH, static pm-rank only (k=8, no model, SCL list_size=1) -- frozen-bit t±2 windows")
    fig.tight_layout()

    out = HERE / "frozen_window_bler_SC.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
