#!/usr/bin/env python3
"""Merged single-figure plot of the Eb/N0 sweep: BLER vs SNR for every
method on one axis. Color and linestyle are both keyed by checkpoint t and
SHARED between static-pm-rank and trace-learned at the same t, so a same-t
pair is directly comparable by color; marker shape instead carries method
identity (triangle=static-pm-rank, square=trace-learned) consistently
across every t. Every series is also legend-labelled, so identity is never
color-alone."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.txt"

GRAY = "#52514e"

ALL_T = [8, 32, 64, 80, 96]
TRACE_T = [8, 32, 64]

# One categorical hue + linestyle per checkpoint t, shared across methods.
T_COLOR = {
    8: "#2a78d6",    # blue
    32: "#eb6834",   # orange
    64: "#1baf7a",   # aqua
    80: "#e87ba4",   # magenta
    96: "#4a3aa7",   # violet
}
# linestyle (and marker) now carry method identity, consistent across every
# t; only color varies by t.
METHOD_LINESTYLE = {"static_pm_rank": "-", "trace_learned": "--"}
METHOD_MARKER = {"static_pm_rank": "^", "trace_learned": "s"}


def series_style(method, t):
    if method == "full_ped":
        return dict(color="black", linestyle="--", marker="None", linewidth=1.2,
                    label="full PED (M=64)", zorder=5)
    if method == "full_ped_m8":
        return dict(color=GRAY, linestyle="--", marker="D", markersize=7,
                    linewidth=2.0, label="full PED (M=8)", zorder=4)
    label = "static pm-rank" if method == "static_pm_rank" else "trace-learned"
    return dict(color=T_COLOR[t], linestyle=METHOD_LINESTYLE[method], marker=METHOD_MARKER[method],
                markersize=7, linewidth=1.6, label=f"{label} (k=8, t={t})", zorder=3)


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    series = defaultdict(list)  # (method, t) -> [(snr, bler)]
    for r in rows:
        t = r["t"] if r["t"] == "all" else int(r["t"])
        series[(r["method"], t)].append((float(r["snr"]), float(r["bler"])))

    fig, ax = plt.subplots(figsize=(10, 7.5))

    order = [("full_ped", "all"), ("full_ped_m8", "all")]
    for t in ALL_T:
        order.append(("static_pm_rank", t))
    for t in TRACE_T:
        order.append(("trace_learned", t))

    for method, t in order:
        pts = sorted(series[(method, t)])
        if not pts:
            continue
        style = series_style(method, t if t != "all" else None)
        snrs = [s for s, _ in pts]
        blers = [b for _, b in pts]
        ax.plot(snrs, blers, **style)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("m7t10_m64: full-PED (M=64, M=8) vs static-pm-rank vs trace-learned (k=8), SCL L=4")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9, ncol=2, loc="lower left", framealpha=0.9)

    fig.tight_layout()
    out = HERE / "ebn0_sweep_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
