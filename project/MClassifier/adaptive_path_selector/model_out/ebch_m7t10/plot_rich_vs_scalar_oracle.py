#!/usr/bin/env python3
"""Quick SNR=2.0-3.5 comparison: scalar m (original design), internal_SC_APS
(m_prefix, weighted + rich trace_learnt-style input), and the real-pickBest
oracle, all at t=64. Data is hardcoded from the individual streaming runs
(smaller, targeted sweep for fast inspection, not the full 8-point trial)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

# (snr, bler)
SCALAR_M = [(2.0, 0.177054), (2.25, 0.114416), (2.5, 0.0707014), (2.75, 0.0400513),
            (3.0, 0.0212965), (3.5, 0.00505449)]
INTERNAL_SC_APS = [(2.0, 0.154943), (2.25, 0.109794), (2.5, 0.057006), (2.75, 0.035073),
                    (3.0, 0.017681), (3.5, 0.00404989)]
ORACLE = [(2.0, 0.0849473), (2.25, 0.0529269), (2.5, 0.0280143), (2.75, 0.0145734),
          (3.0, 0.00706504), (3.5, 0.0013919)]

# pm_vs_model_out_results.txt -- OLD proxy metric (basin-membership fail
# condition), NOT real pickBest. Kept visually distinct (dotted, different
# markers) so it's never confused with the real-pickBest curves above.
PM_ADAPTIVE = [(2.0, 0.171167), (2.25, 0.114500), (2.5, 0.071333), (2.75, 0.042000),
               (3.0, 0.022300), (3.5, 0.004828)]
PM_ORACLE = [(2.0, 0.082167), (2.25, 0.049167), (2.5, 0.029333), (2.75, 0.015167),
             (3.0, 0.008000), (3.5, 0.001334)]
MO_ADAPTIVE = [(2.0, 0.169333), (2.25, 0.108000), (2.5, 0.071667), (2.75, 0.038500),
               (3.0, 0.022400), (3.5, 0.004727)]
MO_ORACLE = [(2.0, 0.082500), (2.25, 0.044500), (2.5, 0.028500), (2.75, 0.013000),
             (3.0, 0.007600), (3.5, 0.001295)]


def main():
    fig, ax = plt.subplots(figsize=(10, 7.5))
    for label, pts, color, marker in (
        ("scalar m (real pickBest)", SCALAR_M, "#a31f8a", "^"),
        ("internal_SC_APS (real pickBest)", INTERNAL_SC_APS, "#1f8a7a", "o"),
        ("oracle (real pickBest)", ORACLE, "black", "D"),
    ):
        s, b = zip(*pts)
        ax.plot(s, b, color=color, marker=marker, markersize=6, linewidth=2, label=label, zorder=3)

    for label, pts, color, marker in (
        ("pm: adaptive (proxy metric)", PM_ADAPTIVE, "#d6622a", "^"),
        ("pm: oracle (proxy metric)", PM_ORACLE, "#2a78d6", "D"),
        ("model_out: adaptive (proxy metric)", MO_ADAPTIVE, "#8a3fa3", "^"),
        ("model_out: oracle (proxy metric)", MO_ORACLE, "#1fa34a", "D"),
    ):
        s, b = zip(*pts)
        ax.plot(s, b, color=color, marker=marker, markersize=5, linewidth=1.3, linestyle=":",
                 label=label, zorder=1, alpha=0.75)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER (real pickBest)")
    ax.set_title("t=64: scalar m vs internal_SC_APS vs oracle, real pickBest (solid) "
                 "vs pm_vs_model_out proxy metric (dotted) -- SNR 2.0-3.5")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, ncol=2, loc="lower left")
    fig.tight_layout()
    out = HERE / "rich_vs_scalar_oracle.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
