#!/usr/bin/env python3
"""Plot trial1_results.txt as BLER vs Eb/N0, one combined figure, over
whatever SNR points are currently available. Four families: full_ped
(genuine M-branch ensemble, no pruning), static_pm_rank / trace_learned
(AED M=64 pruned to k, model-free / learned), pure_scl (no AED at all,
list size L). Color = family (sequential shade by k/L/M within family),
marker/linestyle = family identity, so no series relies on color alone."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "trial1_results.txt"

FULL_PED_M = [1, 4, 32, 64]
STATIC_K = [1, 4, 8, 32]
TRACE_K = [1, 4, 8, 32]
SCL_L = [4, 16, 32, 128]


def _ramp(cmap_name, n):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(0.35 + 0.55 * i / max(n - 1, 1)) for i in range(n)]


FULL_PED_COLORS = dict(zip(FULL_PED_M, _ramp("Blues", len(FULL_PED_M))))
STATIC_COLORS = dict(zip(STATIC_K, _ramp("Reds", len(STATIC_K))))
TRACE_COLORS = dict(zip(TRACE_K, _ramp("Greens", len(TRACE_K))))
SCL_COLORS = dict(zip(SCL_L, _ramp("Purples", len(SCL_L))))


def series_style(config, method, k):
    if method == "full_ped":
        return dict(color=FULL_PED_COLORS[k], linestyle="-", marker="o", markersize=6,
                    linewidth=1.6, label=f"full PED (M={k}, no pruning)")
    if method == "static_pm_rank":
        return dict(color=STATIC_COLORS[k], linestyle="-", marker="^", markersize=6,
                    linewidth=1.6, label=f"AED M=64 prune to k={k} (static)")
    if method == "trace_learned":
        return dict(color=TRACE_COLORS[k], linestyle="--", marker="s", markersize=6,
                    linewidth=1.6, label=f"AED M=64 prune to k={k} (trace-learned)")
    return dict(color=SCL_COLORS[k], linestyle=":", marker="D", markersize=6,
                linewidth=1.6, label=f"pure SCL, L={k} (no AED)")


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    series = defaultdict(list)  # (method, k) -> [(snr, bler)]
    for r in rows:
        series[(r["method"], int(r["k"]))].append((float(r["snr"]), float(r["bler"])))

    order = [("full_ped", m) for m in FULL_PED_M]
    order += [("static_pm_rank", k) for k in STATIC_K]
    order += [("trace_learned", k) for k in TRACE_K]
    # pure_scl excluded -- mclass_bler_sim was forcing use_AED=true regardless
    # of the ini, so every "pure_scl" run so far actually ran AED M=64 + a
    # big SCL list, not plain SCL. Fixed in the C++ source; re-add once the
    # trial reruns with the corrected binary.

    fig, ax = plt.subplots(figsize=(10, 7.5))
    for method, k in order:
        pts = sorted(series.get((method, k), []))
        if not pts:
            continue
        style = series_style(None, method, k)
        ax.plot([s for s, _ in pts], [b for _, b in pts], **style)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("m7t10_m64 trial 1: BLER vs Eb/N0 (currently available points)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, ncol=2, loc="lower left", framealpha=0.9)
    fig.tight_layout()

    out = HERE / "trial1_snr.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
