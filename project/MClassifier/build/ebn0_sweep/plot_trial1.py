#!/usr/bin/env python3
"""Plot trial1_results.txt: BLER vs path/list budget k, one subplot per SNR.
Four families: full_ped (genuine M-branch AED ensemble, no pruning, M=k),
static_pm_rank / trace_learned (prune M=64 down to k, model-free / learned),
and pure_scl (no AED at all, plain SCL with list size L=k)."""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "trial1_results.txt"

STYLE = {
    "full_ped": dict(color="#2a78d6", marker="o", linestyle="-", label="full PED (M=k, no pruning)"),
    "static_pm_rank": dict(color="#e34948", marker="^", linestyle="-", label="AED M=64, prune to k (static pm-rank)"),
    "trace_learned": dict(color="#008300", marker="s", linestyle="--", label="AED M=64, prune to k (trace-learned)"),
    # pure_scl excluded -- mclass_bler_sim was forcing use_AED=true regardless
    # of the ini, so every "pure_scl" run so far actually ran AED M=64 + a
    # big SCL list, not plain SCL. Fixed in the C++ source; re-add once the
    # trial reruns with the corrected binary.
}


def main():
    with open(RESULTS, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    snrs = sorted({float(r["snr"]) for r in rows})
    series = defaultdict(lambda: defaultdict(list))  # snr -> method -> [(k, bler)]
    for r in rows:
        series[float(r["snr"])][r["method"]].append((int(r["k"]), float(r["bler"])))

    fig, axes = plt.subplots(1, len(snrs), figsize=(6.5 * len(snrs), 5.5), sharey=True)
    if len(snrs) == 1:
        axes = [axes]

    for ax, snr in zip(axes, snrs):
        for method, style in STYLE.items():
            pts = sorted(series[snr].get(method, []))
            if not pts:
                continue
            ks = [k for k, _ in pts]
            blers = [b for _, b in pts]
            ax.plot(ks, blers, **style)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("path / list budget k")
        ax.set_title(f"Eb/N0 = {snr} dB")
        ax.grid(True, which="both", alpha=0.3)

    axes[0].set_ylabel("BLER")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("m7t10_m64 trial 1: BLER vs path budget -- prune-a-big-ensemble vs run-a-small-ensemble vs bigger-SCL-list")
    fig.tight_layout()

    out = HERE / "trial1.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
