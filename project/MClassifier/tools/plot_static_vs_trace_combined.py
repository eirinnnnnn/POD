#!/usr/bin/env python3
"""Combined k-vs-BLER and k-vs-basin-hit plot across all trace checkpoints:
one color per checkpoint t, solid = learned trace classifier, dashed =
model-free static path-metric ranking (pm_min sort, no training)."""
import csv
import subprocess

import matplotlib.pyplot as plt
import numpy as np

DATASET = "data/m7t10_m64_snr2p0_test_metric_3k.csv"
KS = [1, 2, 4, 8, 16, 32, 64]
STEPS = [8, 16, 32, 64, 128]

# User-requested colors. Note: red/orange (t=64/32) fail the palette
# validator's adjacent-pair CVD and normal-vision checks (ΔE well below
# target) -- compensating with distinct marker shapes per t below.
COLORS = {
    128: "#2a78d6",  # blue
    64: "#e34948",   # red
    32: "#eb6834",   # orange
    16: "#008300",   # green
    8: "#4a3aa7",    # purple
}
MARKERS = {128: "o", 64: "^", 32: "D", 16: "s", 8: "v"}

LEARNED_PRED = {
    8: "build/m7t10_m64_snr2p0_trace_upto_t8_basin_h64_e20_top64_3k.csv",
    16: "build/m7t10_m64_snr2p0_trace_upto_t16_basin_h64_e20_top64_3k.csv",
    32: "build/m7t10_m64_snr2p0_trace_upto_t32_basin_h64_e20_top64_3k.csv",
    64: "build/m7t10_m64_snr2p0_trace_upto_t64_basin_h64_e20_top64_3k.csv",
    128: "build/m7t10_m64_snr2p0_trace_basin_h64_e20_top64_3k.csv",
}
STATIC_PRED = {t: f"build/m7t10_m64_snr2p0_static_pm_rank_t{t}_top64_3k.csv" for t in STEPS}
MAJORITY_PRED = "build/m7t10_m64_snr2p0_majority_train_top64_3k.csv"


def load_ground_truth():
    rows = list(csv.DictReader(open(DATASET, encoding="utf-8")))
    M = len([c for c in rows[0] if c.startswith("metric_") and c[7:].isdigit()])
    metrics = np.array([[float(r[f"metric_{i}"]) for i in range(M)] for r in rows])
    basin = np.abs(metrics - metrics.min(axis=1, keepdims=True)) <= 1e-12
    return M, basin


def curve(pred_path, M, basin):
    pred_rows = list(csv.DictReader(open(pred_path, encoding="utf-8")))
    pred_order = np.array([[int(r[f"idx{i}"]) for i in range(M)] for r in pred_rows])
    blers, hits = [], []
    for k in KS:
        out = subprocess.check_output(
            ["./build/mclass_eval", "--dataset", DATASET, "--pred", pred_path, "--topk", str(k)],
            text=True,
        )
        vals = dict(line.split(",") for line in out.strip().splitlines())
        top = pred_order[:, :k]
        hit = basin[np.arange(basin.shape[0])[:, None], top].any(axis=1).mean()
        blers.append(float(vals["topk_bler"]))
        hits.append(hit)
    return blers, hits


def main():
    M, basin = load_ground_truth()

    majority_bler, majority_hit = curve(MAJORITY_PRED, M, basin)
    learned = {t: curve(LEARNED_PRED[t], M, basin) for t in STEPS}
    static = {t: curve(STATIC_PRED[t], M, basin) for t in STEPS}

    # --- BLER plot ---
    plt.figure(figsize=(7.5, 4.5))
    plt.semilogy(KS, majority_bler, color="#8a8a86", linestyle=":", marker="o",
                 markersize=4, linewidth=1, label="majority")
    for t in STEPS:
        c, m = COLORS[t], MARKERS[t]
        plt.semilogy(KS, learned[t][0], color=c, linestyle="-", marker=m,
                     markersize=4, linewidth=1, label=f"trace t={t}")
        plt.semilogy(KS, static[t][0], color=c, linestyle="--", marker=m,
                     markersize=4, linewidth=1, markerfacecolor="none", label=f"static pm_rank t={t}")
    plt.xscale("log", base=2)
    plt.xticks(KS, [str(k) for k in KS])
    plt.xlabel(f"decoded AED branches k out of M={M}")
    plt.ylabel("BLER")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=7.5, ncol=2, loc="upper right")
    plt.title("Learned trace classifier vs. model-free path-metric ranking")
    plt.tight_layout()
    plt.savefig("build/m7t10_m64_snr2p0_static_vs_trace_combined_bler.png", dpi=200)
    print("wrote build/m7t10_m64_snr2p0_static_vs_trace_combined_bler.png")

    # --- basin-hit plot ---
    plt.figure(figsize=(7.5, 4.5))
    plt.plot(KS, majority_hit, color="#8a8a86", linestyle=":", marker="o",
              markersize=4, linewidth=1, label="majority")
    for t in STEPS:
        c, m = COLORS[t], MARKERS[t]
        plt.plot(KS, learned[t][1], color=c, linestyle="-", marker=m,
                  markersize=4, linewidth=1, label=f"trace t={t}")
        plt.plot(KS, static[t][1], color=c, linestyle="--", marker=m,
                  markersize=4, linewidth=1, markerfacecolor="none", label=f"static pm_rank t={t}")
    plt.xscale("log", base=2)
    plt.xticks(KS, [str(k) for k in KS])
    plt.xlabel(f"predicted top-k branches out of M={M}")
    plt.ylabel("best-basin hit rate")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=7.5, ncol=2, loc="lower right")
    plt.title("Learned trace classifier vs. model-free path-metric ranking")
    plt.tight_layout()
    plt.savefig("build/m7t10_m64_snr2p0_static_vs_trace_combined_basin_hit.png", dpi=200)
    print("wrote build/m7t10_m64_snr2p0_static_vs_trace_combined_basin_hit.png")


if __name__ == "__main__":
    main()
