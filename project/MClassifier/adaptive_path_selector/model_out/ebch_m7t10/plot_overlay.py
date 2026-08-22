#!/usr/bin/env python3
"""Overlay pm's and model_out's adaptive-path-selector curves (t=64) on one
BLER-vs-Eb/N0 plot + one mean-m plot. Mirrors
adaptive_path_selector/pm/ebch_m7t10/plot_adaptive_sweep.py's logic (same
score() definition: fail iff any_basin==0 or used_m < true
m_required_basin) for both sides, but pm and model_out each get their OWN
oracle curve -- m_required_basin is defined relative to whichever ranking
produced the sorted feature vector (pm_min order vs trace-prob order), so
the two oracles are not the same ground-truth number, just the same KIND
of upper bound for each ranking scheme.

Re-runnable at any point mid-sweep: only plots SNR points where model_out's
sweep_snr*.csv has reached its FULL target sample count (an in-progress
file for the current SNR point is silently skipped, not partially
plotted, so an early call never shows a truncated/misleadingly-noisy
point)."""
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"
PM_MODEL_PT = PM_DIR / "data" / "model_t64_snr2.pt"
MO_MODEL_PT = HERE / "data" / "model_t64_snr2.pt"

# SNR -> target sample count for model_out's sweep (must match
# run_sweep_gen.py's SNR_SAMPLES exactly, so a mid-write file is detected
# as incomplete and skipped rather than plotted early/misleadingly).
MO_TARGET_SAMPLES = {
    2.00: 6_000, 2.25: 6_000, 2.50: 6_000, 2.75: 6_000,
    3.00: 10_000, 3.50: 1_000_000, 4.00: 1_000_000, 4.50: 1_000_000,
}


def load_eval(path, feature_marker):
    rows = list(csv.DictReader(open(path)))
    if rows and any(v is None for v in rows[-1].values()):
        rows = rows[:-1]
    cols = [c for c in rows[0].keys() if feature_marker in c]
    cols.sort(key=lambda c: int(c.rsplit("_", 1)[1]))
    X = np.array([[float(r[c]) for c in cols] for r in rows], dtype=np.float32)
    y = np.array([int(r["m_required_basin"]) for r in rows], dtype=np.float64)
    any_basin = np.array([int(r["any_basin_exists"]) for r in rows], dtype=np.float64)
    return X, y, any_basin, len(rows)


def predict(ckpt, X):
    model = nn.Sequential(
        nn.Linear(ckpt["M"], ckpt["hidden"]), nn.ReLU(),
        nn.Linear(ckpt["hidden"], ckpt["hidden"]), nn.ReLU(),
        nn.Linear(ckpt["hidden"], 1),
    )
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    Xt = np.log1p(X) if ckpt.get("log_input", True) else X
    Xn = (Xt - ckpt["x_mean"]) / ckpt["x_std"]
    with torch.no_grad():
        pred = model(torch.from_numpy(Xn.astype(np.float32))).reshape(-1).numpy()
    pred = pred * ckpt["y_std"] + ckpt["y_mean"]
    return np.clip(np.ceil(pred), 1, ckpt["M"])


def score(used_m, y, any_basin):
    fail = (any_basin == 0) | (used_m < y)
    return float(used_m.mean()), float(fail.mean())


def snr_from_name(f):
    m = re.search(r"sweep_snr(\d+)\.csv", f.name)
    digits = m.group(1)
    return float(digits[0] + "." + digits[1:])


def collect(sweep_dir, ckpt_path, feature_marker, target_samples=None):
    files = sorted(sweep_dir.glob("sweep_snr*.csv"))
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    per_snr = []
    for f in files:
        snr = snr_from_name(f)
        X, y, any_basin, n = load_eval(f, feature_marker)
        if target_samples is not None and n < target_samples.get(snr, 0):
            print(f"  skipping SNR={snr} ({sweep_dir.name}): only {n} samples so far "
                  f"(target {target_samples.get(snr)}), still in progress")
            continue
        pred_m = predict(ckpt, X)
        am, ab = score(pred_m, y, any_basin)
        oracle_used = np.where(y == 0, ckpt["M"], y)
        om, ob = score(oracle_used, y, any_basin)
        per_snr.append((snr, am, ab, om, ob))
    return sorted(per_snr), ckpt["M"]


def write_results_txt(path, pm_snr, mo_snr):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["side", "snr", "method", "mean_m", "bler"])
        for side, rows in (("pm", pm_snr), ("model_out", mo_snr)):
            for snr, am, ab, om, ob in rows:
                w.writerow([side, snr, "adaptive", f"{am:.6f}", f"{ab:.6f}"])
                w.writerow([side, snr, "oracle", f"{om:.6f}", f"{ob:.6f}"])
    print(f"wrote {path}")


def main():
    print("pm sweep:")
    pm_snr, pm_M = collect(PM_DIR / "data", PM_MODEL_PT, "_pm_min_sorted_")
    print("model_out sweep:")
    mo_snr, mo_M = collect(HERE / "data", MO_MODEL_PT, "_sorted_", MO_TARGET_SAMPLES)

    write_results_txt(HERE / "pm_vs_model_out_results.txt", pm_snr, mo_snr)

    if not mo_snr:
        print("no complete model_out SNR points yet -- nothing to overlay")
        return

    fig, ax = plt.subplots(figsize=(10.5, 7))
    if pm_snr:
        s, am, ab, om, ob = zip(*pm_snr)
        ax.plot(s, ob, color="#2a78d6", linestyle="-", marker="o", markersize=6,
                linewidth=2, label="pm: oracle (ground truth required m)", zorder=3)
        ax.plot(s, ab, color="#d6622a", linestyle="--", marker="^", markersize=6,
                linewidth=1.6, label="pm: adaptive path selector", zorder=2)
    s, am, ab, om, ob = zip(*mo_snr)
    ax.plot(s, ob, color="#1fa34a", linestyle="-", marker="o", markersize=6,
            linewidth=2, label="model_out: oracle (ground truth required m)", zorder=3)
    ax.plot(s, ab, color="#a31f8a", linestyle="--", marker="^", markersize=6,
            linewidth=1.6, label="model_out: adaptive path selector", zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    n_mo = len(mo_snr)
    ax.set_title(f"pm vs model_out adaptive path selector (t=64) -- model_out has "
                 f"{n_mo}/8 SNR points so far")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = HERE / "pm_vs_model_out_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")

    fig2, ax2 = plt.subplots(figsize=(10.5, 7))
    if pm_snr:
        s, am, ab, om, ob = zip(*pm_snr)
        ax2.plot(s, om, color="#2a78d6", marker="o", label="pm: oracle")
        ax2.plot(s, am, color="#d6622a", marker="^", label="pm: adaptive")
    s, am, ab, om, ob = zip(*mo_snr)
    ax2.plot(s, om, color="#1fa34a", marker="o", label="model_out: oracle")
    ax2.plot(s, am, color="#a31f8a", marker="^", label="model_out: adaptive")
    ax2.set_xlabel("Eb/N0 (dB)")
    ax2.set_ylabel("mean branches used")
    ax2.set_title("Average pruning width used")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    out2 = HERE / "pm_vs_model_out_meanm.png"
    fig2.savefig(out2, dpi=200)
    print(f"wrote {out2}")


if __name__ == "__main__":
    main()
