#!/usr/bin/env python3
"""BLER-vs-Eb/N0 plot for the adaptive path selector: snr2/t64/hidden32
model vs. fixed static-pm-rank(m=FIXED_M) vs. the oracle (true
m_required_basin), over whichever sweep_snr*.csv eval files currently
exist. Re-runnable at any point mid-sweep -- only plots what's landed so
far. FIXED_M is not a hardcoded constant: it's set to the adaptive model's
OWN mean_m, averaged across the SNR sweep -- so the fixed-m baseline uses
exactly the same average compute budget as the adaptive model, making the
BLER comparison an apples-to-apples "same average m, does adapting
per-sample help" test rather than an arbitrary fixed choice."""
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
MODEL_PT = HERE / "data" / "model_t64_snr2.pt"


def load_eval(path):
    rows = list(csv.DictReader(open(path)))
    if rows and any(v is None for v in rows[-1].values()):
        rows = rows[:-1]  # file still being written -- drop a possibly-truncated trailing row
    cols = [c for c in rows[0].keys() if "_pm_min_sorted_" in c]
    cols.sort(key=lambda c: int(c.rsplit("_", 1)[1]))
    X = np.array([[float(r[c]) for c in cols] for r in rows], dtype=np.float32)
    y = np.array([int(r["m_required_basin"]) for r in rows], dtype=np.float64)
    any_basin = np.array([int(r["any_basin_exists"]) for r in rows], dtype=np.float64)
    return X, y, any_basin


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


def main():
    files = sorted(HERE.glob("data/sweep_snr*.csv"))
    if not files:
        print("no sweep_snr*.csv files yet")
        return
    ckpt = torch.load(MODEL_PT, map_location="cpu", weights_only=False)

    # pass 1: load each SNR point's data and the adaptive model's own
    # predictions, so FIXED_M can be set to their cross-SNR average before
    # scoring the fixed-m baseline in pass 2.
    per_snr = []
    for f in files:
        m = re.search(r"sweep_snr(\d+)\.csv", f.name)
        digits = m.group(1)
        snr = float(digits[0] + "." + digits[1:])
        X, y, any_basin = load_eval(f)
        pred_m = predict(ckpt, X)
        am, ab = score(pred_m, y, any_basin)
        per_snr.append((snr, y, any_basin, am, ab))

    FIXED_M = round(float(np.mean([am for _, _, _, am, _ in per_snr])))
    print(f"FIXED_M = adaptive model's own mean_m averaged across the sweep, rounded = {FIXED_M}")

    # fixed-m reference curves: the adaptive model's own average, plus two
    # oracle-informed fixed widths -- m=1 (the oracle's MAJORITY/mode value
    # at every SNR checked, 72-100% of samples) and m=7 (the oracle's MEAN
    # at the hardest SNR=2.00 point, pulled up by the long right tail of
    # hard samples). These aren't oracle curves themselves (no per-sample
    # adaptation) -- they're fixed-width baselines whose width happens to
    # be informed by oracle statistics, kept here for direct comparison.
    fixed_refs = [
        (f"adaptive path avg (m={FIXED_M})", FIXED_M, "#555555", "-.", "s"),
        ("oracle majority (m=1)", 1, "#8a8a2a", ":", "D"),
        ("oracle avg (m=7)", 7, "#8a2a6a", ":", "P"),
    ]

    # pass 2: score every fixed-m baseline and the oracle
    snrs, adaptive_bler, adaptive_m = [], [], []
    oracle_bler, oracle_m = [], []
    fixed_results = {label: {"bler": [], "m": []} for label, _, _, _, _ in fixed_refs}
    M = ckpt["M"]
    for snr, y, any_basin, am, ab in per_snr:
        oracle_used = np.where(y == 0, M, y)
        om, ob = score(oracle_used, y, any_basin)

        snrs.append(snr)
        adaptive_bler.append(ab); adaptive_m.append(am)
        oracle_bler.append(ob); oracle_m.append(om)

        line = f"SNR={snr:.2f}  adaptive(m={am:.2f})={ab:.4f}  oracle(m={om:.2f})={ob:.4f}"
        for label, mval, *_ in fixed_refs:
            fm, fb = score(np.full_like(y, mval), y, any_basin)
            fixed_results[label]["bler"].append(fb)
            fixed_results[label]["m"].append(fm)
            line += f"  {label}={fb:.4f}"
        print(line)

    order = np.argsort(snrs)
    snrs = np.array(snrs)[order]
    adaptive_bler = np.array(adaptive_bler)[order]
    oracle_bler = np.array(oracle_bler)[order]
    adaptive_m = np.array(adaptive_m)[order]
    oracle_m = np.array(oracle_m)[order]
    for label in fixed_results:
        fixed_results[label]["bler"] = np.array(fixed_results[label]["bler"])[order]
        fixed_results[label]["m"] = np.array(fixed_results[label]["m"])[order]

    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    ax.plot(snrs, oracle_bler, color="#2a78d6", linestyle="-", marker="o", markersize=6,
            linewidth=2, label="oracle (ground truth required m)", zorder=3)
    ax.plot(snrs, adaptive_bler, color="#d6622a", linestyle="--", marker="^", markersize=6,
            linewidth=1.6, label="adaptive path selector", zorder=2)
    for label, mval, color, ls, marker in fixed_refs:
        ax.plot(snrs, fixed_results[label]["bler"], color=color, linestyle=ls, marker=marker,
                 markersize=6, linewidth=1.4, label=label, zorder=2)
    ax.set_yscale("log")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("BLER")
    ax.set_title("Adaptive path selector (t=64): oracle vs adaptive vs fixed-m")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = HERE / "adaptive_sweep_bler.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")

    fig2, ax2 = plt.subplots(figsize=(9.5, 6.5))
    ax2.plot(snrs, oracle_m, color="#2a78d6", marker="o", label="oracle (ground truth required m)")
    ax2.plot(snrs, adaptive_m, color="#d6622a", marker="^", label="adaptive path selector")
    for label, mval, color, ls, marker in fixed_refs:
        ax2.axhline(mval, color=color, linestyle=ls, label=label)
    ax2.set_xlabel("Eb/N0 (dB)")
    ax2.set_ylabel("mean branches used")
    ax2.set_title("Average pruning width used")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    out2 = HERE / "adaptive_sweep_meanm.png"
    fig2.savefig(out2, dpi=200)
    print(f"wrote {out2}")


if __name__ == "__main__":
    main()
