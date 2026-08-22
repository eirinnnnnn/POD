#!/usr/bin/env python3
"""Cross-checkpoint generalization test: score each side's FROZEN
t=64-trained adaptive-m selector (data/model_t64_snr2.pt) against eval
sets generated at OTHER checkpoints (data/cross_t/crosst_t{T}_snr{S}.csv,
from run_cross_t_eval.py), without any retraining. Compares model_out
(input = trace-model probability, hypothesized to generalize across t
since it's a bounded, calibrated quantity) against pm (input = raw
pm_min, whose scale/distribution shifts with t, expected to degrade
faster). t=64 itself is the in-distribution sanity-check point."""
import csv
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

T_LIST = [8, 16, 32, 40, 48, 64, 96, 122]
SNR_SAMPLES = {
    2.00: 6_000, 2.25: 6_000, 2.50: 6_000, 2.75: 6_000,
    3.00: 10_000, 3.50: 1_000_000, 4.00: 1_000_000, 4.50: 1_000_000,
}
SNR_LIST = sorted(SNR_SAMPLES)


def load_eval(path, feature_marker):
    rows = list(csv.DictReader(open(path)))
    if rows and any(v is None for v in rows[-1].values()):
        rows = rows[:-1]
    if not rows:
        return None, None, None
    cols = [c for c in rows[0].keys() if feature_marker in c]
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


def run_side(sweep_dir, ckpt_path, feature_marker):
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    rows = []
    for t in T_LIST:
        for snr in SNR_LIST:
            f = sweep_dir / "data" / "cross_t" / f"crosst_t{t}_snr{int(snr*100)}.csv"
            log = sweep_dir / "data" / "cross_t" / f"crosst_t{t}_snr{int(snr*100)}_gen.log"
            if not f.exists():
                print(f"  missing {f}, skipping")
                continue
            # Completion is signaled by the generator's own "wrote ..." line
            # in its log, not by a fixed row-count threshold: with
            # --target-errors early-stopping, a genuinely-finished point can
            # have far fewer rows than the nominal --samples cap (e.g.
            # 74k instead of 1M), so comparing against SNR_SAMPLES would
            # wrongly treat it as still in progress.
            if not (log.exists() and "wrote " in log.read_text()):
                print(f"  t={t} SNR={snr}: no completion marker yet, still in progress, skipping")
                continue
            X, y, any_basin = load_eval(f, feature_marker)
            if y is None:
                continue
            pred_m = predict(ckpt, X)
            mean_m, bler = score(pred_m, y, any_basin)
            rows.append((t, snr, mean_m, bler))
    return rows


def main():
    print("pm cross-t (frozen t=64 selector):")
    pm_rows = run_side(PM_DIR, PM_MODEL_PT, "_pm_min_sorted_")
    print("model_out cross-t (frozen t=64 selector):")
    mo_rows = run_side(HERE, MO_MODEL_PT, "_sorted_")

    out_txt = HERE / "cross_t_results.txt"
    with open(out_txt, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["side", "eval_t", "snr", "mean_m", "bler"])
        for side, rows in (("pm", pm_rows), ("model_out", mo_rows)):
            for t, snr, mean_m, bler in rows:
                w.writerow([side, t, snr, f"{mean_m:.6f}", f"{bler:.6f}"])
    print(f"wrote {out_txt}")

    # vs-t view (one line per side, x=eval checkpoint): one figure PER SNR
    # point actually available, named cross_t_bler_snr{tag}.png (tag =
    # int(snr*100), matching the sweep_snr{tag}.csv convention used
    # elsewhere), so every SNR gets its own file instead of cramming a
    # fixed subset into shared panels.
    colors = {"pm": "#d6622a", "model_out": "#a31f8a"}
    for snr in SNR_LIST:
        has_data = any(s == snr for _, s, _, _ in pm_rows) or any(s == snr for _, s, _, _ in mo_rows)
        if not has_data:
            continue
        fig, ax = plt.subplots(figsize=(9, 6.5))
        for side, rows in (("pm", pm_rows), ("model_out", mo_rows)):
            pts = sorted((t, bler) for t, s, m, bler in rows if s == snr)
            if not pts:
                continue
            ts, blers = zip(*pts)
            ax.plot(ts, blers, marker="o", color=colors[side], label=f"{side} (frozen t=64 selector)")
        ax.axvline(64, color="gray", linestyle=":", linewidth=1, label="training t (in-distribution)")
        ax.set_yscale("log")
        ax.set_xlabel("eval checkpoint t")
        ax.set_ylabel("BLER (using frozen selector's predicted m)")
        ax.set_title(f"Cross-checkpoint generalization, SNR={snr} -- t=64-trained selector evaluated at other t")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        tag = int(round(snr * 100))
        out_png = HERE / f"cross_t_bler_snr{tag}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)
        print(f"wrote {out_png}")

    # "usual" BLER-vs-Eb/N0 plot: one line per eval checkpoint t, two
    # panels (pm | model_out), matching the convention used everywhere
    # else in this project (e.g. t_sweep_m8_bler.png).
    fig3, axes3 = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
    t_colors = dict(zip(T_LIST, plt.get_cmap("viridis")(np.linspace(0.1, 0.9, len(T_LIST)))))
    for ax, (side, rows) in zip(axes3, (("pm", pm_rows), ("model_out", mo_rows))):
        for t in T_LIST:
            pts = sorted((snr, bler) for tt, snr, m, bler in rows if tt == t)
            if not pts:
                continue
            snrs, blers = zip(*pts)
            style = "-" if t == 64 else "--"
            ax.plot(snrs, blers, linestyle=style, marker="o", markersize=5,
                     color=t_colors[t], label=f"eval t={t}" + (" (train)" if t == 64 else ""))
        ax.set_yscale("log")
        ax.set_xlabel("Eb/N0 (dB)")
        ax.set_title(f"{side}: frozen t=64 selector")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="lower left")
    axes3[0].set_ylabel("BLER (using frozen selector's predicted m)")
    fig3.suptitle("Cross-checkpoint generalization -- BLER vs Eb/N0, frozen t=64 selector evaluated at each t")
    fig3.tight_layout()
    out_png3 = HERE / "cross_t_bler_vs_snr.png"
    fig3.savefig(out_png3, dpi=200)
    print(f"wrote {out_png3}")

    fig2, axes2 = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
    for ax, snr in zip(axes2, VS_T_SNRS):
        for side, rows in (("pm", pm_rows), ("model_out", mo_rows)):
            pts = sorted((t, m) for t, s, m, bler in rows if s == snr)
            if not pts:
                continue
            ts, ms = zip(*pts)
            ax.plot(ts, ms, marker="o", color=colors[side], label=f"{side} (frozen t=64 selector)")
        ax.axvline(64, color="gray", linestyle=":", linewidth=1, label="training t (in-distribution)")
        ax.set_xlabel("eval checkpoint t")
        ax.set_title(f"SNR={snr}")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    axes2[0].set_ylabel("mean branches used")
    fig2.suptitle("Cross-checkpoint generalization: mean pruning width used")
    fig2.tight_layout()
    out_png2 = HERE / "cross_t_meanm.png"
    fig2.savefig(out_png2, dpi=200)
    print(f"wrote {out_png2}")


if __name__ == "__main__":
    main()
