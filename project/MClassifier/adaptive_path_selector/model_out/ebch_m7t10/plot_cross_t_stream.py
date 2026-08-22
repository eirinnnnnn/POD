#!/usr/bin/env python3
"""Plot the streaming cross-t generalization results (cross_t_stream_results.txt
from both sides' run_cross_t_stream.py) -- no CSV replay, just read the
small aggregate results files directly. One figure per SNR point, x=eval
checkpoint t, y=BLER, both sides overlaid -- matching the established
cross_t_bler_snr{tag}.png convention. Nothing else."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"
PM_RESULTS = PM_DIR / "cross_t_stream_results.txt"
MO_RESULTS = HERE / "cross_t_stream_results.txt"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
COLORS = {"pm": "#d6622a", "model_out": "#a31f8a"}


def load(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [(int(r["eval_t"]), float(r["snr"]), float(r["mean_m"]), float(r["bler"])) for r in rows]


def main():
    pm_rows = load(PM_RESULTS)
    mo_rows = load(MO_RESULTS)
    print(f"pm: {len(pm_rows)} points, model_out: {len(mo_rows)} points")

    if not mo_rows and not pm_rows:
        print("no data yet")
        return

    for snr in SNR_LIST:
        has_data = any(abs(s - snr) < 1e-6 for _, s, _, _ in pm_rows) or \
                   any(abs(s - snr) < 1e-6 for _, s, _, _ in mo_rows)
        if not has_data:
            continue
        fig, ax = plt.subplots(figsize=(9, 6.5))
        for side, rows in (("pm", pm_rows), ("model_out", mo_rows)):
            pts = sorted((t, bler) for t, s, m, bler in rows if abs(s - snr) < 1e-6)
            if not pts:
                continue
            ts, blers = zip(*pts)
            ax.plot(ts, blers, marker="o", color=COLORS[side], label=f"{side} (frozen t=64 selector)")
        ax.axvline(64, color="gray", linestyle=":", linewidth=1, label="training t (in-distribution)")
        ax.set_yscale("log")
        ax.set_xlabel("eval checkpoint t")
        ax.set_ylabel("BLER (real pickBest-among-top-predicted-m, streaming)")
        ax.set_title(f"Cross-checkpoint generalization (streaming), SNR={snr}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        tag = int(round(snr * 100))
        out = HERE / f"cross_t_bler_snr{tag}.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
