#!/usr/bin/env python3
"""Same x=eval checkpoint t, y=BLER, one-file-per-SNR convention as
plot_cross_t_stream.py, but sourced from generalizability.txt (the clean
t in {8,32,64,96}, target-errors=500/1000 trial) instead of
cross_t_stream_results.txt."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"
PM_RESULTS = PM_DIR / "generalizability.txt"
MO_RESULTS = HERE / "generalizability.txt"

COLORS = {"pm": "#d6622a", "model_out": "#a31f8a"}


def load(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [(int(r["eval_t"]), float(r["snr"]), float(r["bler"])) for r in rows]


def main():
    pm_rows = load(PM_RESULTS)
    mo_rows = load(MO_RESULTS)
    print(f"pm: {len(pm_rows)} points, model_out: {len(mo_rows)} points")

    if not mo_rows and not pm_rows:
        print("no data yet")
        return

    snr_list = sorted(set(s for _, s, _ in pm_rows) | set(s for _, s, _ in mo_rows))

    for snr in snr_list:
        fig, ax = plt.subplots(figsize=(9, 6.5))
        for side, rows in (("pm", pm_rows), ("model_out", mo_rows)):
            pts = sorted((t, bler) for t, s, bler in rows if abs(s - snr) < 1e-6)
            if not pts:
                continue
            ts, blers = zip(*pts)
            ax.plot(ts, blers, marker="o", color=COLORS[side], label=f"{side} (frozen t=64 selector)")
        ax.axvline(64, color="gray", linestyle=":", linewidth=1, label="training t (in-distribution)")
        ax.set_yscale("log")
        ax.set_xlabel("eval checkpoint t")
        ax.set_ylabel("BLER (real pickBest-among-top-predicted-m, streaming)")
        ax.set_title(f"Cross-checkpoint generalization (generalizability.txt trial), SNR={snr}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        tag = int(round(snr * 100))
        out = HERE / f"cross_t_bler_gen_snr{tag}.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
