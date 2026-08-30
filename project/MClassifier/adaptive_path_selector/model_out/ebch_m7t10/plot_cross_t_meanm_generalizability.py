#!/usr/bin/env python3
"""x=eval checkpoint t, y=mean_m, one file per SNR present in
generalizability.txt, both sides overlaid -- mean_m companion to
plot_cross_t_generalizability.py's BLER plots. Also overlays, wherever a
matching (side, t, snr) 12k-sample training CSV exists, the oracle mean
m_required_basin from that training data -- the true minimal list width
needed, vs. mean_m which is what the selector actually predicted/used."""
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"
MO_DIR = HERE
PM_RESULTS = PM_DIR / "generalizability.txt"
MO_RESULTS = MO_DIR / "generalizability.txt"

COLORS = {"pm": "#d6622a", "model_out": "#a31f8a"}

SNR_TAG_MAP = {"snr2p0": 2.0, "snr2": 2.0, "snr3": 3.0, "snr4": 4.0}


def load(path):
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [(int(r["eval_t"]), float(r["snr"]), float(r["mean_m"])) for r in rows]


def load_oracle(data_dir):
    # returns {(t, snr): oracle_mean_m} from any train_t{t}_{snrtag}_12k.csv
    # with an m_required_basin column present
    out = {}
    for f in sorted(data_dir.glob("train_t*_12k.csv")):
        m = re.match(r"train_t(\d+)_(snr2p0|snr2|snr3|snr4)_12k\.csv", f.name)
        if not m:
            continue
        t = int(m.group(1))
        snr = SNR_TAG_MAP[m.group(2)]
        with open(f, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            if "m_required_basin" not in (reader.fieldnames or []):
                continue
            vals = []
            for row in reader:
                v = row.get("m_required_basin")
                if v is None or v == "":
                    continue
                try:
                    vals.append(float(v))
                except ValueError:
                    continue
        if not vals:
            continue
        key = (t, snr)
        # prefer the larger sample (more reliable estimate) if duplicate tags collide
        if key not in out or len(vals) > out[key][1]:
            out[key] = (sum(vals) / len(vals), len(vals))
    return {k: v[0] for k, v in out.items()}


def main():
    pm_rows = load(PM_RESULTS)
    mo_rows = load(MO_RESULTS)
    pm_oracle = load_oracle(PM_DIR / "data")
    mo_oracle = load_oracle(MO_DIR / "data")
    print(f"pm: {len(pm_rows)} points, model_out: {len(mo_rows)} points")
    print(f"pm oracle points: {pm_oracle}")
    print(f"model_out oracle points: {mo_oracle}")

    if not mo_rows and not pm_rows:
        print("no data yet")
        return

    snr_list = sorted(set(s for _, s, _ in pm_rows) | set(s for _, s, _ in mo_rows))

    for snr in snr_list:
        fig, ax = plt.subplots(figsize=(9, 6.5))
        for side, rows, oracle in (("pm", pm_rows, pm_oracle), ("model_out", mo_rows, mo_oracle)):
            pts = sorted((t, m) for t, s, m in rows if abs(s - snr) < 1e-6)
            if pts:
                ts, means = zip(*pts)
                ax.plot(ts, means, marker="o", color=COLORS[side], label=f"{side} used mean_m")
            opts = sorted((t, v) for (t, s), v in oracle.items() if abs(s - snr) < 1e-6)
            if opts:
                ots, ovals = zip(*opts)
                ax.scatter(ots, ovals, marker="*", s=180, color=COLORS[side],
                           edgecolors="black", linewidths=0.8, zorder=5,
                           label=f"{side} oracle m_required_basin (training data)")
        ax.axvline(64, color="gray", linestyle=":", linewidth=1, label="training t (in-distribution)")
        ax.set_xlabel("eval checkpoint t")
        ax.set_ylabel("mean m (used vs. oracle-required)")
        ax.set_title(f"Cross-checkpoint mean_m vs oracle (generalizability.txt trial), SNR={snr}")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7.5)
        fig.tight_layout()
        tag = int(round(snr * 100))
        out = HERE / f"cross_t_meanm_snr{tag}.png"
        fig.savefig(out, dpi=200)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
