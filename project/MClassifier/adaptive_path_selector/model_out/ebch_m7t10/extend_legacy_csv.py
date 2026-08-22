#!/usr/bin/env python3
"""Extend legacy (old CSV-pipeline) points that currently have fewer than
1000 errors, up to --target-errors 1000 (capped at 1,000,000 samples).
Reuses the EXACT seed formula run_cross_t_eval.py/run_t96_eval.py used
(channel-seed=975001+t*10+i, message-seed=875001+t*10+i, i=SNR_LIST
index), so re-running with a higher --samples/--target-errors reproduces
the same initial sample sequence and genuinely continues past where the
old run stopped, rather than starting over with different noise. Writes
to the SAME old CSV path (overwrite in place), then re-scores via the
same predict()/score() proxy the old plot_cross_t.py used (m_required_basin,
via the frozen t=64 selector), updating cross_t_stream_results.txt."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
RESULTS = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
GAIN_DRILL = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/gain_drill")
MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")

TRACE_MODELS = {
    8: RESULTS / "clean_trace_upto_t8_weights.txt",
    16: RESULTS / "clean_trace_upto_t16_weights.txt",
    32: RESULTS / "clean_trace_upto_t32_weights.txt",
    40: GAIN_DRILL / "t40_weights.txt",
    48: GAIN_DRILL / "t48_weights.txt",
    64: RESULTS / "clean_trace_upto_t64_weights.txt",
    96: GAIN_DRILL / "t96_weights.txt",
    122: GAIN_DRILL / "t122_weights.txt",
}
SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
TARGET_ERRORS = 1000
SAMPLES_CAP = 1_000_000
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]

# Explicit, PER-SIDE whitelist: only rows that are ALREADY legacy (old CSV
# pipeline, old 975001/875001 seed formula) and currently under 1000
# errors. Does NOT include fresh new-tool (985001/885001-seeded) rows --
# those have a different seed sequence, so "continuing" them under the OLD
# formula would just be an unrelated fresh run, silently replacing a real
# result with an old-proxy one. pm and model_out diverged in which points
# are legacy vs fresh (pm raced ahead live further than model_out did), so
# these lists are NOT the same -- verified against each side's own file
# before use.
TARGETS_BY_SIDE = {
    "model_out": [
        (16, 2.75), (32, 2.75), (32, 3.00), (40, 2.25), (40, 2.50), (40, 2.75),
        (40, 3.00), (48, 3.00), (64, 2.00), (64, 3.00), (96, 2.00), (96, 2.25),
        (96, 2.50), (96, 2.75), (96, 3.00), (96, 3.50), (96, 4.00), (96, 4.50),
        (122, 2.00), (122, 3.00),
    ],
    "pm": [
        (48, 3.00), (64, 2.00), (64, 3.00), (96, 2.00), (96, 2.25), (96, 2.50),
        (96, 2.75), (96, 3.00), (96, 3.50), (96, 4.00), (96, 4.50),
        (122, 2.00), (122, 3.00),
    ],
}


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


def current_errors(out_txt, t, snr):
    if not out_txt.exists():
        return None
    with open(out_txt, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return float(r["errors"])
    return None


def upsert_row(out_txt, row):
    existing = []
    if out_txt.exists():
        with open(out_txt, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (str(r["eval_t"]), f"{float(r['snr']):.2f}")
    existing = [r for r in existing if key(r) != key(row)]
    existing.append(row)
    with open(out_txt, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", choices=["pm", "model_out"], required=True)
    args = ap.parse_args()

    if args.side == "model_out":
        sim = BUILD / "mclass_adaptive_m_dataset_model_out"
        selector_pt = MO_DIR / "data" / "model_t64_snr2.pt"
        marker = "_sorted_"
        out_txt = MO_DIR / "cross_t_stream_results.txt"
        csv_dir = MO_DIR / "data" / "cross_t"
    else:
        sim = BUILD / "mclass_adaptive_m_dataset"
        selector_pt = PM_DIR / "data" / "model_t64_snr2.pt"
        marker = "_pm_min_sorted_"
        out_txt = PM_DIR / "cross_t_stream_results.txt"
        csv_dir = PM_DIR / "data" / "cross_t"
    ckpt = torch.load(selector_pt, map_location="cpu", weights_only=False)

    for t, snr in TARGETS_BY_SIDE[args.side]:
            i = SNR_LIST.index(snr)
            cur_err = current_errors(out_txt, t, snr)
            if cur_err is None or cur_err >= TARGET_ERRORS:
                print(f"{args.side} t={t} SNR={snr}: already >=1000 or missing, skipping", flush=True)
                continue
            tag = int(round(snr * 100))
            csv_path = csv_dir / f"crosst_t{t}_snr{tag}.csv"
            cs = 975001 + t * 10 + i
            ms = 875001 + t * 10 + i
            cmd = [str(sim), "-ini", str(INI), "--snr", str(snr), "--samples", str(SAMPLES_CAP),
                   "--target-errors", str(TARGET_ERRORS),
                   "--channel-seed", str(cs), "--message-seed", str(ms), "--out", str(csv_path)]
            if args.side == "model_out":
                cmd += ["--trace-model", str(TRACE_MODELS[t])]
            else:
                cmd += ["--checkpoint", str(t)]
            print(f"=== {args.side} t={t} SNR={snr} (was {cur_err:.0f} errors) ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)

            X, y, any_basin = load_eval(csv_path, marker)
            if y is None:
                print(f"  no rows written for t={t} SNR={snr}, skipping", flush=True)
                continue
            pred_m = predict(ckpt, X)
            mean_m, bler = score(pred_m, y, any_basin)
            n = len(y)
            upsert_row(out_txt, {"eval_t": t, "snr": snr, "samples": n, "errors": round(bler * n),
                                  "mean_m": f"{mean_m:.6f}", "bler": f"{bler:.6f}"})
            print(f"{args.side} t={t} SNR={snr} extended: n={n} mean_m={mean_m:.3f} bler={bler:.6f} "
                  f"[legacy proxy metric]", flush=True)
    print(f"EXTEND_LEGACY_{args.side.upper()}_DONE", flush=True)


if __name__ == "__main__":
    main()
