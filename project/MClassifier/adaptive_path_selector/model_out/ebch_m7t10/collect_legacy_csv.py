#!/usr/bin/env python3
"""One-off salvage: convert already-completed old-style per-sample CSVs
(from the retired mclass_adaptive_m_dataset/_model_out + plot_cross_t.py
pipeline) into the new lightweight cross_t_stream_results.txt format, so
run_cross_t_stream.py's resume logic skips re-simulating checkpoints we
already have real data for. Uses the SAME predict()/score() replay as the
old plot_cross_t.py (m_required_basin-presence proxy, not the new
streaming tool's true pickBest-correctness) -- these converted rows are a
close but not bit-identical approximation to what the streaming tool would
compute; kept only as a resume seed for checkpoints the live streaming
run hasn't reached yet, never overwriting an existing (better, true-metric)
row. Only inserts keys not already present at write time, to stay safe
against a concurrently-running live streaming driver."""
import csv
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
PM_DIR = HERE.parent.parent / "pm" / "ebch_m7t10"

SIDES = {
    "pm": (PM_DIR, PM_DIR / "weights_t64_snr2.txt", "_pm_min_sorted_"),
    "model_out": (HERE, HERE / "weights_t64_snr2p0.txt", "_sorted_"),
}
SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]


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


def already_have(out_txt, t, snr):
    if not out_txt.exists():
        return False
    with open(out_txt, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return True
    return False


def append_row(out_txt, row):
    existing = []
    if out_txt.exists():
        with open(out_txt, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (str(r["eval_t"]), f"{float(r['snr']):.2f}")
    have = {key(r) for r in existing}
    if key(row) not in have:
        existing.append(row)
    with open(out_txt, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    for side_name, (side_dir, model_pt_dir_hint, marker) in SIDES.items():
        # data/model_t64_snr2.pt lives in <side_dir>/data/
        ckpt_path = side_dir / "data" / "model_t64_snr2.pt"
        ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        out_txt = side_dir / "cross_t_stream_results.txt"
        cross_t_dir = side_dir / "data" / "cross_t"

        for t_dir_t in [8, 16, 32, 40, 48, 64, 96, 122]:
            for snr in SNR_LIST:
                tag = int(round(snr * 100))
                csv_path = cross_t_dir / f"crosst_t{t_dir_t}_snr{tag}.csv"
                log_path = cross_t_dir / f"crosst_t{t_dir_t}_snr{tag}_gen.log"
                if not csv_path.exists() or not log_path.exists():
                    continue
                if "wrote " not in log_path.read_text():
                    continue
                if already_have(out_txt, t_dir_t, snr):
                    continue
                X, y, any_basin = load_eval(csv_path, marker)
                if y is None:
                    continue
                pred_m = predict(ckpt, X)
                mean_m, bler = score(pred_m, y, any_basin)
                n = len(y)
                # Re-check right before writing (best-effort race guard
                # against the concurrently-running live streaming driver).
                if already_have(out_txt, t_dir_t, snr):
                    continue
                append_row(out_txt, {"eval_t": t_dir_t, "snr": snr, "samples": n,
                                      "errors": round(bler * n), "mean_m": f"{mean_m:.6f}",
                                      "bler": f"{bler:.6f}"})
                print(f"{side_name} t={t_dir_t} SNR={snr}: salvaged from old CSV "
                      f"(n={n}, mean_m={mean_m:.3f}, bler={bler:.6f}) [legacy proxy metric]", flush=True)
    print("LEGACY_COLLECT_DONE", flush=True)


if __name__ == "__main__":
    main()
