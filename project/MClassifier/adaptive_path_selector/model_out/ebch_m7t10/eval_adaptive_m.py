#!/usr/bin/env python3
"""Evaluate trained adaptive-m models (and the oracle) on held-out eval sets."""
import csv
import sys
import numpy as np
import torch
import torch.nn as nn


def load_eval(path):
    rows = list(csv.DictReader(open(path)))
    cols = [c for c in rows[0].keys() if "_sorted_" in c]
    cols.sort(key=lambda c: int(c.rsplit("_", 1)[1]))
    M = len(cols)
    X = np.array([[float(r[c]) for c in cols] for r in rows], dtype=np.float32)
    y = np.array([int(r["m_required_basin"]) for r in rows], dtype=np.float64)
    any_basin = np.array([int(r["any_basin_exists"]) for r in rows], dtype=np.float64)
    return X, y, any_basin, M


def predict(ckpt_path, X):
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
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
    M = ckpt["M"]
    return np.clip(np.ceil(pred), 1, M)


def score(used_m, y, any_basin, M):
    fail = (any_basin == 0) | (used_m < y)
    return used_m.mean(), fail.mean()


def main():
    eval_snrs = [2, 3, 4]
    model_tags = ["snr2", "snr3", "snr4", "mixed"]
    print(f"{'eval_snr':<10}{'model':<10}{'mean_m':<10}{'BLER':<10}")
    for esnr in eval_snrs:
        X, y, any_basin, M = load_eval(f"data/eval_t64_snr{esnr}_6k.csv")
        for tag in model_tags:
            used_m = predict(f"data/model_t64_{tag}.pt", X)
            mean_m, bler = score(used_m, y, any_basin, M)
            print(f"{esnr:<10}{tag:<10}{mean_m:<10.3f}{bler:<10.4f}")
        # oracle: use the TRUE m_required_basin directly (0 -> use all M, unrecoverable anyway)
        oracle_m = np.where(y == 0, M, y)
        mean_m, bler = score(oracle_m, y, any_basin, M)
        print(f"{esnr:<10}{'oracle':<10}{mean_m:<10.3f}{bler:<10.4f}")
        print()


if __name__ == "__main__":
    main()
