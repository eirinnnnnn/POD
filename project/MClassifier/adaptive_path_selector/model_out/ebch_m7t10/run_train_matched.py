#!/usr/bin/env python3
"""Train matched-t APS selectors (t in {8,32,96}) on both sides, using the
same train_adaptive_m.py hyperparameters/format as the existing t=64
selectors (weights_t64_snr2p0.txt) -- default hidden=64/epochs=40/lr=1e-3,
model_out with --no-log-input (trace_prob is already bounded [0,1]), pm
with the default log1p(pm_min) transform."""
import subprocess
import sys
from pathlib import Path

MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")
T_LIST = [8, 32, 96]


def run(cmd, tag):
    print(f"=== {tag} ===", flush=True)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout[-3000:], flush=True)
    if result.returncode != 0:
        print(f"FAILED {tag}", flush=True)
        sys.exit(1)


def main():
    for t in T_LIST:
        data = MO_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        out = MO_DIR / f"weights_t{t}_snr2p0.txt"
        cmd = ["python3", str(MO_DIR / "train_adaptive_m.py"), "--data", str(data), "--out", str(out),
               "--no-log-input"]
        run(cmd, f"model_out t={t}")

        data = PM_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        out = PM_DIR / f"weights_t{t}_snr2p0.txt"
        cmd = ["python3", str(PM_DIR / "train_adaptive_m.py"), "--data", str(data), "--out", str(out)]
        run(cmd, f"pm t={t}")
    print("TRAIN_MATCHED_DONE", flush=True)


if __name__ == "__main__":
    main()
