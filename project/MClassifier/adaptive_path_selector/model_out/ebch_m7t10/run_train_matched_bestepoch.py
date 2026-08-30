#!/usr/bin/env python3
"""Retrain t in {8,32,64,96} on both sides with the new best-epoch model
selection added to train_adaptive_m.py, for comparison against the
existing final-epoch weights. Writes to weights_t{T}_snr2p0_bestepoch.txt
(does NOT touch the production weights_t{T}_snr2p0.txt files, which are
in active use by running background trials)."""
import subprocess
import sys
from pathlib import Path

MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")
T_LIST = [8, 32, 64, 96]


def run(cmd, tag, logf):
    print(f"=== {tag} ===", flush=True)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(logf, "a") as f:
        f.write(f"=== {tag} ===\n{result.stdout}\n")
    print(result.stdout[-1500:], flush=True)
    if result.returncode != 0:
        print(f"FAILED {tag}", flush=True)
        sys.exit(1)


def main():
    logf = "/tmp/train_matched_bestepoch.log"
    for t in T_LIST:
        data = MO_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        out = MO_DIR / f"weights_t{t}_snr2p0_bestepoch.txt"
        cmd = ["python3", str(MO_DIR / "train_adaptive_m.py"), "--data", str(data), "--out", str(out),
               "--no-log-input"]
        run(cmd, f"model_out t={t}", logf)

        data = PM_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        out = PM_DIR / f"weights_t{t}_snr2p0_bestepoch.txt"
        cmd = ["python3", str(PM_DIR / "train_adaptive_m.py"), "--data", str(data), "--out", str(out)]
        run(cmd, f"pm t={t}", logf)
    print("TRAIN_MATCHED_BESTEPOCH_DONE", flush=True)


if __name__ == "__main__":
    main()
