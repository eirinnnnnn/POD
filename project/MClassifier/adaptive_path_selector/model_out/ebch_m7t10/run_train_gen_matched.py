#!/usr/bin/env python3
"""Generate matched-t training sets (12k samples, SNR=2.0) for the
mismatched-vs-matched APS comparison: for t in {8,32,96}, on both sides,
mirroring the exact schedule/format used to build train_t64_snr2p0_12k.csv
(single SNR=2.0, 12000 samples, m_required_basin/any_basin_exists label
columns already emitted by the current dataset-tool build)."""
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
RESULTS = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")

TRACE_MODELS = {
    8: RESULTS / "clean_trace_upto_t8_weights.txt",
    32: RESULTS / "clean_trace_upto_t32_weights.txt",
    96: Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/gain_drill/t96_weights.txt"),
}
T_LIST = [8, 32, 96]
SNR = 2.0
SAMPLES = 12000


def has_basin_col(path):
    if not path.exists():
        return False
    with open(path) as f:
        header = f.readline()
    return "m_required_basin" in header.split(",")


def complete(path, n_needed):
    if not path.exists():
        return False
    if not has_basin_col(path):
        return False
    with open(path) as f:
        lines = sum(1 for _ in f)
    return lines - 1 >= n_needed


def main():
    for t in T_LIST:
        # model_out
        sim = BUILD / "mclass_adaptive_m_dataset_model_out"
        out = MO_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        cs = 945001 + t
        ms = 845001 + t
        if complete(out, SAMPLES):
            print(f"model_out t={t}: already complete, skipping", flush=True)
        else:
            cmd = [str(sim), "-ini", str(INI), "--snr", str(SNR), "--samples", str(SAMPLES),
                   "--trace-model", str(TRACE_MODELS[t]),
                   "--channel-seed", str(cs), "--message-seed", str(ms), "--out", str(out)]
            print(f"=== model_out t={t} ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED model_out t={t}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)
            print(f"wrote {out}", flush=True)

        # pm
        sim = BUILD / "mclass_adaptive_m_dataset"
        out = PM_DIR / "data" / f"train_t{t}_snr2p0_12k.csv"
        if complete(out, SAMPLES):
            print(f"pm t={t}: already complete, skipping", flush=True)
        else:
            cmd = [str(sim), "-ini", str(INI), "--snr", str(SNR), "--samples", str(SAMPLES),
                   "--checkpoint", str(t),
                   "--channel-seed", str(cs), "--message-seed", str(ms), "--out", str(out)]
            print(f"=== pm t={t} ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED pm t={t}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)
            print(f"wrote {out}", flush=True)
    print("TRAIN_GEN_MATCHED_DONE", flush=True)


if __name__ == "__main__":
    main()
