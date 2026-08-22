#!/usr/bin/env python3
"""Generate model_out's sweep_snr*.csv eval sets, mirroring
adaptive_path_selector/pm/ebch_m7t10's sample-count schedule exactly
(6k/6k/6k/6k/10k at SNR<=3.00, 1M each at 3.50/4.00/4.50 -- pm needed that
many at high SNR to get a stable m_required_basin/any_basin_exists
estimate once basin-misses get rare). Uses the shipped t=64 trace model as
--trace-model, matching pm's checkpoint=64. No resume support (matches
pm's own generator, which doesn't have one either) -- if interrupted,
rerun; already-written files are skipped."""
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
SIM = BUILD / "mclass_adaptive_m_dataset_model_out"
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
TRACE_MODEL = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/clean_trace_upto_t64_weights.txt")
HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "data"

SNR_SAMPLES = [
    (2.00, 6_000), (2.25, 6_000), (2.50, 6_000), (2.75, 6_000),
    (3.00, 10_000), (3.50, 1_000_000), (4.00, 1_000_000), (4.50, 1_000_000),
]


def tag(snr):
    digits = f"{snr:.2f}".replace(".", "")
    return digits


def main():
    for i, (snr, n) in enumerate(SNR_SAMPLES):
        out = OUTDIR / f"sweep_snr{tag(snr)}.csv"
        log = OUTDIR / f"sweep_snr{tag(snr)}_gen.log"
        if out.exists():
            with open(out) as f:
                lines = sum(1 for _ in f)
            if lines - 1 >= n:
                print(f"SNR={snr} already have {lines-1} samples, skipping", flush=True)
                continue
        cs = 965001 + i
        ms = 865001 + i
        cmd = [str(SIM), "-ini", str(INI), "--snr", str(snr), "--samples", str(n),
               "--trace-model", str(TRACE_MODEL),
               "--channel-seed", str(cs), "--message-seed", str(ms),
               "--out", str(out)]
        print(f"=== SNR={snr} n={n} ===", flush=True)
        with open(log, "w") as lf:
            result = subprocess.run(cmd, cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            print(f"FAILED at SNR={snr}, see {log}", flush=True)
            sys.exit(1)
        print(f"wrote {out}", flush=True)
    print("ALL_SWEEP_DONE", flush=True)


if __name__ == "__main__":
    main()
