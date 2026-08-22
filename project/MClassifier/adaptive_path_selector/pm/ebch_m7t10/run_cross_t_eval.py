#!/usr/bin/env python3
"""Cross-checkpoint generalization test for pm: generate eval sets at t in
{8,16,32,40,48,64,122}, over the FULL usual SNR grid (matching
run_sweep_gen.py's schedule -- 6k/6k/6k/6k/10k samples at SNR<=3.00, 1M
each at 3.50/4.00/4.50), using pm's own mclass_adaptive_m_dataset (raw
sorted checkpoint pm_min, no trace model). These get scored later by the
FROZEN t=64-trained adaptive-m selector (data/model_t64_snr2.pt) --
pm_min is an unbounded, monotonically-growing path metric whose
scale/distribution shifts a lot with t, so the hypothesis is this frozen
selector should degrade much faster across checkpoints than model_out's
does. t=64 itself is included as a self-consistency sanity check.
Resumable: skips any (t,snr) whose file already has >= the target sample
count."""
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
SIM = BUILD / "mclass_adaptive_m_dataset"
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "data" / "cross_t"

T_LIST = [8, 16, 32, 40, 48, 64, 96, 122]
SNR_SAMPLES = [
    (2.00, 6_000), (2.25, 6_000), (2.50, 6_000), (2.75, 6_000),
    (3.00, 10_000), (3.50, 1_000_000), (4.00, 1_000_000), (4.50, 1_000_000),
]


def main():
    for t in T_LIST:
        for i, (snr, samples) in enumerate(SNR_SAMPLES):
            out = OUTDIR / f"crosst_t{t}_snr{int(snr*100)}.csv"
            log = OUTDIR / f"crosst_t{t}_snr{int(snr*100)}_gen.log"
            if out.exists():
                with open(out) as f:
                    n = sum(1 for _ in f) - 1
                if n >= samples:
                    print(f"t={t} SNR={snr} already have {n} samples, skipping", flush=True)
                    continue
            cs = 975001 + t * 10 + i
            ms = 875001 + t * 10 + i
            cmd = [str(SIM), "-ini", str(INI), "--snr", str(snr), "--samples", str(samples),
                   "--target-errors", "100",
                   "--checkpoint", str(t),
                   "--channel-seed", str(cs), "--message-seed", str(ms),
                   "--out", str(out)]
            print(f"=== pm t={t} SNR={snr} ===", flush=True)
            with open(log, "w") as lf:
                result = subprocess.run(cmd, cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}, see {log}", flush=True)
                sys.exit(1)
            print(f"wrote {out}", flush=True)
    print("CROSS_T_PM_DONE", flush=True)


if __name__ == "__main__":
    main()
