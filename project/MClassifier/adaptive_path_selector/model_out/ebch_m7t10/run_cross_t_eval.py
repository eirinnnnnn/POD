#!/usr/bin/env python3
"""Cross-checkpoint generalization test for model_out: generate eval sets
at t in {8,16,32,40,48,64,122}, over the FULL usual SNR grid (matching
run_sweep_gen.py's schedule -- 6k/6k/6k/6k/10k samples at SNR<=3.00, 1M
each at 3.50/4.00/4.50), using each t's OWN shipped trace_learnt model as
--trace-model. These get scored later by the FROZEN t=64-trained
adaptive-m selector (model_t64_snr2.pt) -- the hypothesis is that since
its input is a calibrated [0,1] probability (not pm's raw,
t-scale-dependent pm_min), one selector trained at t=64 might generalize
across checkpoints without retraining. t=64 itself is included as a
self-consistency sanity check. Resumable: skips any (t,snr) whose file
already has >= the target sample count."""
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
SIM = BUILD / "mclass_adaptive_m_dataset_model_out"
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
RESULTS = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
GAIN_DRILL = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/gain_drill")
HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "data" / "cross_t"

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
SNR_SAMPLES = [
    (2.00, 6_000), (2.25, 6_000), (2.50, 6_000), (2.75, 6_000),
    (3.00, 10_000), (3.50, 1_000_000), (4.00, 1_000_000), (4.50, 1_000_000),
]


def main():
    for t, model_path in TRACE_MODELS.items():
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
                   "--trace-model", str(model_path),
                   "--channel-seed", str(cs), "--message-seed", str(ms),
                   "--out", str(out)]
            print(f"=== model_out t={t} SNR={snr} ===", flush=True)
            with open(log, "w") as lf:
                result = subprocess.run(cmd, cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}, see {log}", flush=True)
                sys.exit(1)
            print(f"wrote {out}", flush=True)
    print("CROSS_T_MODEL_OUT_DONE", flush=True)


if __name__ == "__main__":
    main()
