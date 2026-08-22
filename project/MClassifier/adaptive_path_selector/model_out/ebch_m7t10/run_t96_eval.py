#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
SIM = BUILD / "mclass_adaptive_m_dataset_model_out"
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
MODEL = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/gain_drill/t96_weights.txt")
OUTDIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10/data/cross_t")
SNR_SAMPLES = [(2.00,6000),(2.25,6000),(2.50,6000),(2.75,6000),(3.00,10000),(3.50,1000000),(4.00,1000000),(4.50,1000000)]
t = 96
for i, (snr, samples) in enumerate(SNR_SAMPLES):
    out = OUTDIR / f"crosst_t{t}_snr{int(snr*100)}.csv"
    log = OUTDIR / f"crosst_t{t}_snr{int(snr*100)}_gen.log"
    if out.exists():
        with open(out) as f:
            n = sum(1 for _ in f) - 1
        if n >= samples:
            print(f"t={t} SNR={snr} already have {n} samples, skipping", flush=True)
            continue
    cs = 975001 + t*10 + i
    ms = 875001 + t*10 + i
    cmd = [str(SIM), "-ini", str(INI), "--snr", str(snr), "--samples", str(samples),
           "--target-errors", "100",
           "--trace-model", str(MODEL), "--channel-seed", str(cs), "--message-seed", str(ms),
           "--out", str(out)]
    print(f"=== model_out t=96 SNR={snr} ===", flush=True)
    with open(log, "w") as lf:
        r = subprocess.run(cmd, cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(f"FAILED SNR={snr}, see {log}", flush=True); sys.exit(1)
    print(f"wrote {out}", flush=True)
print("T96_MODEL_OUT_DONE", flush=True)
