#!/usr/bin/env python3
"""Cross-checkpoint generalization test for pm, streaming version:
mclass_adaptive_m_bler_sim computes mean_m/BLER directly (decode -> raw
pm_min ranking -> frozen t=64 adaptive-m selector -> pickBest -> check
correctness), no per-sample CSV. Replaces run_cross_t_eval.py +
run_t96_eval.py + plot_cross_t.py's replay step. Resumable via each
point's own --resume-file; safe to kill and rerun."""
import subprocess
import sys
import csv
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
SIM = BUILD / "mclass_adaptive_m_bler_sim"
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "data" / "cross_t_stream"
SELECTOR = HERE / "weights_t64_snr2.txt"

T_LIST = [8, 16, 32, 40, 48, 64, 96, 122]
SNR_SAMPLES = [
    (2.00, 30_000), (2.25, 50_000), (2.50, 80_000), (2.75, 150_000),
    (3.00, 300_000), (3.50, 700_000), (4.00, 1_500_000), (4.50, 1_000_000),
]
TARGET_ERRORS = 100
OUT_TXT = HERE / "cross_t_stream_results.txt"
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]


def parse_summary(stdout):
    vals = {}
    for line in stdout.splitlines():
        if "," in line and not line.startswith("["):
            k_, v_ = line.rsplit(",", 1)
            try:
                vals[k_] = float(v_)
            except ValueError:
                pass
    return vals


def append_row(row):
    # Live-computed rows always win over anything already present (e.g. a
    # legacy-CSV-salvaged placeholder written by collect_legacy_csv.py
    # while this point was still mid-simulation) -- replace, don't skip.
    existing = []
    if OUT_TXT.exists():
        with open(OUT_TXT, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (str(r["eval_t"]), f"{float(r['snr']):.2f}")
    existing = [r for r in existing if key(r) != key(row)]
    existing.append(row)
    with open(OUT_TXT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def already_done(t, snr):
    if not OUT_TXT.exists():
        return False
    with open(OUT_TXT, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return True
    return False


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for t in T_LIST:
        for i, (snr, cap) in enumerate(SNR_SAMPLES):
            if already_done(t, snr):
                print(f"t={t} SNR={snr} already done, skipping", flush=True)
                continue
            cs = 985001 + t * 10 + i
            ms = 885001 + t * 10 + i
            resume = OUTDIR / f"t{t}_snr{int(snr*100)}.state"
            cmd = [str(SIM), "-ini", str(INI), "--snr", str(snr), "--samples", str(cap),
                   "--target-errors", str(TARGET_ERRORS),
                   "--checkpoint", str(t), "--selector-model", str(SELECTOR),
                   "--channel-seed", str(cs), "--message-seed", str(ms),
                   "--monitor-every", "20000",
                   "--resume-file", str(resume), "--resume-every", "20000"]
            print(f"=== pm t={t} SNR={snr} ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)
            vals = parse_summary(result.stdout)
            append_row({"eval_t": t, "snr": snr, "samples": vals["samples"], "errors": vals["errors"],
                        "mean_m": vals["mean_m"], "bler": vals["bler"]})
            print(f"t={t} SNR={snr} done: samples={vals['samples']:.0f} errors={vals['errors']:.0f} "
                  f"mean_m={vals['mean_m']:.3f} bler={vals['bler']:.6f}", flush=True)
    print("STREAM_PM_DONE", flush=True)


if __name__ == "__main__":
    main()
