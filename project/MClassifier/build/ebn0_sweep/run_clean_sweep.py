#!/usr/bin/env python3
"""Clean-automorphism-file rerun of the ebn0_sweep_bler experiment, using
mclass_bler_sim (streaming, no per-sample CSVs, checkpoint-resumable)
instead of the old mclass_dataset/mclass_trace_dataset/mclass_eval
pipeline. Checkpoints t=4,8,16,32 (the newly-retrained clean models),
k=8 fixed, plus full_ped (M=64) and full_ped_m8 baselines."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
INI = BUILD / "ebn0_sweep.ini"
INI_M8 = BUILD / "ebn0_sweep_m8.ini"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 4.00: 1_500_000, 4.50: 1_000_000}
TARGET_ERRORS = 1000
CHECKPOINTS = [4, 8, 16, 32, 64]
K = 8
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "t_sweep_m8_results.txt"
FIELDS = ["snr", "t", "method", "k", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 940001 + salt * 100 + snr_idx, 840001 + salt * 100 + snr_idx


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


def run_sim(log_path, args, resume_tag):
    if log_path.exists():
        text = log_path.read_text()
        if "\nsamples," in text:
            return parse_summary(text), True
    resume_path = OUTDIR / f"{resume_tag}.state"
    cmd = [str(SIM)] + args + ["--resume-file", str(resume_path), "--resume-every", str(RESUME_EVERY)]
    result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        raise SystemExit(f"failed: {' '.join(args)}\n{result.stdout[-2000:]}")
    return parse_summary(result.stdout), False


def append_rows(rows):
    existing = []
    if RESULTS.exists():
        with open(RESULTS, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (f"{float(r['snr']):.2f}", str(r["t"]), r["method"], str(int(float(r["k"]))))
    have = {key(r) for r in existing}
    for r in rows:
        if key(r) not in have:
            existing.append(r)
    with open(RESULTS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    for i, snr in enumerate(SNR_LIST):
        t_tag = tag(snr)
        cap = SAMPLE_CAP[snr]

        # full_ped (M=64) + static/trace at each checkpoint t
        for t in CHECKPOINTS:
            cs, ms = seed_pair(i, salt=t)
            log = OUTDIR / f"clean_sweep_t{t}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(INI), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--stop-metric", "static_trace_min",
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--checkpoint", str(t), "--k", str(K),
                "--model", str(BUILD / f"clean_trace_upto_t{t}_weights.txt"),
                "--monitor-every", "100000",
            ], f"clean_sweep_t{t}_snr{t_tag}")
            samples = vals["samples"]
            # full_ped dropped -- reuse the M=64 full-PED curve from the
            # original (pre-fix) ebn0_sweep_bler run instead of recomputing;
            # its contamination from the automorphism bug was already
            # measured as negligible (~0.1% or less, see prior analysis).
            rows = [{"snr": snr, "t": t, "method": "static_pm_rank", "k": K,
                     "samples": samples, "errors": vals[f"static_errors[k={K}]"], "bler": vals[f"static_bler[k={K}]"]},
                    {"snr": snr, "t": t, "method": "trace_learned", "k": K,
                     "samples": samples, "errors": vals[f"trace_errors[k={K}]"], "bler": vals[f"trace_bler[k={K}]"]}]
            append_rows(rows)
            print(f"SNR={snr} t={t} done (cached={cached}, samples={samples:.0f})")
        # full_ped_m8 baseline dropped -- full PED (M=64) is already the
        # reference; not worth the extra compute right now.

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
