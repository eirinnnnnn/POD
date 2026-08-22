#!/usr/bin/env python3
"""Static-pm-rank-only (no model) BLER over the usual SNR grid, for two
5-checkpoint windows t in {14..18} and t in {22..26}, centered on the
frozen bits at decode_idx=16 (span 4) and decode_idx=24 (span 4)
respectively. Checkpoint t means positions 0..t-1 are decided, so the
frozen bit at idx=F first enters the metric at t=F+1 -- e.g. for the
idx=16 window, t=16 is "before" and t=17 is "after"."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
RESULTS_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
INI = RESULTS_DIR / "ebn0_sweep.ini"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50]  # capped at 3.5 -- the
# idx=16/idx=24 before/after comparison is already conclusive by 3.5dB,
# not worth the 700k-1.5M-sample cost of 4.0/4.5 for these two windows.
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000}
TARGET_ERRORS = 1000
CHECKPOINTS = [14, 15, 16, 17, 18, 22, 23, 24, 25, 26]
K = 8
RESUME_EVERY = 20_000

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "frozen_window_static_results.txt"
FIELDS = ["snr", "t", "method", "k", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 950001 + salt * 100 + snr_idx, 850001 + salt * 100 + snr_idx


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
        for t in CHECKPOINTS:
            cs, ms = seed_pair(i, t)
            log = OUTDIR / f"frozen_window_t{t}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(INI), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--checkpoint", str(t), "--k-list", str(K),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"frozen_window_t{t}_snr{t_tag}")
            samples = vals["samples"]
            rows = [{"snr": snr, "t": t, "method": "static_pm_rank", "k": K,
                     "samples": samples, "errors": vals[f"static_errors[k={K}]"], "bler": vals[f"static_bler[k={K}]"]}]
            append_rows(rows)
            print(f"SNR={snr} t={t} done (cached={cached}, samples={samples:.0f})", flush=True)

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
