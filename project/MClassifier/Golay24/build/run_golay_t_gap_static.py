#!/usr/bin/env python3
"""Golay24 static-pm-rank probe across the tail gap t in {21,22,23}, filling
in between the t=20 sweep point and the t=24 sanity check (where ranking by
the checkpoint collapses onto full_ped since t=24 IS the final decode
position, n=24). No --model is passed (static_pm_rank is computed
unconditionally by the sim). Rows are appended to BOTH the _single and
_cumulative results files, since static_pm_rank at a given checkpoint
doesn't depend on which model file (if any) happens to be loaded
alongside it."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
INI = BUILD / "golay_sweep_m16_L2.ini"
OUTDIR = BUILD / "sweep"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 45_000, 2.50: 55_000, 2.75: 85_000,
              3.00: 120_000, 3.50: 300_000, 4.00: 800_000, 4.50: 2_000_000}
TARGET_ERRORS = 1000
T_LIST = [21, 22, 23]
M = 2
RESUME_EVERY = 20_000

RESULTS_SINGLE = OUTDIR / "golay_t_sweep_m2_single_results.txt"
RESULTS_CUMULATIVE = OUTDIR / "golay_t_sweep_m2_cumulative_results.txt"
FIELDS = ["snr", "method", "t", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, t):
    return 970001 + t * 100 + snr_idx, 870001 + t * 100 + snr_idx


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


def append_rows(results_path, rows):
    existing = []
    if results_path.exists():
        with open(results_path, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (f"{float(r['snr']):.2f}", r["method"], str(r["t"]), str(int(float(r["m"]))))
    have = {key(r) for r in existing}
    for r in rows:
        if key(r) not in have:
            existing.append(r)
    with open(results_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for T in T_LIST:
        for i, snr in enumerate(SNR_LIST):
            t_tag = tag(snr)
            cap = SAMPLE_CAP[snr]
            cs, ms = seed_pair(i, T)
            log = OUTDIR / f"golay_t_sweep_m2_t{T}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(INI), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--checkpoint", str(T),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--m-list", str(M),
                "--monitor-every", "50000",
            ], f"golay_t_sweep_m2_t{T}_snr{t_tag}")
            samples = vals["samples"]
            row = [{"snr": snr, "method": "static_pm_rank", "t": T, "m": M,
                    "samples": samples, "errors": vals[f"static_errors[m={M}]"],
                    "bler": vals[f"static_bler[m={M}]"]}]
            append_rows(RESULTS_SINGLE, row)
            append_rows(RESULTS_CUMULATIVE, row)
            print(f"SNR={snr} t={T} done (cached={cached}, samples={samples:.0f})")

    print(f"updated {RESULTS_SINGLE} and {RESULTS_CUMULATIVE}")


if __name__ == "__main__":
    main()
