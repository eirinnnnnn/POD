#!/usr/bin/env python3
"""Complexity-matched PED_M-SCL4 reference for the t=40 and t=122 lines
added to t_sweep_m8_bler.png. M chosen so a full M-branch decode through
checkpoint t=128 costs the same as pruning 64->8 at checkpoint t:
M(t) = (64*t + 8*(128-t)) / 128  (k=8 fixed, matches the static/trace k).
t=40 -> M~25.5 -> M=27 (ebn0_sweep_m27.ini, already built/used in gain_drill).
t=122 -> M~61.4 -> M=61 (ebn0_sweep_m61.ini, new)."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
RESULTS_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 4.00: 1_500_000, 4.50: 1_000_000}
TARGET_ERRORS = 500
T_TO_M = {40: 27, 122: 61}
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "t_sweep_m8_PED_complexity_match_results.txt"
FIELDS = ["snr", "t", "M", "method", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 970001 + salt * 100 + snr_idx, 870001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", str(r["t"]), str(int(float(r["M"]))))
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
        for t, M in T_TO_M.items():
            cs, ms = seed_pair(i, salt=M)
            log = OUTDIR / f"t_sweep_ped_match_t{t}_M{M}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(RESULTS_DIR / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"t_sweep_ped_match_t{t}_M{M}_snr{t_tag}")
            append_rows([{"snr": snr, "t": t, "M": M, "method": "ped_m_scl4_matched",
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} t={t} M={M} done (cached={cached}, samples={vals['samples']:.0f})", flush=True)

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
