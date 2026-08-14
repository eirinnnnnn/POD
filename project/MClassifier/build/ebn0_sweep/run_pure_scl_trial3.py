#!/usr/bin/env python3
"""Second pure-SCL (no AED ensemble) diversity trial, L in {45, 55, 95, 149},
same setup as run_pure_scl_trial2.py (clean automorphism file, bug-fixed
mclass_bler_sim, same SNR grid / sample caps for comparability)."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 4.00: 1_500_000, 4.50: 1_000_000}
TARGET_ERRORS = 100
L_LIST = [45, 55, 95, 149]
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "m_sweep_t8_SCL_baseline_results.txt"
FIELDS = ["snr", "L", "method", "samples", "errors", "bler"]


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
    key = lambda r: (f"{float(r['snr']):.2f}", str(int(float(r["L"]))))
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
        for L in L_LIST:
            cs, ms = seed_pair(i, salt=L)
            log = OUTDIR / f"pure_scl3_L{L}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_scl{L}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"pure_scl3_L{L}_snr{t_tag}")
            append_rows([{"snr": snr, "L": L, "method": "pure_scl",
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} L={L} done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
