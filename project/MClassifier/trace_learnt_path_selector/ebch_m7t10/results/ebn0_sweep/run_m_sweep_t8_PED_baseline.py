#!/usr/bin/env python3
"""PED baselines for the m_sweep_t8 trial, two cases, both plotted on
m_sweep_t8_bler.png:

1. PED_M-SCL4 (kind=L_matched): "complexity match" -- genuine M-branch
   automorphism ensemble, SCL list_size=4, M chosen so M*4 equals a target
   total-path-complexity budget L (L in {45,55,95,149,200} -> M in
   {11,14,24,37,50}; L=200/M=50 is the complexity-matched baseline for the
   m_sweep_t8 trial's m=48 point).
2. PED_m-SCL4 (kind=m_matched): "m match" -- genuine m-branch automorphism
   ensemble, SCL list_size=4, for m in {1,4,16,32} -- the exact same
   pruning widths used in the m_sweep_t8 trial itself (there compared as
   top-m picks out of a pruned M=64 ensemble; here run as a genuine
   m-branch ensemble instead).
"""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 4.00: 1_500_000, 4.50: 1_000_000}
TARGET_ERRORS = 500
L_TO_M = {45: 11, 55: 14, 95: 24, 149: 37, 200: 50}
M_MATCH_MAP = {1: 1, 4: 4, 16: 16, 32: 32, 48: 48}
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "m_sweep_t8_PED_baseline_results.txt"
FIELDS = ["snr", "kind", "L", "M", "method", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 990001 + salt * 100 + snr_idx, 890001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", r["kind"], str(int(float(r["M"]))))
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

        for L, M in L_TO_M.items():
            cs, ms = seed_pair(i, salt=M)
            log = OUTDIR / f"m_sweep_ped_Lmatched_M{M}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"m_sweep_ped_Lmatched_M{M}_snr{t_tag}")
            append_rows([{"snr": snr, "kind": "L_matched", "L": L, "M": M, "method": "ped_m_scl4",
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} L_matched M={M} (L={L}) done (cached={cached}, samples={vals['samples']:.0f})")

        for m, M in M_MATCH_MAP.items():
            cs, ms = seed_pair(i, salt=100 + M)
            log = OUTDIR / f"m_sweep_ped_mmatched_M{M}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"m_sweep_ped_mmatched_M{M}_snr{t_tag}")
            append_rows([{"snr": snr, "kind": "m_matched", "L": "", "M": M, "method": "ped_m_scl4",
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} m_matched M={M} done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
