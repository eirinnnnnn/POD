#!/usr/bin/env python3
"""Low-SNR refine pass for the m_sweep_t8_PED_baseline trial: SNR
2.00-3.00dB only, target-errors bumped 500->1000 to smooth out the
noisiest, cheapest region of both the L_matched (M=11,14,24,37) and
m_matched (M=1,4,16,32) curves. Unlike the main driver's append_rows, this
REPLACES any existing row for the same (snr, kind, M) key rather than
skipping it. Meant to be run only after the main sweep/baseline trials
finish, to avoid competing for CPU."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000, 3.00: 300_000}
TARGET_ERRORS = 1000
L_TO_M = {45: 11, 55: 14, 95: 24, 149: 37}
M_LIST = [1, 4, 16, 32]
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "m_sweep_t8_PED_baseline_results.txt"
FIELDS = ["snr", "kind", "L", "M", "method", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 960001 + salt * 100 + snr_idx, 860001 + salt * 100 + snr_idx


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


def replace_rows(rows):
    existing = []
    if RESULTS.exists():
        with open(RESULTS, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (f"{float(r['snr']):.2f}", r["kind"], str(int(float(r["M"]))))
    new_keys = {key(r) for r in rows}
    existing = [r for r in existing if key(r) not in new_keys]
    existing.extend(rows)
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
            log = OUTDIR / f"m_sweep_ped_Lmatched_M{M}_snr{t_tag}_refine1000.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"m_sweep_ped_Lmatched_M{M}_snr{t_tag}_refine1000")
            replace_rows([{"snr": snr, "kind": "L_matched", "L": L, "M": M, "method": "ped_m_scl4",
                           "samples": vals["samples"], "errors": vals["full_ped_errors"],
                           "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} L_matched M={M} (L={L}) refined (cached={cached}, samples={vals['samples']:.0f})")

        for M in M_LIST:
            cs, ms = seed_pair(i, salt=100 + M)
            log = OUTDIR / f"m_sweep_ped_mmatched_M{M}_snr{t_tag}_refine1000.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"m_sweep_ped_mmatched_M{M}_snr{t_tag}_refine1000")
            replace_rows([{"snr": snr, "kind": "m_matched", "L": "", "M": M, "method": "ped_m_scl4",
                           "samples": vals["samples"], "errors": vals["full_ped_errors"],
                           "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} m_matched M={M} refined (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
