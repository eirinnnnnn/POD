#!/usr/bin/env python3
"""PED_M-SCL4 baseline for the t_sweep_m8_SCL_baseline trial: genuine
M-branch automorphism ensemble (not pruned from M=64), each branch decoded
with SCL list_size=4, where M is chosen so M*4 ~= the pure-SCL L value it's
compared against (L in {65,68,76,96} -> M in {16,17,19,24}). Plotted
alongside pure_scl(L) and the AED-64-pruned static/trace curves on
t_sweep_m8_bler.png."""
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
L_TO_M = {32: 8, 65: 16, 68: 17, 76: 19, 96: 24, 144: 36}
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "t_sweep_m8_PED_baseline_results.txt"
FIELDS = ["snr", "L", "M", "method", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 980001 + salt * 100 + snr_idx, 880001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", str(int(float(r["M"]))))
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
            log = OUTDIR / f"t_sweep_ped_M{M}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(BUILD / f"ebn0_sweep_m{M}.ini"), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"t_sweep_ped_M{M}_snr{t_tag}")
            append_rows([{"snr": snr, "L": L, "M": M, "method": "ped_m_scl4",
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} M={M} (L={L}) done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
