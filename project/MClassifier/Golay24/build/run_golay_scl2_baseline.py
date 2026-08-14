#!/usr/bin/env python3
"""Pure-SCL(L=2) no-AED baseline, added to the golay_t_sweep_m2 trial's
results file alongside pure_scl8/pure_sc."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
INI = BUILD / "golay_scl2_pure.ini"
OUTDIR = BUILD / "sweep"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 45_000, 2.50: 55_000, 2.75: 85_000,
              3.00: 120_000, 3.50: 300_000, 4.00: 800_000, 4.50: 2_000_000}
TARGET_ERRORS = 1000
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "golay_t_sweep_m2_results.txt"
FIELDS = ["snr", "method", "t", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx):
    return 770001 + snr_idx, 670001 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", r["method"], str(r["t"]), str(int(float(r["m"]))))
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
        cs, ms = seed_pair(i)
        log = OUTDIR / f"golay_pure_scl2_m2trial_snr{t_tag}.log"
        vals, cached = run_sim(log, [
            "-ini", str(INI), "--snr", str(snr),
            "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
            "--checkpoint", "8",
            "--channel-seed", str(cs), "--message-seed", str(ms),
            "--m-list", "1",
            "--monitor-every", "50000",
        ], f"golay_pure_scl2_m2trial_snr{t_tag}")
        append_rows([{"snr": snr, "method": "pure_scl2", "t": "", "m": 1,
                      "samples": vals["samples"], "errors": vals["full_ped_errors"],
                      "bler": vals["full_ped_bler"]}])
        print(f"SNR={snr} pure_scl2 done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
