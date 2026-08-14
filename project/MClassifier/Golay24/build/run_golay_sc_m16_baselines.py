#!/usr/bin/env python3
"""Two no-AED reference baselines for the golay_sc_m16 trial, at the same
SNR region [0.0, 0.5, 3.0]: pure SCL(L=8) as a rough near-ML reference, and
pure SC (list_size=1) as the weak-decoder floor. Both written into
golay_sc_m16_results.txt alongside run_golay_sc_m16_sweep.py's output."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
OUTDIR = BUILD / "sweep"

CONFIGS = [
    ("pure_scl8", BUILD / "golay_scl8.ini"),
    ("pure_sc", BUILD / "golay_sc_pure.ini"),
]

SNR_LIST = [0.0, 0.5, 3.0]
SAMPLE_CAP = {0.0: 5_000, 0.5: 10_000, 3.0: 100_000}
TARGET_ERRORS = 500
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "golay_sc_m16_results.txt"
FIELDS = ["snr", "method", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 720001 + salt * 100 + snr_idx, 620001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", r["method"], str(int(float(r["m"]))))
    have = {key(r) for r in existing}
    for r in rows:
        if key(r) not in have:
            existing.append(r)
    with open(RESULTS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for i, snr in enumerate(SNR_LIST):
        t_tag = tag(snr)
        cap = SAMPLE_CAP[snr]
        for salt, (method, ini) in enumerate(CONFIGS):
            cs, ms = seed_pair(i, salt)
            log = OUTDIR / f"golay_{method}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(ini), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--m-list", "1",
                "--monitor-every", "100000",
            ], f"golay_{method}_snr{t_tag}")
            append_rows([{"snr": snr, "method": method, "m": 1,
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} {method} done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
