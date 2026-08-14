#!/usr/bin/env python3
"""Golay24 static-pm-rank sweep, SC-base/M=16 variant: base per-branch
decoder weakened to plain SC (list_size=1) and the automorphism ensemble
shrunk to M=16, specifically to surface early-path-selection gain -- at
list_size=4/M=64 (see run_golay_m_sweep.py), branches mostly agree at
SNR>=2.0dB since they use genuine code automorphisms (branches that
converge to the correct codeword necessarily report identical metrics),
so pruning to m<M loses nothing until BLER is high enough for branches to
actually disagree. SNR region narrowed to [0.0, 0.5, 3.0] -- the first two
probe the disagreement regime, 3.0 kept as a reference point already known
(from the M=64/SCL4 run) to show near-total overlap."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
INI = BUILD / "golay_sc_sweep.ini"
OUTDIR = BUILD / "sweep"

SNR_LIST = [0.0, 0.5, 3.0]
SAMPLE_CAP = {0.0: 5_000, 0.5: 10_000, 3.0: 100_000}
TARGET_ERRORS = 500
M_LIST = [1, 4, 8, 12]
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "golay_sc_m16_results.txt"
FIELDS = ["snr", "method", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx):
    return 710001 + snr_idx, 610001 + snr_idx


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
        cs, ms = seed_pair(i)
        log = OUTDIR / f"golay_sc_m16_snr{t_tag}.log"
        vals, cached = run_sim(log, [
            "-ini", str(INI), "--snr", str(snr),
            "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
            "--channel-seed", str(cs), "--message-seed", str(ms),
            "--m-list", ",".join(str(m) for m in M_LIST),
            "--monitor-every", "100000",
        ], f"golay_sc_m16_snr{t_tag}")
        samples = vals["samples"]
        rows = [{"snr": snr, "method": "full_ped", "m": 16,
                 "samples": samples, "errors": vals["full_ped_errors"], "bler": vals["full_ped_bler"]}]
        for m in M_LIST:
            rows.append({"snr": snr, "method": "static_pm_rank", "m": m,
                         "samples": samples, "errors": vals[f"static_errors[m={m}]"],
                         "bler": vals[f"static_bler[m={m}]"]})
        append_rows(rows)
        print(f"SNR={snr} done (cached={cached}, samples={samples:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
