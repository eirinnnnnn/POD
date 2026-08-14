#!/usr/bin/env python3
"""Golay24 static-pm-rank BLER-vs-Eb/N0 sweep (Phase 1: static-pm-rank only,
no trace-learned -- see MClassifier/Golay24 for context). Analogous to the
eBCH m7t10 t_sweep_m8/m_sweep_t8 trials, but with a single dimension (m,
the AED branch-pruning width out of M=64) since Golay24's decoder has no
checkpoint/trace introspection."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
INI = BUILD / "golay_sweep.ini"
OUTDIR = BUILD / "sweep"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
# Scaled down from the eBCH caps: Golay24 decodes ~80 samples/sec (every
# sample fully decodes all M=64 branches regardless of m), vs eBCH's much
# higher throughput -- keeps total sweep runtime to roughly an hour.
SAMPLE_CAP = {2.00: 3_000, 2.25: 5_000, 2.50: 8_000, 2.75: 12_000,
              3.00: 20_000, 3.50: 40_000, 4.00: 80_000, 4.50: 150_000}
TARGET_ERRORS = 500
M_LIST = [1, 4, 8, 16, 32]
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "golay_m_sweep_results.txt"
FIELDS = ["snr", "method", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx):
    return 700001 + snr_idx, 600001 + snr_idx


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
        log = OUTDIR / f"golay_sweep_snr{t_tag}.log"
        vals, cached = run_sim(log, [
            "-ini", str(INI), "--snr", str(snr),
            "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
            "--channel-seed", str(cs), "--message-seed", str(ms),
            "--m-list", ",".join(str(m) for m in M_LIST),
            "--monitor-every", "100000",
        ], f"golay_sweep_snr{t_tag}")
        samples = vals["samples"]
        rows = [{"snr": snr, "method": "full_ped", "m": 64,
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
