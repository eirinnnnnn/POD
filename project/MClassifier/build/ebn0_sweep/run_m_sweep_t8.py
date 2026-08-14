#!/usr/bin/env python3
"""'M choose m' early-path-selecting-PED sweep: fixed checkpoint t=8, fixed
SCL list_size L=4 (ebn0_sweep.ini as-is), vary the pruning width m (formerly
called k) in {1, 4, 16, 32}, static-pm-rank vs trace-learned. Same SNR grid
and sample caps as run_clean_sweep.py for direct comparability; full_ped
(M=64) reference is reused from that run rather than recomputed.
"""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
INI = BUILD / "ebn0_sweep.ini"
MODEL_T8 = BUILD / "clean_trace_upto_t8_weights.txt"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 4.00: 1_500_000, 4.50: 1_000_000}
TARGET_ERRORS = 1000
T = 8
M_LIST = [1, 4, 16, 32, 48]
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "m_sweep_t8_results.txt"
FIELDS = ["snr", "t", "method", "m", "samples", "errors", "bler"]


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


def append_rows(rows):
    existing = []
    if RESULTS.exists():
        with open(RESULTS, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (f"{float(r['snr']):.2f}", str(r["t"]), r["method"], str(int(float(r["m"]))))
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
        cs, ms = seed_pair(i, salt=T)
        log = OUTDIR / f"m_sweep_t{T}_snr{t_tag}.log"
        vals, cached = run_sim(log, [
            "-ini", str(INI), "--snr", str(snr),
            "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
            "--stop-metric", "static_trace_min",
            "--channel-seed", str(cs), "--message-seed", str(ms),
            "--checkpoint", str(T), "--k-list", ",".join(str(m) for m in M_LIST),
            "--model", str(MODEL_T8), "--monitor-every", "100000",
        ], f"m_sweep_t{T}_snr{t_tag}")
        samples = vals["samples"]
        rows = []
        for m in M_LIST:
            rows.append({"snr": snr, "t": T, "method": "static_pm_rank", "m": m,
                         "samples": samples, "errors": vals[f"static_errors[k={m}]"],
                         "bler": vals[f"static_bler[k={m}]"]})
            rows.append({"snr": snr, "t": T, "method": "trace_learned", "m": m,
                         "samples": samples, "errors": vals[f"trace_errors[k={m}]"],
                         "bler": vals[f"trace_bler[k={m}]"]})
        append_rows(rows)
        print(f"SNR={snr} t={T} m-sweep done (cached={cached}, samples={samples:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
