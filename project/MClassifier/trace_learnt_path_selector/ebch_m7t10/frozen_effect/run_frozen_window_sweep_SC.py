#!/usr/bin/env python3
"""Same as run_frozen_window_sweep.py (static-pm-rank-only, no model, two
5-checkpoint windows centered on the frozen bits at decode_idx=16 and
decode_idx=24), but with list_size=1 (plain SC per AED branch) instead
of the original list_size=4. Writes frozen_window_static_results_SC.txt."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
RESULTS_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
INI = RESULTS_DIR / "ebn0_sweep_SC.ini"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 3.75, 4.00, 4.25, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 50_000, 2.50: 80_000, 2.75: 150_000,
              3.00: 300_000, 3.50: 700_000, 3.75: 1_000_000, 4.00: 1_500_000,
              4.25: 1_000_000, 4.50: 1_000_000}
TARGET_ERRORS = 1000
TARGET_ERRORS_OVERRIDE = {4.25: 500, 4.50: 500}  # lower target at the two
# priciest SNR points -- 1000 was taking many hours per checkpoint there
CHECKPOINTS = [14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28]
K = 8
RESUME_EVERY = 20_000

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "frozen_window_static_results_SC.txt"
FIELDS = ["snr", "t", "method", "k", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    # distinct seed family from the L=4 trial (950001/850001+...) so the
    # two never collide even though they share the same channel/message
    # RNG formula shape.
    return 1950001 + salt * 100 + snr_idx, 1850001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", str(r["t"]), r["method"], str(int(float(r["k"]))))
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
        target_errors = TARGET_ERRORS_OVERRIDE.get(snr, TARGET_ERRORS)
        for t in CHECKPOINTS:
            cs, ms = seed_pair(i, t)
            log = OUTDIR / f"frozen_window_SC_t{t}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(INI), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(target_errors),
                "--checkpoint", str(t), "--k-list", str(K),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--monitor-every", "100000",
            ], f"frozen_window_SC_t{t}_snr{t_tag}")
            samples = vals["samples"]
            rows = [{"snr": snr, "t": t, "method": "static_pm_rank", "k": K,
                     "samples": samples, "errors": vals[f"static_errors[k={K}]"], "bler": vals[f"static_bler[k={K}]"]}]
            append_rows(rows)
            print(f"SNR={snr} t={t} done (cached={cached}, samples={samples:.0f})", flush=True)

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
