#!/usr/bin/env python3
"""Trial 1 driver.

1A: fix t=64, sweep k in {1,4,8,32} for static-pm-rank and trace-learned
    (one mclass_bler_sim pass per SNR gives all k's at once), plus genuine
    PED_k baselines (k=1,4,32; PED8/PED64 already exist from the main sweep).
2B: pure-SCL list-size diversity reference, no AED ensemble, L in
    {4,16,32,128}.

SNR grid and error target are deliberately lighter than the main sweep
(coarser grid, target 200 errors) for a faster first pass -- see
EXPERIMENT_NOTES.txt. Resumable: each run is skipped if its log already has
a completed summary block.
"""
import csv
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
OUTDIR = BUILD / "ebn0_sweep"
SIM = BUILD / "mclass_bler_sim"
MODEL_T64 = BUILD / "trace_upto_t64_weights.txt"

SNR_LIST = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
TARGET_ERRORS = 200
SAMPLE_CAP = {2.0: 30_000, 2.5: 60_000, 3.0: 150_000, 3.5: 400_000, 4.0: 1_000_000, 4.5: 3_000_000}
K_LIST = [1, 4, 8, 32]
RESUME_EVERY = 20_000  # mclass_bler_sim checkpoints its running state at this cadence,
                       # so a kill mid-run loses at most this many samples, not the whole config

RESULTS = OUTDIR / "trial1_results.txt"
FIELDS = ["trial", "snr", "config", "method", "k", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr, salt):
    i = SNR_LIST.index(snr)
    return 920001 + salt * 100 + i, 820001 + salt * 100 + i


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


def run_sim(log_path, args):
    if log_path.exists():
        text = log_path.read_text()
        if "\nsamples," in text:
            return parse_summary(text), True
    resume_path = log_path.with_suffix(".state")
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
    # Normalize types before comparing -- CSV-read rows have every field as a
    # str, freshly-appended rows have snr as float / k as int, so an
    # unnormalized key never matches across separate script invocations and
    # silently fails to dedup (this bit us: exact triplicate rows after
    # three restarts before this fix).
    key = lambda r: (r["trial"], f"{float(r['snr']):.2f}", r["config"], r["method"], str(int(float(r["k"]))))
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

    for snr in SNR_LIST:
        t_tag = tag(snr)
        cap = SAMPLE_CAP[snr]
        cs, ms = seed_pair(snr, salt=1)

        # --- Trial 1A: k-sweep at t=64, static + trace-learned ---
        log = OUTDIR / f"trial1a_snr{t_tag}.log"
        vals, cached = run_sim(log, [
            "-ini", "ebn0_sweep.ini", "--snr", str(snr),
            "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
            "--stop-metric", "static_trace_min",
            "--channel-seed", str(cs), "--message-seed", str(ms),
            "--checkpoint", "64", "--k-list", ",".join(str(k) for k in K_LIST),
            "--model", str(MODEL_T64), "--monitor-every", "50000",
        ])
        samples = vals["samples"]
        rows = [{"trial": "1A", "snr": snr, "config": "full_ped_m64", "method": "full_ped", "k": 64,
                 "samples": samples, "errors": vals["full_ped_errors"], "bler": vals["full_ped_bler"]}]
        for k in K_LIST:
            rows.append({"trial": "1A", "snr": snr, "config": "prune", "method": "static_pm_rank", "k": k,
                         "samples": samples, "errors": vals[f"static_errors[k={k}]"],
                         "bler": vals[f"static_bler[k={k}]"]})
            rows.append({"trial": "1A", "snr": snr, "config": "prune", "method": "trace_learned", "k": k,
                         "samples": samples, "errors": vals[f"trace_errors[k={k}]"],
                         "bler": vals[f"trace_bler[k={k}]"]})
        append_rows(rows)
        print(f"SNR={snr} trial1A done (cached={cached}, samples={samples:.0f})")

        # --- Trial 1A baselines: genuine PED_k for k=1,4,32 ---
        for i, m in enumerate((1, 4, 32)):
            cs2, ms2 = seed_pair(snr, salt=2 + i)
            log = OUTDIR / f"trial1a_ped{m}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", f"ebn0_sweep_m{m}.ini", "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs2), "--message-seed", str(ms2),
                "--monitor-every", "50000",
            ])
            append_rows([{"trial": "1A", "snr": snr, "config": f"full_ped_m{m}", "method": "full_ped", "k": m,
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} PED{m} done (cached={cached}, samples={vals['samples']:.0f})")

        # --- Trial 1B: pure SCL list-size diversity, no AED ---
        # L=128 dropped -- disproportionately expensive (2h38min alone at
        # 2.0dB/30k-sample cap vs ~26min for everything else combined) for a
        # first trend-inspection pass; already have 2.0/2.5/3.0 data for it.
        for i, L in enumerate((4, 16, 32)):
            cs3, ms3 = seed_pair(snr, salt=10 + i)
            log = OUTDIR / f"trial1b_scl{L}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", f"ebn0_sweep_scl{L}.ini", "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--channel-seed", str(cs3), "--message-seed", str(ms3),
                "--monitor-every", "50000",
            ])
            append_rows([{"trial": "1B", "snr": snr, "config": f"scl{L}_m1", "method": "pure_scl", "k": L,
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} SCL{L} done (cached={cached}, samples={vals['samples']:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
