#!/usr/bin/env python3
"""Real-pickBest streaming re-measurement of the static-width PED baseline
(fixed m = round(oracle's own mean_m at that SNR), same value used in
plot_overlay.py's proxy-metric version) -- target-errors=1000, samples
capped at 100k, t=64, full SNR grid, both sides. Writes static_ped.txt."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
RESULTS = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")
TRACE_MODEL_T64 = RESULTS / "clean_trace_upto_t64_weights.txt"

# static m per SNR = round(oracle mean_m), same values already computed in
# pm_vs_model_out_results.txt (identical for both sides after rounding)
STATIC_M = {
    2.00: 8, 2.25: 5, 2.50: 3, 2.75: 2, 3.00: 2, 3.50: 1, 4.00: 1, 4.50: 1,
}
TARGET_ERRORS = 1000
SAMPLES_CAP = 100_000
FIELDS = ["eval_t", "snr", "static_m", "samples", "errors", "mean_m", "bler"]


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


def upsert_row(out_txt, row):
    existing = []
    if out_txt.exists():
        with open(out_txt, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: f"{float(r['snr']):.2f}"
    existing = [r for r in existing if key(r) != key(row)]
    existing.append(row)
    with open(out_txt, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", choices=["pm", "model_out"], required=True)
    args = ap.parse_args()

    if args.side == "model_out":
        sim = BUILD / "mclass_adaptive_m_bler_sim_model_out"
        out_txt = MO_DIR / "static_ped.txt"
        state_dir = MO_DIR / "data" / "static_ped"
    else:
        sim = BUILD / "mclass_adaptive_m_bler_sim"
        out_txt = PM_DIR / "static_ped.txt"
        state_dir = PM_DIR / "data" / "static_ped"
    state_dir.mkdir(parents=True, exist_ok=True)

    for i, (snr, m) in enumerate(sorted(STATIC_M.items())):
        cs = 1605001 + i
        ms = 1705001 + i
        resume = state_dir / f"t64_snr{int(snr*100)}.state"
        cmd = [str(sim), "-ini", str(INI), "--snr", str(snr), "--samples", str(SAMPLES_CAP),
               "--target-errors", str(TARGET_ERRORS), "--static-m", str(m),
               "--channel-seed", str(cs), "--message-seed", str(ms),
               "--monitor-every", "20000", "--resume-file", str(resume), "--resume-every", "20000"]
        if args.side == "model_out":
            cmd += ["--trace-model", str(TRACE_MODEL_T64)]
        else:
            cmd += ["--checkpoint", "64"]
        print(f"=== {args.side} t=64 SNR={snr} static_m={m} ===", flush=True)
        result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            print(f"FAILED SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
            sys.exit(1)
        vals = parse_summary(result.stdout)
        upsert_row(out_txt, {"eval_t": 64, "snr": snr, "static_m": m, "samples": vals["samples"],
                              "errors": vals["errors"], "mean_m": vals["mean_m"], "bler": vals["bler"]})
        print(f"{args.side} SNR={snr}: samples={vals['samples']:.0f} errors={vals['errors']:.0f} "
              f"bler={vals['bler']:.6f}", flush=True)
    print(f"STATIC_PED_{args.side.upper()}_DONE", flush=True)


if __name__ == "__main__":
    main()
