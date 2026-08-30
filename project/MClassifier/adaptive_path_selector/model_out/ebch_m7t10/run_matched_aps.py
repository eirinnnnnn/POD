#!/usr/bin/env python3
"""Evaluate the matched-t APS selectors (t in {8,32,96}, each trained on
its own checkpoint's data) at their own matching checkpoint, across the
same SNR grid as generalizability.txt (target-errors=500), so the result
can be overlaid against the t=64-trained (mismatched, for t!=64) selector
already in generalizability.txt. Writes matched_aps.txt on each side --
same streaming, no-CSV approach as run_generalizability.py."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build")
INI = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results/ebn0_sweep.ini")
RESULTS = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/results")
GAIN_DRILL = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/trace_learnt_path_selector/ebch_m7t10/gain_drill")
MO_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/model_out/ebch_m7t10")
PM_DIR = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/adaptive_path_selector/pm/ebch_m7t10")

TRACE_MODELS = {
    8: RESULTS / "clean_trace_upto_t8_weights.txt",
    32: RESULTS / "clean_trace_upto_t32_weights.txt",
    96: GAIN_DRILL / "t96_weights.txt",
}
T_LIST = [8, 32, 96]
SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
TARGET_ERRORS = 500
SAMPLES_CAP = 1_000_000
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]


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


def already_done(out_txt, t, snr):
    if not out_txt.exists():
        return False
    with open(out_txt, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return True
    return False


def upsert_row(out_txt, row):
    existing = []
    if out_txt.exists():
        with open(out_txt, newline="", encoding="utf-8") as f:
            existing = list(csv.DictReader(f))
    key = lambda r: (str(r["eval_t"]), f"{float(r['snr']):.2f}")
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
        out_txt = MO_DIR / "matched_aps.txt"
        state_dir = MO_DIR / "data" / "matched_aps"
    else:
        sim = BUILD / "mclass_adaptive_m_bler_sim"
        out_txt = PM_DIR / "matched_aps.txt"
        state_dir = PM_DIR / "data" / "matched_aps"
    state_dir.mkdir(parents=True, exist_ok=True)

    for t in T_LIST:
        selector = (MO_DIR if args.side == "model_out" else PM_DIR) / f"weights_t{t}_snr2p0.txt"
        for i, snr in enumerate(SNR_LIST):
            if already_done(out_txt, t, snr):
                print(f"{args.side} t={t} SNR={snr}: already in matched_aps.txt, skipping", flush=True)
                continue
            cs = 1105001 + t * 10 + i
            ms = 1005001 + t * 10 + i
            resume = state_dir / f"t{t}_snr{int(snr*100)}.state"
            cmd = [str(sim), "-ini", str(INI), "--snr", str(snr), "--samples", str(SAMPLES_CAP),
                   "--target-errors", str(TARGET_ERRORS), "--selector-model", str(selector),
                   "--channel-seed", str(cs), "--message-seed", str(ms),
                   "--monitor-every", "20000", "--resume-file", str(resume), "--resume-every", "20000"]
            if args.side == "model_out":
                cmd += ["--trace-model", str(TRACE_MODELS[t])]
            else:
                cmd += ["--checkpoint", str(t)]
            print(f"=== {args.side} t={t} SNR={snr} (matched selector) ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)
            vals = parse_summary(result.stdout)
            upsert_row(out_txt, {"eval_t": t, "snr": snr, "samples": vals["samples"], "errors": vals["errors"],
                                  "mean_m": vals["mean_m"], "bler": vals["bler"]})
            print(f"{args.side} t={t} SNR={snr} done: samples={vals['samples']:.0f} errors={vals['errors']:.0f} "
                  f"mean_m={vals['mean_m']:.3f} bler={vals['bler']:.6f}", flush=True)
    print(f"MATCHED_APS_{args.side.upper()}_DONE", flush=True)


if __name__ == "__main__":
    main()
