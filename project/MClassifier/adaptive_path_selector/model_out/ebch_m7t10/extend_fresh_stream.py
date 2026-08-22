#!/usr/bin/env python3
"""Extend fresh (new streaming-tool, 985001/885001-seeded) points that
currently have fewer than 1000 errors, up to --target-errors 1000 (capped
at 1,000,000 samples). Reuses the SAME seed formula run_cross_t_stream.py
uses, so re-running reproduces the same initial sample sequence and
continues naturally. Uses the REAL streaming tool (pickBest correctness),
consistent with how these rows were originally produced."""
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
    16: RESULTS / "clean_trace_upto_t16_weights.txt",
    32: RESULTS / "clean_trace_upto_t32_weights.txt",
    40: GAIN_DRILL / "t40_weights.txt",
    48: GAIN_DRILL / "t48_weights.txt",
    64: RESULTS / "clean_trace_upto_t64_weights.txt",
    96: GAIN_DRILL / "t96_weights.txt",
    122: GAIN_DRILL / "t122_weights.txt",
}
SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
TARGET_ERRORS = 1000
SAMPLES_CAP = 1_000_000
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]

TARGETS_BY_SIDE = {
    "model_out": [
        (8, 2.00), (8, 2.25), (8, 2.50), (8, 2.75), (8, 3.00), (8, 3.50), (8, 4.00),
        (16, 3.50), (16, 4.50), (40, 4.00),
    ],
    "pm": [
        (8, 2.00), (8, 2.25), (8, 2.50), (8, 2.75), (8, 3.00), (8, 3.50), (8, 4.00), (8, 4.50),
        (16, 2.00), (16, 2.25), (16, 2.50), (16, 2.75), (16, 3.00), (16, 3.50), (16, 4.00), (16, 4.50),
        (32, 2.00), (32, 2.25), (32, 2.50), (32, 2.75), (32, 3.00), (32, 3.50), (32, 4.00), (32, 4.50),
        (40, 2.00), (40, 2.25), (40, 2.50), (40, 2.75), (40, 3.00), (40, 3.50), (40, 4.00), (40, 4.50),
        (48, 2.25), (48, 2.50), (48, 2.75), (48, 3.50), (48, 4.00),
    ],
}


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
    key = lambda r: (str(r["eval_t"]), f"{float(r['snr']):.2f}")
    existing = [r for r in existing if key(r) != key(row)]
    existing.append(row)
    with open(out_txt, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def current_errors(out_txt, t, snr):
    if not out_txt.exists():
        return None
    with open(out_txt, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return float(r["errors"])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", choices=["pm", "model_out"], required=True)
    args = ap.parse_args()

    if args.side == "model_out":
        sim = BUILD / "mclass_adaptive_m_bler_sim_model_out"
        selector = MO_DIR / "weights_t64_snr2p0.txt"
        out_txt = MO_DIR / "cross_t_stream_results.txt"
        state_dir = MO_DIR / "data" / "cross_t_stream"
    else:
        sim = BUILD / "mclass_adaptive_m_bler_sim"
        selector = PM_DIR / "weights_t64_snr2.txt"
        out_txt = PM_DIR / "cross_t_stream_results.txt"
        state_dir = PM_DIR / "data" / "cross_t_stream"
    state_dir.mkdir(parents=True, exist_ok=True)

    for t, snr in TARGETS_BY_SIDE[args.side]:
        i = SNR_LIST.index(snr)
        cur_err = current_errors(out_txt, t, snr)
        if cur_err is None or cur_err >= TARGET_ERRORS:
            print(f"{args.side} t={t} SNR={snr}: already >=1000 or missing, skipping", flush=True)
            continue
        cs = 985001 + t * 10 + i
        ms = 885001 + t * 10 + i
        resume = state_dir / f"t{t}_snr{int(snr*100)}.state"
        cmd = [str(sim), "-ini", str(INI), "--snr", str(snr), "--samples", str(SAMPLES_CAP),
               "--target-errors", str(TARGET_ERRORS), "--selector-model", str(selector),
               "--channel-seed", str(cs), "--message-seed", str(ms),
               "--monitor-every", "20000", "--resume-file", str(resume), "--resume-every", "20000"]
        if args.side == "model_out":
            cmd += ["--trace-model", str(TRACE_MODELS[t])]
        else:
            cmd += ["--checkpoint", str(t)]
        print(f"=== {args.side} t={t} SNR={snr} (was {cur_err:.0f} errors) ===", flush=True)
        result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            print(f"FAILED t={t} SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
            sys.exit(1)
        vals = parse_summary(result.stdout)
        upsert_row(out_txt, {"eval_t": t, "snr": snr, "samples": vals["samples"], "errors": vals["errors"],
                              "mean_m": vals["mean_m"], "bler": vals["bler"]})
        print(f"{args.side} t={t} SNR={snr} extended: samples={vals['samples']:.0f} errors={vals['errors']:.0f} "
              f"mean_m={vals['mean_m']:.3f} bler={vals['bler']:.6f}", flush=True)
    print(f"EXTEND_FRESH_{args.side.upper()}_DONE", flush=True)


if __name__ == "__main__":
    main()
