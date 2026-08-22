#!/usr/bin/env python3
"""Fill exactly the (t,snr) gaps NOT already present in cross_t_results.txt
(the data behind cross_t_bler_vs_snr.png), restricted to t in
{8,16,32,48,64,96}, using the new streaming tool (real pickBest
correctness) with --target-errors 1000 and a 1,000,000-sample cap. Writes
into cross_t_stream_results.txt (separate from the old proxy-metric file,
since the two measure different things -- see conversation). Run for both
sides via --side {pm,model_out}."""
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
    48: GAIN_DRILL / "t48_weights.txt",
    64: RESULTS / "clean_trace_upto_t64_weights.txt",
    96: GAIN_DRILL / "t96_weights.txt",
}
T_LIST = [8, 16, 32, 48, 64, 96]
SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLES_CAP = 1_000_000
TARGET_ERRORS = 1000

OLD_RESULTS_MO = MO_DIR / "cross_t_results.txt"  # old proxy-metric file (backs cross_t_bler_vs_snr.png)
FIELDS = ["eval_t", "snr", "samples", "errors", "mean_m", "bler"]


def old_covered(side, t, snr):
    if not OLD_RESULTS_MO.exists():
        return False
    with open(OLD_RESULTS_MO, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["side"] == side and int(r["eval_t"]) == t and abs(float(r["snr"]) - snr) < 1e-6:
                return True
    return False


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


def append_row(out_txt, row):
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


def already_done(out_txt, t, snr):
    if not out_txt.exists():
        return False
    with open(out_txt, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if str(r["eval_t"]) == str(t) and abs(float(r["snr"]) - snr) < 1e-6:
                return True
    return False


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

    for t in T_LIST:
        for i, snr in enumerate(SNR_LIST):
            if old_covered(args.side, t, snr):
                print(f"{args.side} t={t} SNR={snr}: already in cross_t_bler_vs_snr.png data, skipping", flush=True)
                continue
            if already_done(out_txt, t, snr):
                print(f"{args.side} t={t} SNR={snr}: already streamed this session, skipping", flush=True)
                continue
            cs = 995001 + t * 10 + i
            ms = 895001 + t * 10 + i
            resume = state_dir / f"gapfill_t{t}_snr{int(snr*100)}.state"
            cmd = [str(sim), "-ini", str(INI), "--snr", str(snr), "--samples", str(SAMPLES_CAP),
                   "--target-errors", str(TARGET_ERRORS), "--selector-model", str(selector),
                   "--channel-seed", str(cs), "--message-seed", str(ms),
                   "--monitor-every", "20000", "--resume-file", str(resume), "--resume-every", "20000"]
            if args.side == "model_out":
                cmd += ["--trace-model", str(TRACE_MODELS[t])]
            else:
                cmd += ["--checkpoint", str(t)]
            print(f"=== {args.side} t={t} SNR={snr} ===", flush=True)
            result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode != 0:
                print(f"FAILED t={t} SNR={snr}:\n{result.stdout[-2000:]}", flush=True)
                sys.exit(1)
            vals = parse_summary(result.stdout)
            append_row(out_txt, {"eval_t": t, "snr": snr, "samples": vals["samples"], "errors": vals["errors"],
                                  "mean_m": vals["mean_m"], "bler": vals["bler"]})
            print(f"{args.side} t={t} SNR={snr} done: samples={vals['samples']:.0f} errors={vals['errors']:.0f} "
                  f"mean_m={vals['mean_m']:.3f} bler={vals['bler']:.6f}", flush=True)
    print(f"GAP_FILL_{args.side.upper()}_DONE", flush=True)


if __name__ == "__main__":
    main()
