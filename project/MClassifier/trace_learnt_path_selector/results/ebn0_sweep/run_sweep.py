#!/usr/bin/env python3
"""Driver for the Eb/N0 sweep experiment: full-PED (all M=64 paths) vs
static-pm-rank vs trace-learned (k=8), at checkpoints t=8,32,64, across
SNR=2.00..4.50dB (step 0.25). Regenerates a fresh Monte-Carlo dataset per
SNR point (adaptive stop at 1000 full-ensemble block errors, capped at
200k samples), then predicts + evals each method.

Resumable: every file-producing step is skipped if its output already
exists, and SNR points already fully recorded in results.txt are skipped
entirely -- safe to re-run after an interruption."""
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

MCLASSIFIER = Path(__file__).resolve().parents[2]
BUILD = MCLASSIFIER / "build"
OUTDIR = BUILD / "ebn0_sweep"
INI = BUILD / "ebn0_sweep.ini"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
CHECKPOINTS = [8, 32, 64]
CHECKPOINT_LADDER = [8, 16, 32, 64, 128]
K = 8
TARGET_ERRORS = 1000
SAMPLE_CAP = 100000
RESULTS_FIELDS = ["snr", "t", "method", "k", "samples", "block_errors", "bler"]

MODEL_FOR_T = {
    t: BUILD / f"m7t10_m64_snr2p0_trace_upto_t{t}_basin_h64_e20.pt"
    for t in CHECKPOINTS
}

PY = sys.executable


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def run(cmd, log):
    log.write("+ " + " ".join(str(c) for c in cmd) + "\n")
    log.flush()
    t0 = time.time()
    result = subprocess.run(cmd, cwd=BUILD, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    dt = time.time() - t0
    log.write(result.stdout)
    log.write(f"[exit={result.returncode}, {dt:.1f}s]\n\n")
    log.flush()
    if result.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(str(c) for c in cmd)}")
    return result.stdout, dt


def run_cached(cmd, log, output_path, side_suffixes=()):
    """Run cmd unless output_path already exists. The command must take its
    output path as one argument (matched by equality to str(output_path));
    that argument is redirected to a .tmp path and atomically renamed into
    place only after a clean exit, so a kill/crash mid-write can never leave
    a truncated file that looks done. side_suffixes: extra suffixes appended
    to the output path for sibling files the same command derives from it
    (e.g. ".meta.json"), renamed alongside using the same tmp/final pairing."""
    if output_path.exists():
        log.write(f"[skip, exists] {output_path}\n")
        log.flush()
        return

    tmp_path = output_path.with_name(output_path.name + ".tmp")
    cmd = [str(tmp_path) if c == str(output_path) else c for c in cmd]

    tmp_path.unlink(missing_ok=True)
    side_pairs = [(Path(str(tmp_path) + suf), Path(str(output_path) + suf)) for suf in side_suffixes]
    for tmp_side, _ in side_pairs:
        tmp_side.unlink(missing_ok=True)

    run(cmd, log)

    os.replace(tmp_path, output_path)
    for tmp_side, final_side in side_pairs:
        os.replace(tmp_side, final_side)


def parse_eval(stdout):
    vals = {}
    for line in stdout.splitlines():
        if "," in line:
            k_, v_ = line.split(",", 1)
            vals[k_] = v_
    return vals


def load_existing_rows(results_path):
    if not results_path.exists():
        return [], set()
    with open(results_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    done_snr = set()
    for snr in SNR_LIST:
        count = sum(1 for r in rows if float(r["snr"]) == snr)
        if count == 1 + 2 * len(CHECKPOINTS):  # full_ped + (static+trace) x checkpoints
            done_snr.add(snr)
    return rows, done_snr


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results_path = OUTDIR / "results.txt"
    log_path = OUTDIR / "sweep.log"

    rows, done_snr = load_existing_rows(results_path)

    with open(log_path, "a", encoding="utf-8") as log:
        for i, snr in enumerate(SNR_LIST):
            if snr in done_snr:
                log.write(f"\n===== SNR={snr} already complete, skipping =====\n")
                log.flush()
                continue

            t_tag = tag(snr)
            log.write(f"\n===== SNR={snr} ({t_tag}) =====\n")
            channel_seed = 900001 + i
            message_seed = 800001 + i

            base_csv = OUTDIR / f"snr{t_tag}_metric.csv"
            meta_path = Path(str(base_csv) + ".meta.json")
            run_cached([
                str(BUILD / "mclass_dataset"), "-ini", str(INI),
                "--snr", str(snr),
                "--target-errors", str(TARGET_ERRORS),
                "--samples", str(SAMPLE_CAP),
                "--channel-seed", str(channel_seed),
                "--message-seed", str(message_seed),
                "--out", str(base_csv),
            ], log, base_csv, side_suffixes=[".meta.json"])

            meta = json.loads(meta_path.read_text())
            samples = meta["samples"]
            block_errors = meta["block_errors"]

            static_pred = {}
            trace_pred = {}
            for t in CHECKPOINTS:
                # Each "upto_t{t}" model was trained on the cumulative set of
                # checkpoints <= t (e.g. upto_t32 = [8,16,32]), not just t
                # itself -- match that exactly or input_dim mismatches.
                checkpoints_upto = [c for c in CHECKPOINT_LADDER if c <= t]
                trace_csv_t = OUTDIR / f"snr{t_tag}_trace_t{t}.csv"
                run_cached([
                    str(BUILD / "mclass_trace_dataset"), "-ini", str(INI),
                    "--dataset", str(base_csv),
                    "--checkpoints", ",".join(str(c) for c in checkpoints_upto),
                    "--out", str(trace_csv_t),
                ], log, trace_csv_t)

                sp = OUTDIR / f"snr{t_tag}_static_t{t}.csv"
                run_cached([PY, str(MCLASSIFIER / "architectures/static_pm_rank/static_path_metric_rank.py"),
                            "--trace", str(trace_csv_t), "--step", str(t), "--out", str(sp)], log, sp)
                static_pred[t] = sp

                tp = OUTDIR / f"snr{t_tag}_trace_t{t}_k{K}.csv"
                run_cached([PY, str(MCLASSIFIER / "architectures/trace/predict_mclass_trace_ranker.py"),
                            "--trace", str(trace_csv_t), "--model", str(MODEL_FOR_T[t]),
                            "--topk", str(K), "--out", str(tp)], log, tp)
                trace_pred[t] = tp

            # full-PED baseline: k=64 (all M paths), t-independent, reuse any full ranking file
            out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                          "--pred", str(static_pred[CHECKPOINTS[0]]), "--topk", "64"], log)
            vals = parse_eval(out)
            rows.append({"snr": snr, "t": "all", "method": "full_ped", "k": 64,
                         "samples": samples, "block_errors": block_errors,
                         "bler": vals["topk_bler"]})

            for t in CHECKPOINTS:
                out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                              "--pred", str(static_pred[t]), "--topk", str(K)], log)
                vals = parse_eval(out)
                rows.append({"snr": snr, "t": t, "method": "static_pm_rank", "k": K,
                             "samples": samples, "block_errors": block_errors,
                             "bler": vals["topk_bler"]})

                out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                              "--pred", str(trace_pred[t]), "--topk", str(K)], log)
                vals = parse_eval(out)
                rows.append({"snr": snr, "t": t, "method": "trace_learned", "k": K,
                             "samples": samples, "block_errors": block_errors,
                             "bler": vals["topk_bler"]})

            # write results incrementally so progress is inspectable mid-run
            with open(results_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=RESULTS_FIELDS)
                writer.writeheader()
                writer.writerows(rows)

            log.write(f"SNR={snr} done: samples={samples} block_errors={block_errors}\n")

    print(f"wrote {results_path}")


if __name__ == "__main__":
    main()
