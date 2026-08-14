#!/usr/bin/env python3
"""Re-run the two event-starved SNR points (4.0, 4.5 dB) targeting 50
full-ensemble block errors instead of the 100k-sample cap used in the main
sweep. Only checkpoint t=64 is regenerated (full-PED is t-independent and
gets refreshed for free from the same bigger base dataset; t=8/t=32 rows in
results.csv are left as-is from the original run). Outputs use an "_e50"
suffix so the original 100k-cap artifacts aren't touched."""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_sweep import (  # noqa: E402
    BUILD, INI, MCLASSIFIER, OUTDIR, MODEL_FOR_T, CHECKPOINT_LADDER, K,
    PY, run, run_cached, parse_eval,
)

TARGETS = {
    # snr: (target_errors, sample_cap, channel_seed, message_seed)
    4.00: (50, 400_000, 901007, 801007),
    4.50: (50, 4_000_000, 901008, 801008),
}
T = 64


def main():
    results_path = OUTDIR / "results.csv"
    log_path = OUTDIR / "sweep.log"

    with open(results_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    with open(log_path, "a", encoding="utf-8") as log:
        for snr, (target_errors, sample_cap, channel_seed, message_seed) in TARGETS.items():
            t_tag = f"{snr:.2f}".replace(".", "p")
            log.write(f"\n===== sparse rerun SNR={snr} target_errors={target_errors} cap={sample_cap} =====\n")

            base_csv = OUTDIR / f"snr{t_tag}_metric_e50.csv"
            meta_path = Path(str(base_csv) + ".meta.json")
            run_cached([
                str(BUILD / "mclass_dataset"), "-ini", str(INI),
                "--snr", str(snr),
                "--target-errors", str(target_errors),
                "--samples", str(sample_cap),
                "--channel-seed", str(channel_seed),
                "--message-seed", str(message_seed),
                "--out", str(base_csv),
            ], log, base_csv, side_suffixes=[".meta.json"])

            meta = json.loads(meta_path.read_text())
            samples = meta["samples"]
            block_errors = meta["block_errors"]

            checkpoints_upto = [c for c in CHECKPOINT_LADDER if c <= T]
            trace_csv = OUTDIR / f"snr{t_tag}_trace_t{T}_e50.csv"
            run_cached([
                str(BUILD / "mclass_trace_dataset"), "-ini", str(INI),
                "--dataset", str(base_csv),
                "--checkpoints", ",".join(str(c) for c in checkpoints_upto),
                "--out", str(trace_csv),
            ], log, trace_csv)

            sp = OUTDIR / f"snr{t_tag}_static_t{T}_e50.csv"
            run_cached([PY, str(MCLASSIFIER / "architectures/static_pm_rank/static_path_metric_rank.py"),
                        "--trace", str(trace_csv), "--step", str(T), "--out", str(sp)], log, sp)

            tp = OUTDIR / f"snr{t_tag}_trace_t{T}_k{K}_e50.csv"
            run_cached([PY, str(MCLASSIFIER / "architectures/trace/predict_mclass_trace_ranker.py"),
                        "--trace", str(trace_csv), "--model", str(MODEL_FOR_T[T]),
                        "--topk", str(K), "--out", str(tp)], log, tp)

            out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                          "--pred", str(sp), "--topk", "64"], log)
            full_ped_bler = parse_eval(out)["topk_bler"]

            out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                          "--pred", str(sp), "--topk", str(K)], log)
            static_bler = parse_eval(out)["topk_bler"]

            out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                          "--pred", str(tp), "--topk", str(K)], log)
            trace_bler = parse_eval(out)["topk_bler"]

            # replace the full_ped and t=64 rows for this snr; leave t=8/t=32 untouched
            def keep(r):
                same_snr = float(r["snr"]) == snr
                is_full_ped = r["method"] == "full_ped"
                is_t64 = r["t"] == str(T)
                return not (same_snr and (is_full_ped or is_t64))

            rows = [r for r in rows if keep(r)]
            rows.append({"snr": snr, "t": "all", "method": "full_ped", "k": 64,
                         "samples": samples, "block_errors": block_errors, "bler": full_ped_bler})
            rows.append({"snr": snr, "t": T, "method": "static_pm_rank", "k": K,
                         "samples": samples, "block_errors": block_errors, "bler": static_bler})
            rows.append({"snr": snr, "t": T, "method": "trace_learned", "k": K,
                         "samples": samples, "block_errors": block_errors, "bler": trace_bler})

            with open(results_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["snr", "t", "method", "k", "samples", "block_errors", "bler"])
                writer.writeheader()
                writer.writerows(sorted(rows, key=lambda r: (float(r["snr"]), r["method"], str(r["t"]))))

            log.write(f"sparse rerun SNR={snr} done: samples={samples} block_errors={block_errors}\n")

    print(f"updated {results_path}")


if __name__ == "__main__":
    main()
