#!/usr/bin/env python3
"""Extend the static-pm-rank (model-free, k=8) curve to t=80,96,112 to see
where it closes the gap to full-PED M=64 -- no training involved, so this
is cheap to explore before deciding whether new trace-ranker checkpoints
are worth training. Reuses each SNR's existing base metric dataset (the
larger _e50 dataset for 4.0/4.5 dB where available)."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_sweep import BUILD, INI, MCLASSIFIER, OUTDIR, K, PY, SNR_LIST, run, run_cached, parse_eval  # noqa: E402

NEW_T = [80, 96]  # 112 dropped -- was ~5-8% off full-PED same as 96, not worth the extra compute
CHECKPOINT_LADDER = [8, 16, 32, 64, 80, 96, 112, 128]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


SPARSE_SNR = {4.00, 4.50}  # these got a dedicated bigger-sample _e50 rerun; wait for it rather than use the sparse original


def pick_base_csv(snr):
    t_tag = tag(snr)
    e50 = OUTDIR / f"snr{t_tag}_metric_e50.csv"
    # Check the actual data file, not just its .meta.json -- a meta.json can
    # survive (e.g. after an accidental cleanup) while the multi-GB CSV it
    # describes does not, and that must NOT be mistaken for "ready".
    if e50.exists():
        return e50
    if snr in SPARSE_SNR:
        return None  # e50 not available (still in flight, or abandoned) -- don't fall back to the sparse original
    plain = OUTDIR / f"snr{t_tag}_metric.csv"
    if plain.exists():
        return plain
    return None


def main():
    results_path = OUTDIR / "results.txt"
    log_path = OUTDIR / "sweep_extra_t.log"

    with open(results_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    have = {(float(r["snr"]), int(r["t"])) for r in rows
            if r["method"] == "static_pm_rank" and r["t"] != "all"}

    with open(log_path, "a", encoding="utf-8") as log:
        for snr in SNR_LIST:
            t_tag = tag(snr)
            base_csv = pick_base_csv(snr)
            if base_csv is None:
                log.write(f"\n===== SNR={snr}: no base dataset yet, skipping =====\n")
                log.flush()
                continue

            suffix = "_e50" if "_e50" in base_csv.name else ""
            checkpoints_upto = [c for c in CHECKPOINT_LADDER if c <= max(NEW_T)]
            trace_csv = OUTDIR / f"snr{t_tag}_trace_extra_t{suffix}.csv"
            log.write(f"\n===== SNR={snr} extra-t (base={base_csv.name}) =====\n")
            run_cached([
                str(BUILD / "mclass_trace_dataset"), "-ini", str(INI),
                "--dataset", str(base_csv),
                "--checkpoints", ",".join(str(c) for c in checkpoints_upto),
                "--out", str(trace_csv),
            ], log, trace_csv)

            import json
            samples = None
            block_errors = None
            meta_path = Path(str(base_csv) + ".meta.json")
            if meta_path.exists():
                meta = json.loads(meta_path.read_text())
                samples = meta["samples"]
                block_errors = meta["block_errors"]

            for t in NEW_T:
                if (snr, t) in have:
                    continue
                sp = OUTDIR / f"snr{t_tag}_static_t{t}{suffix}.csv"
                run_cached([PY, str(MCLASSIFIER / "architectures/static_pm_rank/static_path_metric_rank.py"),
                            "--trace", str(trace_csv), "--step", str(t), "--out", str(sp)], log, sp)

                out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                              "--pred", str(sp), "--topk", str(K)], log)
                bler = parse_eval(out)["topk_bler"]
                rows.append({"snr": snr, "t": t, "method": "static_pm_rank", "k": K,
                             "samples": samples, "block_errors": block_errors, "bler": bler})

            with open(results_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["snr", "t", "method", "k", "samples", "block_errors", "bler"])
                writer.writeheader()
                writer.writerows(sorted(rows, key=lambda r: (float(r["snr"]), r["method"], str(r["t"]))))

            log.write(f"SNR={snr} extra-t done\n")

    print(f"updated {results_path}")


if __name__ == "__main__":
    main()
