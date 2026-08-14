#!/usr/bin/env python3
"""Genuine full-PED baseline at M=8 (aed_L=8, same automorphism ordering
file truncated to its first 8 entries -- a real, smaller-ensemble decode,
not a top-8-of-64 selection). One row per SNR point, appended into the
main results.txt as method="full_ped_m8", t="all", k=8."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_sweep import BUILD, OUTDIR, SNR_LIST, TARGET_ERRORS, run, run_cached, parse_eval  # noqa: E402

INI_M8 = BUILD / "ebn0_sweep_m8.ini"
SAMPLE_CAP = 100_000


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def write_identity_pred(path, n_samples, m):
    with open(path, "w", encoding="utf-8") as f:
        f.write("sample" + "".join(f",idx{i}" for i in range(m)) + "\n")
        row = "," + ",".join(str(i) for i in range(m))
        for s in range(n_samples):
            f.write(f"{s}{row}\n")


def main():
    results_path = OUTDIR / "results.txt"
    log_path = OUTDIR / "sweep_m8.log"

    with open(results_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    have = {float(r["snr"]) for r in rows if r["method"] == "full_ped_m8"}

    with open(log_path, "a", encoding="utf-8") as log:
        for i, snr in enumerate(SNR_LIST):
            if snr in have:
                log.write(f"\n===== M8 SNR={snr} already complete, skipping =====\n")
                log.flush()
                continue

            t_tag = tag(snr)
            log.write(f"\n===== M8 SNR={snr} ({t_tag}) =====\n")
            channel_seed = 950001 + i
            message_seed = 850001 + i

            base_csv = OUTDIR / f"snr{t_tag}_metric_m8.csv"
            meta_path = Path(str(base_csv) + ".meta.json")
            run_cached([
                str(BUILD / "mclass_dataset"), "-ini", str(INI_M8),
                "--snr", str(snr),
                "--target-errors", str(TARGET_ERRORS),
                "--samples", str(SAMPLE_CAP),
                "--channel-seed", str(channel_seed),
                "--message-seed", str(message_seed),
                "--out", str(base_csv),
            ], log, base_csv, side_suffixes=[".meta.json"])

            import json
            meta = json.loads(meta_path.read_text())
            samples = meta["samples"]
            block_errors = meta["block_errors"]

            pred_csv = OUTDIR / f"snr{t_tag}_identity_m8.csv"
            if not pred_csv.exists():
                write_identity_pred(pred_csv, samples, 8)

            out, _ = run([str(BUILD / "mclass_eval"), "--dataset", str(base_csv),
                          "--pred", str(pred_csv), "--topk", "8"], log)
            bler = parse_eval(out)["topk_bler"]

            rows.append({"snr": snr, "t": "all", "method": "full_ped_m8", "k": 8,
                         "samples": samples, "block_errors": block_errors, "bler": bler})
            with open(results_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["snr", "t", "method", "k", "samples", "block_errors", "bler"])
                writer.writeheader()
                writer.writerows(sorted(rows, key=lambda r: (float(r["snr"]), r["method"], str(r["t"]))))

            log.write(f"M8 SNR={snr} done: samples={samples} block_errors={block_errors}\n")

    print(f"updated {results_path}")


if __name__ == "__main__":
    main()
