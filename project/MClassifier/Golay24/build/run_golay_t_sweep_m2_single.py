#!/usr/bin/env python3
"""Golay24 t-sweep trial: fixed pruning width m=2, checkpoint t in
{4,8,16,20}, static-pm-rank vs trace-learned (one model per t, each trained
at SNR=2.5dB -- see golay_train_t{T}_12k.csv / golay_trace_t{T}_weights.txt),
over the same 8-point SNR grid as the fixed-t=8 trial, target_errors=1000
with sample caps sized to actually reach it (same corrected caps as that
trial). full_ped (M=16) is checkpoint-independent (t only affects which
mid-decode snapshot is used for RANKING, never the final decode/selection),
so it's reused directly from golay_trace_vs_static_results.txt rather than
recomputed. Also includes pure-SCL(L=8) and pure-SC(L=1) no-AED baselines,
same as the fixed-t=8 trial."""
import csv
import subprocess
from pathlib import Path

BUILD = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/Golay24/build")
SIM = BUILD / "golay_mclass_bler_sim"
INI = BUILD / "golay_sweep_m16_L2.ini"
SCL8_INI = BUILD / "golay_scl8.ini"
SC_INI = BUILD / "golay_sc_pure.ini"
OUTDIR = BUILD / "sweep"

SNR_LIST = [2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50]
SAMPLE_CAP = {2.00: 30_000, 2.25: 45_000, 2.50: 55_000, 2.75: 85_000,
              3.00: 120_000, 3.50: 300_000, 4.00: 800_000, 4.50: 2_000_000}
TARGET_ERRORS = 1000
T_LIST = [4, 8, 16, 20]
M = 2
RESUME_EVERY = 20_000

RESULTS = OUTDIR / "golay_t_sweep_m2_single_results.txt"
SOURCE_RESULTS = OUTDIR / "golay_trace_vs_static_results.txt"  # reused for full_ped
FIELDS = ["snr", "method", "t", "m", "samples", "errors", "bler"]


def tag(snr):
    return f"{snr:.2f}".replace(".", "p")


def seed_pair(snr_idx, salt):
    return 760001 + salt * 100 + snr_idx, 660001 + salt * 100 + snr_idx


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
    key = lambda r: (f"{float(r['snr']):.2f}", r["method"], str(r["t"]), str(int(float(r["m"]))))
    have = {key(r) for r in existing}
    for r in rows:
        if key(r) not in have:
            existing.append(r)
    with open(RESULTS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(existing)


def import_full_ped():
    """Reuse full_ped rows from the fixed-t=8 trial's results (checkpoint-independent)."""
    if not SOURCE_RESULTS.exists():
        print("warning: source results for full_ped not found, skipping import")
        return
    with open(SOURCE_RESULTS, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["method"] == "full_ped"]
    out_rows = [{"snr": r["snr"], "method": "full_ped", "t": "", "m": 16,
                 "samples": r["samples"], "errors": r["errors"], "bler": r["bler"]} for r in rows]
    append_rows(out_rows)
    print(f"imported {len(out_rows)} full_ped rows from {SOURCE_RESULTS.name}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    import_full_ped()

    # pure-SCL(L=8) and pure-SC(L=1) baselines, t-independent
    for method, ini, salt in [("pure_scl8", SCL8_INI, 900), ("pure_sc", SC_INI, 901)]:
        for i, snr in enumerate(SNR_LIST):
            t_tag = tag(snr)
            cap = SAMPLE_CAP[snr]
            cs, ms = seed_pair(i, salt)
            log = OUTDIR / f"golay_{method}_m2trial_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(ini), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--checkpoint", "8",  # required arg but irrelevant when use_AED=false (single candidate)
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--m-list", "1",
                "--monitor-every", "50000",
            ], f"golay_{method}_m2trial_snr{t_tag}")
            append_rows([{"snr": snr, "method": method, "t": "", "m": 1,
                          "samples": vals["samples"], "errors": vals["full_ped_errors"],
                          "bler": vals["full_ped_bler"]}])
            print(f"SNR={snr} {method} done (cached={cached}, samples={vals['samples']:.0f})")

    for T in T_LIST:
        model = BUILD / f"golay_trace_t{T}_weights.txt"
        for i, snr in enumerate(SNR_LIST):
            t_tag = tag(snr)
            cap = SAMPLE_CAP[snr]
            cs, ms = seed_pair(i, T)
            log = OUTDIR / f"golay_t_sweep_m2_t{T}_snr{t_tag}.log"
            vals, cached = run_sim(log, [
                "-ini", str(INI), "--snr", str(snr),
                "--samples", str(cap), "--target-errors", str(TARGET_ERRORS),
                "--checkpoint", str(T),
                "--channel-seed", str(cs), "--message-seed", str(ms),
                "--m-list", str(M),
                "--model", str(model),
                "--monitor-every", "50000",
            ], f"golay_t_sweep_m2_t{T}_snr{t_tag}")
            samples = vals["samples"]
            rows = [
                {"snr": snr, "method": "static_pm_rank", "t": T, "m": M,
                 "samples": samples, "errors": vals[f"static_errors[m={M}]"],
                 "bler": vals[f"static_bler[m={M}]"]},
                {"snr": snr, "method": "trace_learned", "t": T, "m": M,
                 "samples": samples, "errors": vals[f"trace_errors[m={M}]"],
                 "bler": vals[f"trace_bler[m={M}]"]},
            ]
            append_rows(rows)
            print(f"SNR={snr} t={T} done (cached={cached}, samples={samples:.0f})")

    print(f"updated {RESULTS}")


if __name__ == "__main__":
    main()
