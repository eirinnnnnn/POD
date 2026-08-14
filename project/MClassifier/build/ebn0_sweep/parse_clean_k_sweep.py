#!/usr/bin/env python3
import csv
import re
from pathlib import Path

LOG = Path("/home/eirin/Polar_eirin_20260107/project/MClassifier/build/clean_k_sweep.log")
OUT = Path(__file__).resolve().parent / "clean_k_sweep_results.txt"

rows = []
t = None
samples = None
full_ped = None
for line in LOG.read_text().splitlines():
    m = re.match(r"=== t=(\d+) ===", line)
    if m:
        t = int(m.group(1))
        continue
    if line.startswith("samples,"):
        samples = int(line.split(",")[1])
        continue
    if line.startswith("full_ped_bler,"):
        full_ped = float(line.split(",")[1])
        rows.append({"t": t, "method": "full_ped", "k": 64, "samples": samples, "bler": full_ped})
        continue
    m = re.match(r"(static|trace)_bler\[k=(\d+)\],([\d.eE+-]+)", line)
    if m:
        method = "static_pm_rank" if m.group(1) == "static" else "trace_learned"
        rows.append({"t": t, "method": method, "k": int(m.group(2)), "samples": samples, "bler": float(m.group(3))})

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["t", "method", "k", "samples", "bler"])
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT} ({len(rows)} rows)")
