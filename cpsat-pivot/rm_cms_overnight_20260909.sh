#!/usr/bin/env bash
# Overnight: CMS stage-1 (RM-supercode realization) across codes and degree bounds.
# Each job: find P_0 with  C_b P_0 subseteq RM(d,6), enumerate a few, verify d,
# report pivot profile + Bhattacharyya Z-sum, write the P_0 matrices.
set -u
cd "$(dirname "$0")"
PY=.venv/bin/python
LOG=rm_cms_overnight.log
: > "$LOG"

run() {                     # gen  d  enum  cap_seconds  tag
  local gen=$1 d=$2 en=$3 cap=$4 tag=$5
  echo "######## $tag :  $(basename "$gen")  d<=$d  enumerate=$en  cap=${cap}s" | tee -a "$LOG"
  timeout "$cap" $PY -u rm_supercode_cms_20260909.py \
      --gen "$gen" --m 6 --d "$d" --enumerate "$en" --threads 6 \
      --max-seconds "$cap" --snr 3.0 \
      --out-prefix "cnc_work/P0_${tag}" >> "$LOG" 2>&1
  echo "---- $tag exit $? at $(date +%H:%M:%S) ----" | tee -a "$LOG"
}

# 1. eBCH[64,16] at the tight bound d=2  (dim RM(2,6)=22 >= 16; CP-SAT could not crack this)
run ../project/POD/eBCH_m6_t11.matrix 2 8 14400 ebch16_d2

# 2. eBCH[64,16] relaxed to d=3  (does the looser bound solve faster / give a better profile?)
run ../project/POD/eBCH_m6_t11.matrix 3 8 7200  ebch16_d3

# 3. BCH[64,36] at its floor d=3  (dim RM(3,6)=42 >= 36) and relaxed d=4
run ../data/BCH_36_64.matrix 3 4 14400 bch36_d3
run ../data/BCH_36_64.matrix 4 4 7200  bch36_d4

echo "######## ALL DONE  $(date) ########" | tee -a "$LOG"
