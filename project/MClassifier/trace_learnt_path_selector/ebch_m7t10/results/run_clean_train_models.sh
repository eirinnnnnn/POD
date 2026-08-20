#!/bin/bash
cd /home/eirin/Polar_eirin_20260107/project/MClassifier/build
for t in 4 8 16 32; do
  echo "=== training t=$t ===" >> clean_train_models.log
  python3 ../architectures/trace/train_mclass_trace_ranker.py \
    --trace clean_trace_12k_upto_t${t}.csv \
    --out clean_trace_upto_t${t}_basin_h64_e20.pt \
    --target basin --hidden 64 --epochs 20 >> clean_train_models.log 2>&1
  echo "exit=$?" >> clean_train_models.log
done
echo "ALL TRAINING DONE" >> clean_train_models.log
