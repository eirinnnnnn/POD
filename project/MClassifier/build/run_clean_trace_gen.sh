#!/bin/bash
cd /home/eirin/Polar_eirin_20260107/project/MClassifier/build
for t in 4 8 16 32; do
  cps=$(python3 -c "print(','.join(str(c) for c in [4,8,16,32] if c<=$t))")
  echo "=== t=$t checkpoints=$cps (train) ===" >> clean_trace_gen_all.log
  ./mclass_trace_dataset -ini ebn0_sweep.ini --dataset clean_train_metric_12k.csv --checkpoints "$cps" --out clean_trace_12k_upto_t${t}.csv >> clean_trace_gen_all.log 2>&1
  echo "=== t=$t checkpoints=$cps (test) ===" >> clean_trace_gen_all.log
  ./mclass_trace_dataset -ini ebn0_sweep.ini --dataset clean_test_metric_3k.csv --checkpoints "$cps" --out clean_trace_3k_upto_t${t}.csv >> clean_trace_gen_all.log 2>&1
done
echo "ALL DONE" >> clean_trace_gen_all.log
