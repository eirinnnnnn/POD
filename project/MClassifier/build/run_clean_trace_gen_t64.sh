#!/bin/bash
cd /home/eirin/Polar_eirin_20260107/project/MClassifier/build
echo "=== t=64 checkpoints=4,8,16,32,64 (train) ===" >> clean_trace_gen_t64.log
./mclass_trace_dataset -ini ebn0_sweep.ini --dataset clean_train_metric_12k.csv --checkpoints "4,8,16,32,64" --out clean_trace_12k_upto_t64.csv >> clean_trace_gen_t64.log 2>&1
echo "=== t=64 checkpoints=4,8,16,32,64 (test) ===" >> clean_trace_gen_t64.log
./mclass_trace_dataset -ini ebn0_sweep.ini --dataset clean_test_metric_3k.csv --checkpoints "4,8,16,32,64" --out clean_trace_3k_upto_t64.csv >> clean_trace_gen_t64.log 2>&1
echo "ALL DONE" >> clean_trace_gen_t64.log
