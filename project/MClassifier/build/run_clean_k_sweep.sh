#!/bin/bash
cd /home/eirin/Polar_eirin_20260107/project/MClassifier/build
for t in 4 8 16 32; do
  echo "=== t=$t ===" >> clean_k_sweep.log
  ./mclass_bler_sim -ini ebn0_sweep.ini --snr 2.0 --samples 50000 --target-errors 1000 \
    --stop-metric static_trace_min \
    --channel-seed $((930000+t)) --message-seed $((830000+t)) \
    --checkpoint $t --k-list 1,2,4,8,16,32,64 \
    --model clean_trace_upto_t${t}_weights.txt \
    --resume-file clean_ksweep_t${t}.state --resume-every 20000 \
    >> clean_k_sweep.log 2>&1
  echo "exit=$?" >> clean_k_sweep.log
done
echo "K-SWEEP DONE" >> clean_k_sweep.log
