#!/bin/bash
# Sweep P_IDX over SA_triC_49_1_5.py, launching each (P_IDX, seed) combo in its
# own detached tmux session so it survives after this script/terminal exits.
set -euo pipefail
cd "$(dirname "$0")"

P_IDX_LIST=(-2 -3 -4 -5)
SEEDS_PER_P=2

for p in "${P_IDX_LIST[@]}"; do
    pabs=${p#-}
    for s in $(seq 0 $((SEEDS_PER_P - 1))); do
        run_dir="triC/run_${pabs}_${s}"
        session="triC_p${pabs}_s${s}"
        seed=$(( pabs * 100 + s ))

        mkdir -p "$run_dir"

        if tmux has-session -t "$session" 2>/dev/null; then
            echo "skip ${session}: session already exists"
            continue
        fi

        tmux new-session -d -s "$session" \
            "uv run python -O SA_triC_49_1_5.py --run-dir ${run_dir} --P_IDX ${p} --seed ${seed} 2>&1 | tee ${run_dir}/stdout.log"
        echo "launched ${session} -> ${run_dir}  (P_IDX=${p}, seed=${seed})"
    done
done

echo
echo "active tmux sessions:"
tmux ls