#!/usr/bin/env bash
# Launch the eq-vs-shear permutation_src trial set.
#   ./run_eq_vs_shear.sh            # all 22, 4 at a time
#   ./run_eq_vs_shear.sh 8          # all 22, 8 at a time
#   ./run_eq_vs_shear.sh 4 SC_eq SCL4_eq   # just these
# Each job writes ./<ini basename>/log.txt (POD appends, so an existing
# log.txt from a previous run must be moved aside first).
set -u
cd "$(dirname "$0")"

JOBS="${1:-4}"
shift || true

DECODERS=(SC SCL4 SCL8 SCL32 SCL128 AED4SC AED8SC AED4SCL8 AED8SCL4 AED16SCL8 AED8SCL16)

if [ "$#" -gt 0 ]; then
    TASKS=("$@")
else
    TASKS=()
    for s in eq shear; do
        for d in "${DECODERS[@]}"; do TASKS+=("${d}_${s}"); done
    done
fi

mkdir -p run_logs
printf '%s\n' "${TASKS[@]}" | xargs -P "$JOBS" -I{} \
    bash -c './POD -ini {}.ini > run_logs/{}.stdout 2>&1; echo "{} exit=$?"'
