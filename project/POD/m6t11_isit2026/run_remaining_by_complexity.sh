#!/usr/bin/env bash
# Run the still-outstanding eq/shear jobs cheapest-first, so the trend across
# decoders shows up early instead of after the eff=128 runs drain.
#
# Tiers are by effective list size (PED order M x list size L), which is what
# dominates runtime.  Each tier finishes before the next starts.
#
# Deliberately excludes the runs that are already complete or already in
# flight from the first launch: SC/SCL4/SCL8/AED4SC/AED8SC on the eq side
# (done) and SCL32_eq/SCL128_eq/AED4SCL8_eq/AED8SCL4_eq (running).
set -u
cd "$(dirname "$0")"

JOBS="${1:-4}"

TIER1=(SC_shear SCL4_shear AED4SC_shear SCL8_shear AED8SC_shear)          # eff 1-8
TIER2=(SCL32_shear AED4SCL8_shear AED8SCL4_shear)                          # eff 32
TIER3=(SCL128_shear AED16SCL8_shear AED8SCL16_shear AED16SCL8_eq AED8SCL16_eq)  # eff 128

mkdir -p run_logs

run_tier() {
    local name="$1"; shift
    echo "=== tier $name: $* ==="
    printf '%s\n' "$@" | xargs -P "$JOBS" -I{} \
        bash -c './POD -ini {}.ini > run_logs/{}.stdout 2>&1; echo "{} exit=$?"'
    echo "=== tier $name done ==="
}

run_tier 1 "${TIER1[@]}"
run_tier 2 "${TIER2[@]}"
run_tier 3 "${TIER3[@]}"
