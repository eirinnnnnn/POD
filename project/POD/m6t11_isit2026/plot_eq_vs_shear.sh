#!/usr/bin/env bash
# Render the eq / shear / overlay figures for the permutation_src trial set.
# Folder names are listed explicitly on purpose: a *_shear glob would also
# pick up the unrelated schi_shear / SChi_shear runs.
#
# Safe to run mid-sweep: each figure is attempted independently and a figure
# with no finished runs yet is skipped rather than aborting the rest.
set -u
cd "$(dirname "$0")"

# eff=128 (SCL128, AED16SCL8, AED8SCL16) is deliberately left out; eff=16
# covers every (M, L) split with M*L = 16.
DECODERS=(
    SC                                  # eff 1
    SCL4 AED4SC                         # eff 4
    SCL8 AED8SC                         # eff 8
    SCL16 AED16SC AED2SCL8 AED4SCL4 AED8SCL2   # eff 16
    SCL32 AED4SCL8 AED8SCL4             # eff 32
)

EQ=(); SHEAR=(); GL=()
for d in "${DECODERS[@]}"; do EQ+=("${d}_eq"); SHEAR+=("${d}_shear"); GL+=("${d}_gl"); done

# ML decoding does not depend on the polar permutation, so MLD is one shared
# reference curve on every figure rather than an eq/shear pair.
# mld_merged/MLD pools this project's MLD log with AED/m6t11_isit2026's
# MLD_16_64 run.  Those use different generator matrices, but the two span
# permutation-equivalent codes (identical weight enumerator, d=24), so ML
# performance is the same and the error/iteration counts can be pooled.
REF=(mld_merged/MLD)

render() {  # render <savepath> <title> <folder>...
    local out="$1" title="$2"; shift 2
    echo "=== $out"
    python3 plot_logs.py "$@" --title "$title" --savepath "$out" \
        || echo "[skip] $out: no finished runs yet"
}

render ./eBCH_m6_t11_perm_eq.png \
    'eBCH(64,16), permutation: equation' "${EQ[@]}" "${REF[@]}"

render ./eBCH_m6_t11_perm_shear.png \
    'eBCH(64,16), permutation: shear $z_1z_3$' "${SHEAR[@]}" "${REF[@]}"

render ./eBCH_m6_t11_perm_gl.png \
    'eBCH(64,16), permutation: GL(6,2) optimum' "${GL[@]}" "${REF[@]}"

render ./eBCH_m6_t11_perm_all.png \
    'eBCH(64,16): equation (solid), shear (dashed), GL optimum (dotted)' "${EQ[@]}" "${SHEAR[@]}" "${GL[@]}" "${REF[@]}"
