#!/usr/bin/env bash
# cnc_solve  --  cube-and-conquer for one plain-DIMACS CNF.
#
#   ./cnc_solve_20260908.sh  FILE.cnf  [PAR]  [march args...]
#
# march_cu splits FILE into cubes; PAR (default 8) iglucose workers each take a
# round-robin slice of the cubes as an incremental (p inccnf) formula.  Prints
#   s SATISFIABLE      -- one cube satisfiable  (+ the winning worker log)
#   s UNSATISFIABLE    -- every cube refuted
#   s UNKNOWN          -- a worker hit its per-run limit without deciding
#
set -u
CNF=$1
PAR=${2:-8}
shift $(( $# >= 2 ? 2 : 1 ))
MARCH_ARGS=("$@")

HERE=$(cd "$(dirname "$0")" && pwd)
DIR=$HERE/tools/CnC
MARCH=$DIR/march_cu/march_cu
IGLU=$DIR/iglucose/core/iglucose
WORK=$(mktemp -d /tmp/cnc.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

t0=$(date +%s)
echo "c [cnc] cubing $CNF  ${MARCH_ARGS[*]}"
"$MARCH" "$CNF" -o "$WORK/cubes" "${MARCH_ARGS[@]}" 2>&1 \
  | grep -E "number of cubes|cutoff|time =" | sed 's/^/c [march] /'
NCUBES=$(grep -c '^a ' "$WORK/cubes" 2>/dev/null || echo 0)
echo "c [cnc] $NCUBES cubes -> $PAR workers"
if [ "$NCUBES" -eq 0 ]; then echo "s UNKNOWN"; exit 2; fi

grep -v '^c' "$CNF" > "$WORK/base.cnf"

pids=()
for ((w=0; w<PAR; w++)); do
  {
    echo "p inccnf"
    cat "$WORK/base.cnf"
    awk -v p="$PAR" -v w="$w" 'NR % p == w' "$WORK/cubes"
  } > "$WORK/f$w.icnf"
  "$IGLU" "$WORK/f$w.icnf" "$WORK/out$w.txt" -verb=0 &
  pids+=($!)
done

RESULT=UNKNOWN
while :; do
  running=0
  for i in "${!pids[@]}"; do kill -0 "${pids[$i]}" 2>/dev/null && running=$((running+1)); done
  if grep -lq '^SAT' "$WORK"/out*.txt 2>/dev/null; then
    RESULT=SATISFIABLE
    for p in "${pids[@]}"; do kill "$p" 2>/dev/null; done
    break
  fi
  nunsat=$(grep -l '^UNSAT' "$WORK"/out*.txt 2>/dev/null | wc -l)
  if [ "$nunsat" -eq "$PAR" ]; then RESULT=UNSATISFIABLE; break; fi
  if [ "$running" -eq 0 ]; then
    # all finished, none SAT, not all cleanly UNSAT -> at least one INDETERMINATE
    if [ "$nunsat" -eq "$PAR" ]; then RESULT=UNSATISFIABLE; else RESULT=UNKNOWN; fi
    break
  fi
  sleep 2
done
wait 2>/dev/null

t1=$(date +%s)
echo "c [cnc] workers: $(grep -h -oE '^(SAT|UNSAT|INDETERMINATE)' "$WORK"/out*.txt 2>/dev/null | sort | uniq -c | tr '\n' ' ')"
echo "c [cnc] elapsed $((t1-t0))s"
echo "s $RESULT"
[ "$RESULT" = SATISFIABLE ] && exit 10
[ "$RESULT" = UNSATISFIABLE ] && exit 20
exit 2
