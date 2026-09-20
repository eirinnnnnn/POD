#!/usr/bin/env bash
# Wait for the current POD sweep to drain, then run the 4.75-6.00 dB extension.
# Results land in <name>_ext/log.txt and are appended to <name>/log.txt afterwards.
cd "$(dirname "$0")"
while [ "$(pgrep -xc POD)" -gt 0 ]; do sleep 60; done
xargs -P 6 -I{} bash -c './POD -ini {}.ini > run_logs/{}.stdout 2>&1; echo "{} exit=$?"' < tasks_ext.txt
echo "=== extension done; merging into base logs ==="
while read -r n; do
  base="${n%_ext}"
  if [ -f "$n/log.txt" ] && [ -f "$base/log.txt" ]; then
      if ! grep -q "SNR =  4.75" "$base/log.txt"; then
          cat "$n/log.txt" >> "$base/log.txt"; echo "  merged $n -> $base"
      else
          echo "  SKIP $base (already has 4.75+)"
      fi
  fi
done < tasks_ext.txt
