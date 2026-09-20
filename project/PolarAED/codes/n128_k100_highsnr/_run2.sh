#!/bin/bash
# Pass 2: extend the sweep to 5.5 and 6.0 dB. Re-uses the SAME ini names so
# main.cpp appends to the existing <run>/log.txt.
cd "$(dirname "$0")"
while pgrep -f 'PolarAED -ini' > /dev/null; do sleep 10; done
python3 _mkini.py 5.5 6.1 800000 100 ""
for f in SC SCL32 AE32SC_LTA AE32SC_UTL AE32SC_PU AE32SC_RND AE4SCL8_UTL; do
  ( ../../build/PolarAED -ini ./${f}.ini > ${f}.hi.stdout 2>&1 ; echo "DONE2 $f" ) &
done
wait
echo "ALL DONE2"
