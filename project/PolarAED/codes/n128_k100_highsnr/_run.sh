#!/bin/bash
cd "$(dirname "$0")"
for f in SC SCL32 AE32SC_LTA AE32SC_UTL AE32SC_PU AE32SC_RND AE4SCL8_UTL; do
  ( ../../build/PolarAED -ini ./${f}.ini > ${f}.stdout 2>&1 ; echo "DONE $f" ) &
done
wait
echo "ALL DONE"
