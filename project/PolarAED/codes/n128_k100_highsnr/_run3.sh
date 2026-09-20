#!/bin/bash
cd "$(dirname "$0")"
for f in AE16SCL2_UTL AE8SCL4_UTL AE2SCL16_UTL; do
  ( ../../build/PolarAED -ini ./${f}.ini > /dev/null 2>&1 ; echo "DONE $f" ) &
done
( ../../../OSD/build/OSD_performance -ini ./OSD1.ini > /dev/null 2>&1 ; echo "DONE OSD1" ) &
wait
echo "ALL DONE"
