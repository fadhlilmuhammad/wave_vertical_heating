#!/bin/bash

# Patching the error due to bugs when creating data at YYYY0101. Bug found on the NCL code and now is corrected.

set -euo pipefail

JOBDIR="/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob"
OUTFILE="${JOBDIR}/dates_s2sfilter_patch.txt"
> "$OUTFILE"

start="1982-01-01"
end="2017-01-01"

d="$start"
while [[ "$(date -d "$d" +%Y%m%d)" -le "$(date -d "$end" +%Y%m%d)" ]]; do
    day_min0=$(date -d "$d" +%Y%m%d)
    day_min1=$(date -d "$d -1 day" +%Y%m%d)
    day_min2=$(date -d "$d -2 day" +%Y%m%d)

    echo "${day_min0} ${day_min1} ${day_min2}" >> "$OUTFILE"

    d=$(date -d "$d +1 year" +%Y-%m-01)
done

echo "Wrote $(wc -l < "$OUTFILE") dates to $OUTFILE"
chmod +x $OUTFILE