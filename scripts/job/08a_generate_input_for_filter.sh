#!/bin/bash
# generate_datelist.sh
#
# Builds dates.txt: one line per month, with
#   day_min0 = first of month
#   day_min1 = day_min0 - 1 day
#   day_min2 = day_min0 - 2 days
# This file is shared by every (ens, wave, var) job — each job loops
# through all 455 lines sequentially inside run_wave_filter.sh.
#
# Run this once before submitting the PBS job:
#   bash generate_datelist.sh

set -euo pipefail

JOBDIR="/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob"
OUTFILE="${JOBDIR}/dates_s2sfilter.txt"
> "$OUTFILE"

start="1981-02-01"
end="2018-12-01"

d="$start"
while [[ "$(date -d "$d" +%Y%m%d)" -le "$(date -d "$end" +%Y%m%d)" ]]; do
    day_min0=$(date -d "$d" +%Y%m%d)
    day_min1=$(date -d "$d -1 day" +%Y%m%d)
    day_min2=$(date -d "$d -2 day" +%Y%m%d)

    echo "${day_min0} ${day_min1} ${day_min2}" >> "$OUTFILE"

    d=$(date -d "$d +1 month" +%Y-%m-01)
done

echo "Wrote $(wc -l < "$OUTFILE") dates to $OUTFILE"
chmod +x $OUTFILE