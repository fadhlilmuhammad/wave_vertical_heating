#!/bin/bash

outfile="/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/year_foraccess.txt"

> "$outfile"

current="1981-02-01"

while [[ "$current" < "2019-01-01" ]]; do

    d0=$(date -d "$current" +%Y%m%d)
    d1=$(date -d "$current -1 day" +%Y%m%d)
    d2=$(date -d "$current -2 day" +%Y%m%d)

    echo "$d2 $var $e" >> "$outfile"
    echo "$d1 $var $e" >> "$outfile"
    echo "$d0 $var $e" >> "$outfile"

    current=$(date -d "$current +1 month" +%Y-%m-%d)

done

chmod +x "$outfile"
