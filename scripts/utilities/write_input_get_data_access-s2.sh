#!/bin/bash

# for var in t q u v w; do
#     for e in e01 e02 e03; do

        outfile="/home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/year.txt"

        > "$outfile"

        current="1982-02-01"
        end="2018-12-01"

        while [[ "$current" < "2019-01-01" ]]; do
            yyyymmdd=$(date -d "$current" +%Y%m%d)

            echo "$yyyymmdd $var $e" >> "$outfile"

            current=$(date -d "$current +1 month" +%Y-%m-%d)
        done

        chmod +x "$outfile"

#     done
# done