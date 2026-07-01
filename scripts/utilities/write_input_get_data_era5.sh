#!/bin/bash

for var in t q u v w; do
    for year in $(seq 1979 2020); do
        echo "$year $var"
    done > /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/getdata_input_${var}.txt
    chmod +x /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/getdata_input_${var}.txt
done