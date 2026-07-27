#!/bin/bash

for year in $(seq 1979 2020); do
    echo "$year"
done > /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/year_input.txt
chmod +x /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/year_input.txt
