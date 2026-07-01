#!/bin/bash

year=$1
var=$2
folder_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/era5/daily/${var}
file=${var}_era5_oper_pl_merge_1deg_daily_${year}.nc

folder_out_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/$var/era5
fileout=${var}_era5_oper_pl_merge_2.5deg_daily_eq_${year}.nc

mkdir -p $folder_out_path
cdo remapcon,r144x73 $folder_path/$file $folder_out_path/$fileout
