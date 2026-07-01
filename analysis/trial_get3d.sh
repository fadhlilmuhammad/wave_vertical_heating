#!/bin/bash

year=1981
var=t
folder_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/era5/daily/t
file=${var}_era5_oper_pl_merge_1deg_daily_${year}.nc

folder_out_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/temp
fileout=${var}_era5_oper_pl_merge_1deg_daily_eq_${year}.nc

cdo sellonlatbox,0,360,-25,25 $folder_path/$file $folder_out_path/$fileout