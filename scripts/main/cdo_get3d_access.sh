#!/bin/bash

yyyymmdd=$1
var=$2
ens=$3

folder_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/access-s2_hindcast_raw/${var}/daily/${ens}
file=da_${var}_${yyyymmdd}_${ens}.nc

case "$var" in
    ua) var="u" ;;
    va) var="v" ;;
    hus) var="q" ;;
    ta) var="t" ;;
    wa) var="w" ;;
esac

grid=/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/grid.txt

folder_out_path=/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/$var/access-s2
fileout=da_${var}_${yyyymmdd}_${ens}_2-5deg.nc

# rm $folder_path/$file $folder_out_path/$fileout
mkdir -p $folder_out_path
cdo -remapcon,$grid $folder_path/$file $folder_out_path/$fileout
