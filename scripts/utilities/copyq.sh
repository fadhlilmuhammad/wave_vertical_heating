#!/bin/bash
#PBS -N copyq
#PBS -P v46
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l wd
#PBS -l mem=150GB
#PBS -l storage=massdata/v46+gdata/v46+scratch/v46
#PBS -l walltime=10:00:00
#PBS -j oe

set -e

v1=adv_q

SRC="/g/data/v46/fm6730/data/access_s2/integrated/${v1}"
TAR="/scratch/v46/fm6730/data/access_s2/access_s2_${v1}_raw.tar.gz"

echo "Starting ${v1}: $(date)"
echo "Source: $SRC"
echo "Tar: $TAR"

# Send to MASSDATA
mdss put "$TAR" fm6730/data/access_s2/

# echo "MASSDATA transfer completed: $(date)"

# Remove temporary tar
# rm -f "$TAR"

echo "Finished ${v1}: $(date)"