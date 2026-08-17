#!/bin/bash
#PBS -N tar_to_massdata
#PBS -P v46
#PBS -q normal
#PBS -l ncpus=4
#PBS -l wd
#PBS -l mem=16GB
#PBS -l storage=massdata/v46+gdata/v46+scratch/v46
#PBS -l walltime=12:00:00
#PBS -j oe

set -e

v1=adv_q
v2=conv_q
v3=q1q2_yanai

SRC_1="/g/data/v46/fm6730/data/access_s2/integrated/$v1"
SRC_2="/g/data/v46/fm6730/data/access_s2/integrated/$v2"
SRC_3="/g/data/v46/fm6730/data/access_s2/integrated/$v3"
# SRC_4="/g/data/v46/fm6730/data/access_s2/integrated/q2"
# SRC_5="/g/data/v46/fm6730/data/access_s2/integrated/vorticity"
# SRC_6="/g/data/v46/fm6730/data/access_s2/integrated/vorticity_budget"

TAR_1="/scratch/v46/fm6730/data/access_s2/access_s2_${v1}_raw.tar.gz"
TAR_2="/scratch/v46/fm6730/data/access_s2/access_s2_${v2}_raw.tar.gz"
TAR_3="/scratch/v46/fm6730/data/access_s2/access_s2_${v3}_raw.tar.gz"
# DEST="/massdata/v46/fm6730.tar"

echo "Starting: $(date)"

# Create tar archive
tar -cvzf "$TAR_1" "$SRC_1"
tar -cvzf "$TAR_1" "$SRC_2"
tar -cvzf "$TAR_1" "$SRC_3"

echo "Tar completed: $(date)"

# Send to MASSDATA
mdss put $TAR_1 fm6730/data/access_s2/
mdss put $TAR_2 fm6730/data/access_s2/
mdss put $TAR_3 fm6730/data/access_s2/

echo "Copy completed: $(date)"

# Remove temporary tar file
rm -f "$TAR_1"
rm -f "$TAR_2"
rm -f "$TAR_3"

echo "Finished: $(date)"