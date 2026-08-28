#!/bin/bash
VARS=(
    # DONE
    # "vorticity"
    # "vorticity_budget"
    # "q2"
    # "q1q2_yanai"

    "adv_q850"
    "conv_q850"

    # NOT DONE, waiting for reviewer
    # "u850"
    # "v850"
    # "rlut"
    # "q850"
)

BASE_DIR="/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/access-s2"  # <-- set this to the directory containing your ens_member subfolders


# ENS_MEMBER=(
#     "e01"
#     "e02"
#     "e03"
#     "e01_tlag_long"
#     "e02_tlag_long"
#     "e03_tlag_long"
#     "e01_prep"
#     "e02_prep"
#     "e03_prep"
#     "e01_anom"
#     "e02_anom"
#     "e03_anom"
#     "ensmean_anom"
#     "ensmean_prep"
# )

# VARS=(
#     "adv_q"
#     "conv_q"
#     "q1q2_yanai"
# )

# scan all subfolders and use their names as ENS_MEMBER values
ENS_MEMBER=()
for var in "${VARS[@]}"; do
    for dir in "$BASE_DIR"/$var/*; do
        ENS_MEMBER+=("$(basename "$dir")")
    done
done

echo "Found ensemble members: ${ENS_MEMBER[@]}"


for var in "${VARS[@]}"; do
    for ens_member in "${ENS_MEMBER[@]}"; do
        echo "${ens_member} ${var}"

        # integrated variables
        # qsub -v ENS_MEMBER="$ens_member",VAR="$var" /home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/compress_tar_formassdata.pbs
        
        # non-integrated variables
        qsub -v ENS_MEMBER="$ens_member",VAR="$var" /home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/compress_tar_formassdata_nonintegrated.pbs
    done
done