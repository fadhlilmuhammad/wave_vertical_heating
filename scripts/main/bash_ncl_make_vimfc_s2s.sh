
module use /g/data/xp65/public/modules/
module load conda/analysis3
module load ncl

ncl /home/563/fm6730/localrepo/wave_vertical_heating/scripts/main/ncl_make_vimfc_s2s.ncl
python3 /home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/compression_s2s_mfc.py