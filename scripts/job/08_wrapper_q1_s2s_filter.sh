for ens in e01 e02 e03; do
    for wave in mjo er td kelvin mrg; do
        for var in q1 q2 qr_eddy; do

            # var=q1 
            # ens=e01
            # wave=mjo
            qsub -v ens="$ens",wave="$wave",var="$var" /home/563/fm6730/localrepo/wave_vertical_heating/scripts/job/08b_ncl_q1_s2s_filter.qsub

        done
    done
done