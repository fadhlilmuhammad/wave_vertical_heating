for var in ua va wa hus ta; do 
    for ens in e01 e02 e03; do
        for f in /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/chunks/chunk_*; do
            qsub -v var="$var",ens="$ens",input_file="$f" \
                /home/563/fm6730/localrepo/wave_vertical_heating/scripts/job/02_cdo_get3d_access.qsub 
        done

    done
done