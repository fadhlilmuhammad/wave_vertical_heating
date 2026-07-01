# for var in ua va wa hus ta; do 
for var in wa; do 
    for ens in e01 e02 e03; do 
        qsub -v var="$var",ens="$ens" /home/563/fm6730/localrepo/wave_vertical_heating/scripts/job/02_cdo_get3d_access.qsub
    done 
done