# for var in u v w q t; do 
for var in q1; do 
    qsub -v var="$var" \
        /home/563/fm6730/localrepo/wave_vertical_heating/scripts/job/03_create_doyclim.qsub 
done