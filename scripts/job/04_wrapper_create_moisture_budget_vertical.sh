for f in /home/563/fm6730/localrepo/wave_vertical_heating/scripts/input_forjob/chunks/chunk_*; do
    for ens in e01 e02 e03; do 
    # for ens in e01; do 
        qsub -v ens="$ens",input_file="$f" \
            /home/563/fm6730/localrepo/wave_vertical_heating/scripts/job/04a_create_moisture_budget_vertical.qsub 
    done
done