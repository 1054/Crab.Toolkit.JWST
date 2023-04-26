#!/bin/bash
# 

usage() {
    echo "Usage: "
    echo "  go-run-screen-jobs-in-parallel.bash go_qsub_script.bash"
    echo "Notes: "
    echo "  This runs the qsub script with screen."
    echo "  The qsub "
}

if [[ $# -lt 1 ]]; then
    usage
    exit
fi

script_file=$1 # ./go_qsub_processing_jwst_imaging_20230422_00h04m35s.bash
chmod +x $script_file

timestamp=$(date +"%Y%m%d_%Hh%Mm%Ss")
process_prefix=batch_go_run_qsub_script
max_concurrent=9
sleep_level=1

for (( i=0; i<9; i++ )); do
    export PBS_ARRAYID=$((i+1))
    sleep 0.5
    sleep_level=1
    while [[ $(screen -ls | grep ".$process_prefix." | wc -l) -gt $max_concurrent ]]; do
        sleep $(awk "BEGIN {print 2**($sleep_level);}")
        if [[ $sleep_level -lt 5 ]]; then
            sleep_level=$((sleep_level+1))
        fi
    done
    echo "logfile log_$process_prefix.$PBS_ARRAYID.txt" > log_$process_prefix.$PBS_ARRAYID.conf
    echo "logfile flush 1" >> log_$process_prefix.$PBS_ARRAYID.conf
    echo "log on" >> log_$process_prefix.$PBS_ARRAYID.conf
    echo "screen -d -m -S \"$process_prefix.$PBS_ARRAYID\" -L -c \"log_$process_prefix.$PBS_ARRAYID.conf\" $script_file"
    screen -d -m -S "$process_prefix.$PBS_ARRAYID" -L -c "log_$process_prefix.$PBS_ARRAYID.conf" $script_file
    sleep 1.5
done

