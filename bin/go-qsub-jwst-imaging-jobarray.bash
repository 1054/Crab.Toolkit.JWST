#!/bin/bash
# 
dataset_names=($(ls -1d jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_*))
if [[ ${#dataset_names[@]} -eq 0 ]]; then
    echo "No dataset dir is found: jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_*"
    echo "under current directory: $(pwd -P)"
    exit 255
fi
crds_context="$CRDS_CONTEXT" # "jwst_0986.pmap"
conda_env="jwstpmap0995" # "base"
concurrent=10
ncpu=1
mem="40gb" # maximum-cores 48 processes will need 68GB
timestamp=$(date +"%Y%m%d_%Hh%Mm%Ss")
currentdir=$(pwd -P)
goscript="go_qsub_processing_jwst_imaging_${timestamp}.bash"
echo "#!/bin/bash" > $goscript
echo "#PBS -N JW${timestamp}" >> $goscript
echo "#PBS -l nodes=1:ppn=${ncpu},mem=${mem},walltime=24:00:00" >> $goscript
echo "#PBS -d ${currentdir}/" >> $goscript
echo "#PBS -o log_processing_jwst_imaging_${timestamp}" >> $goscript
#echo "#PBS -e log_processing_jwst_imaging_${timestamp}_\${PBS_ARRAYID}.err" >> $goscript
echo "#PBS -j oe" >> $goscript # join stdout and stderr
#echo "#PBS -k oe" >> $goscript # keep stdout and stderr on running hostname
#echo "#PBS -m abe" >> $goscript # send email notifications for all, begin, end
echo "#PBS -m n" >> $goscript # do not send emails -- not working
echo "#PBS -t 1-$((${#dataset_names[@]}+1))%${concurrent}" >> $goscript
#
echo "set -e" >> $goscript
if [[ ! -z $crds_context ]]; then
    echo "export CRDS_CONTEXT=\"$crds_context\"" >> $goscript
fi
if [[ ! -z $CONDA_EXE ]] && [[ ! -z $conda_env ]]; then
    echo "source $(dirname $CONDA_EXE)/activate $conda_env" >> $goscript
    echo "export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:/usr/local/lib64:/usr/lib64:/lib64:/usr/lib:/lib:." >> $goscript
fi
#
echo "" >> $goscript
echo "type go-jwst-imaging-stage-1" >> $goscript
echo "type go-jwst-imaging-stage-2" >> $goscript
echo "type go-jwst-imaging-stage-3" >> $goscript
# 
echo "" >> $goscript
echo "dataset_names=( \\" >> $goscript
for (( i=0; i<${#dataset_names[@]}; i++ )); do
    echo "  ${dataset_names[i]} \\" >> $goscript
done
echo ")" >> $goscript
# 
echo "" >> $goscript
echo "if [[ \${PBS_ARRAYID} -le \${#dataset_names[@]} ]]; then" >> $goscript
echo "" >> $goscript
echo "  go-jwst-imaging-stage-1 \${dataset_names[\${PBS_ARRAYID}-1]} --maximum-cores=\"quarter\"" >> $goscript
echo "" >> $goscript
echo "  go-jwst-imaging-stage-2 \${dataset_names[\${PBS_ARRAYID}-1]}" >> $goscript
echo "" >> $goscript
echo "fi" >> $goscript
# 
cat << EOF >> $goscript

if [[ \${PBS_ARRAYID} -gt \${#dataset_names[@]} ]]; then
    waitseconds=\$(awk "BEGIN {print \${#dataset_names[@]}*10;}") # 10 sec per dataset
    echo "sleep \$waitseconds"
    sleep \$waitseconds
    echo "Checking cal files ..."
    echo "go-jwst-imaging-check-cal-files > log_checking_cal_files.txt"
    go-jwst-imaging-check-cal-files > log_checking_cal_files.txt
    while [[ \$(tail -n 1 log_checking_cal_files.txt | grep "All good") -eq 0 ]]; do
        waitseconds2=\$(awk "BEGIN {print \${#dataset_names[@]}*5;}") # 5 sec per dataset
        echo "sleep \$waitseconds2"
        sleep \$waitseconds2
        echo "Checking cal files ..."
        echo "go-jwst-imaging-check-cal-files > log_checking_cal_files.txt"
        go-jwst-imaging-check-cal-files > log_checking_cal_files.txt
    done
    echo "Finally, running go-jwst-imaging-stage-3 ..."
    go-jwst-imaging-stage-3 \${dataset_names[@]}
fi

EOF
# 
echo "Prepared qsub script: $goscript"
while true; do
    read -p "Ready to submit the qsub job? " yn
    case $yn in
        [Yy]* ) echo "Submitting the qsub job!"; qsub $goscript; echo "Job submitted! Exit."; break;;
        [Nn]* ) echo "Not submitting the job! Exit!"; exit;;
        * ) echo "Please answer yes or no.";;
    esac
done





