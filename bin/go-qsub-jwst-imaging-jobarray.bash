#!/bin/bash
# 
dataset_names=($(ls -1d jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_*))
if [[ ${#dataset_names[@]} -eq 0 ]]; then
    echo "No dataset dir is found: jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_*"
    echo "under current directory: $(pwd -P)"
    exit 255
fi
crds_context="" # "jwst_0986.pmap"
conda_env="base"
ncpu=4
mem="40gb"
timestamp=$(date +"%Y%m%d_%Hh%Mm%Ss")
currentdir=$(pwd -P)
goscript="go_qsub_processing_jwst_imaging_${timestamp}.bash"
echo "#!/bin/bash" > $goscript
echo "#PBS -N JW${timestamp}" >> $goscript
echo "#PBS -l nodes=1:ppn=${ncpu},mem=${mem},walltime=24:00:00" >> $goscript
echo "#PBS -d ${currentdir}/" >> $goscript
echo "#PBS -o ${currentdir}/log_processing_jwst_imaging_${timestamp}.txt" >> $goscript
echo "#PBS -e ${currentdir}/log_processing_jwst_imaging_${timestamp}.err" >> $goscript
echo "#PBS -j oe" >> $goscript
echo "#PBS -k oe" >> $goscript
#echo "#PBS -m abe" >> $goscript
echo "#PBS -m n" >> $goscript # do not send emails -- not working
echo "#PBS -t 1-$((${#dataset_names[@]}+1))" >> $goscript
#
echo "set -e" >> $goscript
if [[ ! -z $crds_context ]]; then
    echo "export CRDS_CONTEXT=\"$crds_context\"" >> $goscript
fi
if [[ ! -z $CONDA_PREFIX ]] && [[ ! -z $conda_env ]]; then
    echo "source $CONDA_PREFIX/bin/activate base" >> $goscript
    echo "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:/usr/local/lib64:/usr/lib64:/lib64:/usr/lib:/lib:." >> $goscript
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
echo "  go-jwst-imaging-stage-1 \${dataset_names[\${PBS_ARRAYID}-1]}" >> $goscript
echo "" >> $goscript
echo "  go-jwst-imaging-stage-2 \${dataset_names[\${PBS_ARRAYID}-1]}" >> $goscript
echo "" >> $goscript
echo "fi" >> $goscript
# 
echo "" >> $goscript
echo "if [[ \${PBS_ARRAYID} -gt \${#dataset_names[@]} ]]; then" >> $goscript
echo "" >> $goscript
echo "  echo Then please run:" >> $goscript
echo "  echo go-jwst-imaging-stage-3 \${dataset_names[@]}" >> $goscript
echo "" >> $goscript
echo "fi" >> $goscript
echo "" >> $goscript
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





