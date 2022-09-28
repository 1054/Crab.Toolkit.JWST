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
echo "#PBS -m abe" >> $goscript
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
echo "type go-jwst-imaging-stage-1" >> $goscript
echo "type go-jwst-imaging-stage-2" >> $goscript
echo "type go-jwst-imaging-stage-3" >> $goscript
# 
for (( i=0; i<${#dataset_names[@]}; i++ )); do
    dataset_name="${dataset_names[i]}"
    echo "" >> $goscript
    echo "go-jwst-imaging-stage-1 $dataset_name" >> $goscript
    echo "" >> $goscript
    echo "go-jwst-imaging-stage-2 $dataset_name" >> $goscript
done
echo "" >> $goscript
echo "go-jwst-imaging-stage-3 \\" >> $goscript
for (( i=0; i<${#dataset_names[@]}; i++ )); do
    if [[ $((i+1)) -lt ${#dataset_names[@]} ]]; then
        echo "    ${dataset_names[i]} \\" >> $goscript
    else
        echo "    ${dataset_names[i]}" >> $goscript
    fi
done
echo "" >> $goscript





