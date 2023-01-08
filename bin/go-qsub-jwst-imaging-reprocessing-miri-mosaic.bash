#!/bin/bash
# 
dataset_names=($(ls -1d jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_mirimage))
if [[ ${#dataset_names[@]} -eq 0 ]]; then
    echo "No dataset dir is found: jw*_[0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]_mirimage"
    echo "under current directory: $(pwd -P)"
    exit 255
fi
if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
    echo "Please set conda environment first!"
    exit 255
fi
if [[ -z "$CRDS_CONTEXT" ]]; then
    echo "Please set CRDS_CONTEXT!"
    exit 255
fi
crds_context="$CRDS_CONTEXT" # "jwst_0986.pmap"
conda_env="$CONDA_DEFAULT_ENV" # "jwstpmap1009" # "jwstpmap0995" # "base"
concurrent=20
ncpu=1
mem="50gb" # 
maxcores="none" # "quarter"
timestamp=$(date +"%Y%m%d_%Hh%Mm%Ss")
currentdir=$(pwd -P)
goscript="go_qsub_processing_jwst_imaging_${timestamp}_repreocessing_miri_mosaic.bash"
echo "#!/bin/bash" > $goscript
echo "#PBS -N JW${timestamp}" >> $goscript
echo "#PBS -l nodes=1:ppn=${ncpu},mem=${mem},walltime=48:00:00" >> $goscript
echo "#PBS -d ${currentdir}/" >> $goscript
echo "#PBS -o log_processing_jwst_imaging_${timestamp}_repreocessing_miri_mosaic" >> $goscript
#echo "#PBS -e log_processing_jwst_imaging_${timestamp}_\${PBS_ARRAYID}.err" >> $goscript
echo "#PBS -j oe" >> $goscript # join stdout and stderr
#echo "#PBS -k oe" >> $goscript # keep stdout and stderr on running hostname
#echo "#PBS -m abe" >> $goscript # send email notifications for all, begin, end
echo "#PBS -m n" >> $goscript # do not send emails -- not working
#
echo "set -e" >> $goscript
if [[ ! -z $crds_context ]]; then
    echo "export CRDS_CONTEXT=\"$crds_context\"" >> $goscript
    # also set CRDS_CONTEXT in the conda enviornment activating script
    if [[ ! -z $CONDA_EXE ]] && [[ ! -z $conda_env ]]; then
        if [[ ! -f $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh ]]; then
            if [[ ! -d $CONDA_PREFIX/etc/conda/activate.d ]]; then
                mkdir -p $CONDA_PREFIX/etc/conda/activate.d
            fi
            echo "#!/bin/bash" > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
            chmod +x $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
        fi
        if [[ $(grep CRDS_CONTEXT $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh | wc -l) -eq 0 ]]; then
            echo "export CRDS_CONTEXT=\"$crds_context\"" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
        fi
    fi
fi
if [[ ! -z $CONDA_EXE ]] && [[ ! -z $conda_env ]]; then
    echo "source $(dirname $CONDA_EXE)/activate $conda_env" >> $goscript
    echo "export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:/usr/local/lib64:/usr/lib64:/lib64:/usr/lib:/lib:." >> $goscript
fi
#
echo "" >> $goscript
echo "type go-jwst-imaging-stage-3" >> $goscript
# 
echo "" >> $goscript
echo "dataset_names=( \\" >> $goscript
for (( i=0; i<${#dataset_names[@]}; i++ )); do
    echo "  ${dataset_names[i]} \\" >> $goscript
done
echo ")" >> $goscript
# 
cat << EOF >> $goscript

# Reset the merged dark rate so that the code will combine the dark rates from all possible obs_num and visitnum. 
# In principle we do not need it here now, because \`go-jwst-imaging-stage-3\` with \`--reprocess-miri\` and without
# \`--combine-obsnum\` will already try to find all possible asn files in the sibling directories of the output directory, 
# which means all obs_num and visitnum which has a stage3 asn file created there will be read in. 
# But in case that during the processing of each obs_num and visitnum, some asn files are created later than the time
# the current running process searching for sibling asn files, then the merging of dark rates will be incomplete. 
# So here when we reprocess stage3 with the \`--combine-obsnum\` option, we can re-do the merging of dark rate to make 
# sure all silbing obs_num and visitnum data are considered. 

reset_merged_dark_rate=0

if [[ \$reset_merged_dark_rate -gt 0 ]]; then
    for (( i=0; i<\${#dataset_names[@]}; i++ )); do
        dataset_name="\${dataset_names[i]}"
        if [[ -f "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.fits" ]]; then
            mv "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.fits" \\
               "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.fits.backup"
        fi
        if [[ -f "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.list.txt" ]]; then
            mv "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.list.txt" \\
               "\$dataset_name/calibrated1_rates/merged_other_visits_masked_source_emission_rate.list.txt.backup"
        fi
    done
fi

# Use an absolute reference catalog for the stage3 \`tweakreg\` step for the astrometric alignment. 
# Please prepare your "abs_refcat.fits" FITS-format catalog file in advance if you turn this on!

use_abs_refcat=1

if [[ \$use_abs_refcat -gt 0 ]]; then
    
    echo "*** "
    echo "*** Running: go-jwst-imaging-stage-3 \${dataset_names[@]} --reprocess-miri --run-individual-steps --combine-obsnum --abs-refcat abs_refcat.fits --pixel-scale 0.060"
    echo "*** "
    go-jwst-imaging-stage-3 \${dataset_names[@]} --reprocess-miri --run-individual-steps --combine-obsnum --abs-refcat abs_refcat.fits --pixel-scale 0.060
    
else
    
    echo "*** "
    echo "*** Running: go-jwst-imaging-stage-3 \${dataset_names[@]} --reprocess-miri --run-individual-steps --combine-obsnum --pixel-scale 0.060"
    echo "*** "
    go-jwst-imaging-stage-3 \${dataset_names[@]} --reprocess-miri --run-individual-steps --combine-obsnum --pixel-scale 0.060
    
fi

EOF
# 
echo "Prepared qsub script: $goscript"
while true; do
    read -p "Ready to submit the qsub job? [y/n] " yn
    case $yn in
        [Yy]* ) echo "Submitting the qsub job!"; qsub $goscript; echo "Job submitted! Please check your qstat then!"; break;;
        [Nn]* ) echo "Not submitting the job! Exit!"; exit;;
        * ) echo "Please answer yes or no.";;
    esac
done





