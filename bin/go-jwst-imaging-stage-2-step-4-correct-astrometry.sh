#!/bin/bash
# 
# Correct astrometry for specific data sets with custom wcs correction json 
# files under "util_correct_astrometry_for_specific_data_sets_wcs_files".
# 

# Set necessary system variables
if [[ -z "$CRDS_PATH" ]]; then
    export CRDS_PATH=$HOME/jwst_crds_cache
fi
if [[ -z "$CRDS_SERVER_URL" ]]; then
    export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
fi

# Get script_dir
script_dir=$(dirname "${BASH_SOURCE[0]}")

# Run 
if [[ $# -eq 0 ]]; then
    echo "Please input a JWST imaging cal file, e.g.,"
    echo "    ./this_script jw01837004020_02201_00001_nrca1_cal.fits"
    exit 255
fi
# 
echo $script_dir/util_correct_astrometry_for_specific_data_sets.py $@
$script_dir/util_correct_astrometry_for_specific_data_sets.py $@
# 
if [[ $? -ne 0 ]]; then
    echo "Error?! Please check previous messages."
    exit 255
fi


