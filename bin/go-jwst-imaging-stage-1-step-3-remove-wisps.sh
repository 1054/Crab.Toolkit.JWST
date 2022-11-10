#!/bin/bash
# 
# Fitting and removing wisps using the NIRCam team's wisps templates, 
# see https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-features-and-caveats/nircam-claws-and-wisps
# template data are downloaded from https://stsci.box.com/s/1bymvf1lkrqbdn9rnkluzqk30e8o2bne
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
    echo "Please input a JWST imaging rate file, e.g.,"
    echo "    ./this_script jw01837022001_02201_00002_nrca1_rate.fits"
    echo "Output: "
    echo "    The input rate file will be updated inplace. "
    echo "    A backup will be created with a suffix of \"*_before_removing_wisps.fits\""
    exit 255
fi
if [[ "$1" == *"nrca"* ]] || [[ "$1" != *"nrcb"* ]] || [[ "$1" != *"nrclong"* ]]; then
    echo $script_dir/util_remove_wisps_with_templates.py $@
    $script_dir/util_remove_wisps_with_templates.py $@
    # 
    if [[ $? -ne 0 ]]; then
        echo "Error?! Please check previous messages."
        exit 255
    fi
else
    echo "The input data is not a NIRCam imaging data. No need to remove wisps."
fi


