#!/bin/bash
# 
# Removing horizontal/vertical strips in a NIRCam rate image.
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

# Run remstriping
if [[ $# -eq 0 ]]; then
    echo "Please input a JWST dataset name, e.g.,"
    echo "    ./this_script jw01837022001_02201_00002_nrca1_rate.fits"
    exit 255
fi
if [[ "$1" == *"nrc"* ]] && [[ "$1" != *"miri"* ]]; then
    echo $script_dir/remstriping-dzliu $@ --apply-flat --inplace
    $script_dir/remstriping-dzliu $@ --apply-flat --inplace
else
    echo "The input data is not a NIRCam imaging data. Will not run remstriping."
fi


