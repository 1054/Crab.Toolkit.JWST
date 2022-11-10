#!/bin/bash
# 
# Removing horizontal/vertical/angled strips in a NIRCam rate image.
# 
# 20221002: can remove angle=+-30deg strips using remstriping-dzliu.
# 20221109: renamed the called script. added a prestep global sigma clipping.
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
    echo "Options:"
    echo "    --apply-flat : for MIRI"
    exit 255
fi
if [[ "$1" == *"nrca"* ]] || [[ "$1" == *"nrcb"* ]] || [[ "$1" == *"nrclong"* ]]; then
    echo $script_dir/util_remove_stripes_along_four_angles.py $@
    $script_dir/util_remove_stripes_along_four_angles.py $@
    if [[ $? -ne 0 ]]; then
        echo "Error?! Please check previous messages."
        exit 255
    fi
else
    echo "The input data is not a NIRCam imaging data. No need to remove the stripes (1/f noise)."
fi


