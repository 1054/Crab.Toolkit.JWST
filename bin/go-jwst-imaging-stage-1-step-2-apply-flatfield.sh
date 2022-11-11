#!/bin/bash
# 
# Apply flat field. 
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
    echo "Notes:"
    echo "    We will query CRDS to get flatfield and apply it to the input rate file."
    echo "    A history will be written into the input rate file as well."
    exit 255
fi
# 
echo $script_dir/util_apply_flatfield.py $@
$script_dir/util_apply_flatfield.py $@
# 
if [[ $? -ne 0 ]]; then
    echo "Error?! Please check previous messages."
    exit 255
fi


