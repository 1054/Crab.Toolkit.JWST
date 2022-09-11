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

# Read input
if [[ $# -ne 2 ]]; then
    echo "Please input:"
    echo "    i2d image, e.g. jw01837_obs022_NIRCAM_F090W_i2d.fits"
    echo "    rate image, e.g. jw01837022001_02201_00002_nrca1_rate.fits"
    echo "Output:"
    echo "    source emission masked rate image, e.g. jw01837022001_02201_00002_nrca1_rate_.fits"
    exit 255
fi
input_image_file="$1"
seed_image_file=$(echo "$input_image_file" | perl -p -e 's/\.fits$//g')"_remstriping_galaxy_seed_image.fits"

input_rate_file="$2"
output_rate_file=$(echo "$input_rate_file" | perl -p -e 's/\.fits$//g')"_masked_source_emission.fits"

# Run remstriping
if ([[ "$input_image_file" == *"nrc"* ]] || [[ "$input_image_file" == *"NIRCAM"* ]]) && [[ "$input_image_file" != *"miri"* ]]; then
    echo $script_dir/remstriping "$input_image_file" --apply-flat --no-inplace --output-suffix remstriping --smooth 3.0
    $script_dir/remstriping "$input_image_file" --apply-flat --no-inplace --output-suffix remstriping --smooth 3.0
else
    echo $script_dir/remstriping "$input_image_file" --no-apply-flat --no-inplace --output-suffix remstriping --smooth 3.0
    $script_dir/remstriping "$input_image_file" --no-apply-flat --no-inplace --output-suffix remstriping --smooth 3.0
fi

if [[ ! -f "$seed_image_file" ]]; then
    echo "Error! Failed to produce the output file: \"$seed_image_file\"!"
    exit 255
fi

# Mask source emission for the input rate data
echo $script_dir/util_mask_rate_data_with_seed_image.py \
    "$input_rate_file" \
    "$seed_image_file" \
    "$output_rate_file"
$script_dir/util_mask_rate_data_with_seed_image.py \
    "$input_rate_file" \
    "$seed_image_file" \
    "$output_rate_file"

if [[ ! -f "$output_rate_file" ]]; then
    echo "Error! Failed to produce the output file: \"$output_rate_file\"!"
    exit 255
fi


