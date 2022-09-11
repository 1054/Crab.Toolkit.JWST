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
if [[ $# -lt 2 ]]; then
    echo "Please input:"
    echo "    jwst dataset name, e.g. jw01837022001_02201_00002_nrca1"
    echo "    rate images, e.g. jw01837022001_02201_00002_nrca1_rate.fits jw01837022001_04201_00002_nrca1_rate.fits"
    echo "Output:"
    echo "    merged/averaged source-emission-masked rate image, "
    echo "    write as {dataset_name}/calibrated1_rates/merged_other_visits_mirimage_rate_masked_source_emission.fits"
    echo "Note:"
    echo "    we will exclude the input dataset_name if it is also in the input rate images."
    exit 255
fi
input_dataset_name="$1"

# validate input_dataset_name
$script_dir/go-jwst-parse-dataset-name.py "$input_dataset_name"
if [[ $? -ne 0 ]]; then
    echo "Error! The input dataset name \"$input_dataset_name\" seems incorrect!"
    exit 255
fi

# Read input rate images
input_rate_images=()
for (( i = 2; i <= $#; i++ )); do
    this_rate_image="${!i}"
    if [[ "$this_rate_image" == *"$input_dataset_name"* ]]; then
        continue
    fi
    input_rate_images+=("${!i}")
done

# 
if [[ ${#input_rate_images[@]} -eq 0 ]]; then
    echo "Error! No valid input rate images!"
    exit 255
fi

# check output dir
if [[ ! -d "$input_dataset_name/calibrated1_rates" ]]; then
    echo "Error! Dataset directory does not exist: \"$input_dataset_name/calibrated1_rates\""
    exit 255
fi

# 
output_rate_image="$input_dataset_name/calibrated1_rates/merged_other_visits_mirimage_rate_masked_source_emission.fits"

# 
echo $script_dir/util_merge_source_emission_masked_rate_data.py \
    ${input_rate_images[@]} \
    $output_rate_image
$script_dir/util_merge_source_emission_masked_rate_data.py \
    ${input_rate_images[@]} \
    $output_rate_image


