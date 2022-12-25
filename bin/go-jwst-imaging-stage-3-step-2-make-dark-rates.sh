#!/bin/bash
# 
# Make source emission seed image from the input mosaic image, then mask rate image with that seed, 
# so that the masked rate image is like a darkobs.
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

# Read user input
if [[ $# -lt 2 ]]; then
    echo "Please input:"
    echo "    mosaic_image_.fits"
    echo "    mosaic_image_i2d_asn.json"
    echo "Notes:"
    echo "    The second input JWST association file should contain many cal images, "
    echo "    for which we will parse the names to get rate images and mask out "
    echo "    source emission with the mosaic seed, so that they are like darkobs."
    exit 255
fi
iarg=1
overwrite=0
mosaic_image=""
mosaic_asn=""
while [[ $iarg -le $# ]]; do
    argstr="${!iarg}"
    if [[ "$argstr" == "--overwrite" ]]; then
        overwrite=1
    elif [[ "$mosaic_image"x == ""x ]]; then
        mosaic_image="${!iarg}"
        echo "mosaic_image = \"$mosaic_image\""
    elif [[ "$mosaic_asn"x == ""x ]]; then
        mosaic_asn="${!iarg}"
        echo "mosaic_asn = \"$mosaic_asn\""
    fi
    iarg=$((iarg+1))
done



# this_file_in="$mosaic_image"
# this_file_out=$(echo "$this_file_in" | perl -p -e 's/\.fits$//g')"_remstriping_galaxy_seed_image.fits"
# this_dir_out=$(dirname "$this_file_out")
# if [[ ! -f "$this_file_out" ]] || [[ $overwrite -gt 0 ]]; then
#     echo $script_dir/remstriping "$this_file_in" \
#         --no-apply-flat \
#         --no-inplace \
#         --output-dir "$this_dir_out" \
#         --output-suffix "remstriping" \
#         --smooth 1.0
#     $script_dir/remstriping "$this_file_in" \
#         --no-apply-flat \
#         --no-inplace \
#         --output-dir "$this_dir_out" \
#         --output-suffix "remstriping" \
#         --smooth 1.0
# fi
# if [[ ! -f "$this_file_out" ]]; then
#     echo "Error! Failed to produce the output files: \"$this_file_out\""
#     exit 255
# fi

this_file_in="$mosaic_image"
this_file_out=$(echo "$this_file_in" | perl -p -e 's/\.fits$//g')"_galaxy_seed_image.fits"
this_script_out=$(echo "$this_file_in" | perl -p -e 's/\.fits$//g')"_galaxy_seed_image_script.sh"
this_dir_out=$(dirname "$this_file_out")
if [[ ! -f "$this_file_out" ]] || [[ $overwrite -gt 0 ]]; then
    echo "#!/bin/bash" > "$this_script_out"
    echo "#cd $(pwd)" >> "$this_script_out"
    echo "$script_dir/util_detect_source_and_create_seed_image.py \\" >> "$this_script_out"
    echo "    \"$this_file_in\" \\" >> "$this_script_out"
    echo "    --smooth 1.0 \\" >> "$this_script_out"
    if [[ $overwrite -gt 0 ]]; then
        echo "    --overwrite \\" >> "$this_script_out"
    fi
    echo "    " >> "$this_script_out"
    cat "$this_script_out"
    chmod +x "$this_script_out"
    bash "$this_script_out"
fi
if [[ ! -f "$this_file_out" ]]; then
    echo "Error! Failed to produce the output files: \"$this_file_out\""
    exit 255
fi



seed_image="$this_file_out"

# find all associated cal.fits and rate.fits, mask them with source emission mask
mosaic_dir=$(dirname "$mosaic_asn")
cal_images=($(cat "$mosaic_asn" | grep "expname" | perl -p -e 's/^.*: ["](.*)["].*/\1/g'))
rate_images=()
masked_rate_images=()
for (( k = 0; k < ${#cal_images[@]}; k++ )); do
    cal_image="${cal_images[k]}"
    dataset_name=$(basename "$cal_image" | perl -p -e 's/_cal.fits$//g')
    # get cal_image full absolute path
    if [[ "$cal_image" != *"/"* ]]; then
        cal_image="$mosaic_dir/$cal_image"
    fi
    # fix linked path
    if [[ -L "$cal_image" ]]; then
        cal_image=$(realpath "$cal_image")
    fi
    # fix relative path
    if [[ "$cal_image" == "../"* ]]; then
        cal_image=$(dirname "mosaic_dir")/$(echo "${cal_image}" | perl -p -e 's%^../%%g') # paths in asn are relative to asn dir.
    fi
    rate_image=$(echo "$cal_image" | perl -p -e 's%/calibrated2_cals/%/calibrated1_rates/%g' | perl -p -e 's%_cal.fits$%_rate.fits%g')
    masked_rate=$(echo "$rate_image" | perl -p -e 's%_rate\.fits$%%g')"_masked_source_emission_rate.fits"
    masked_rate_script=$(echo "$rate_image" | perl -p -e 's%_rate.fits$%%g')"_masked_source_emission_rate_script.sh"
    # make source-emission-masked rate image
    if [[ ! -f "$masked_rate" ]] || [[ $overwrite -gt 0 ]]; then
        echo "#!/bin/bash" > "$masked_rate_script"
        echo "#cd $(pwd)" >> "$masked_rate_script"
        echo "$script_dir/util_mask_rate_data_with_seed_image.py \\" >> "$masked_rate_script"
        echo "    \"$rate_image\" \\" >> "$masked_rate_script"
        echo "    \"$seed_image\" \\" >> "$masked_rate_script"
        echo "    \"$masked_rate\" \\" >> "$masked_rate_script"
        if [[ $overwrite -gt 0 ]]; then
            echo "    --overwrite \\" >> "$masked_rate_script"
        fi
        echo "    " >> "$masked_rate_script"
        cat "$masked_rate_script"
        chmod +x "$masked_rate_script"
        bash "$masked_rate_script"
    fi
    if [[ ! -f "$masked_rate" ]]; then
        echo "Error! Failed to produce the output files: \"$masked_rate\""
        exit 255
    fi
    rate_images+=("$rate_image")
    masked_rate_images+=("$masked_rate")
done






