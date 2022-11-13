#!/bin/bash
# 
# Merge a list of masked rate images.
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
if [[ $# -lt 1 ]]; then
    echo "Please input:"
    echo "    mosaic_image_1_i2d_asn.json"
    echo "    mosaic_image_2_i2d_asn.json"
    echo "Notes:"
    echo "    The inputs are all JWST association files, which should contain many cal images, "
    echo "    for which we will parse the names to get the so-called masked rate images, "
    echo "    then merge them into one data, but excluding each self data."
    exit 255
fi
iarg=1
overwrite=0
mosaic_asn_files=()
while [[ $iarg -le $# ]]; do
    argstr="${!iarg}"
    if [[ "$argstr" == "--overwrite" ]]; then
        overwrite=1
    else
        mosaic_asn_files+=("${!iarg}")
        echo "mosaic_asn_files += \"${!iarg}\""
    fi
    iarg=$((iarg+1))
done


# get all cal images in the asn files
multiobs_cal_images=()
multiobs_rate_images=()
multiobs_masked_rate_images=()
multiobs_dataset_names=()
for (( i = 0; i < ${#mosaic_asn_files[@]}; i++ )); do
    temp_mosaic_asn="${mosaic_asn_files[i]}"
    temp_cal_images=($(cat "$temp_mosaic_asn" | grep "expname" | perl -p -e 's/^.*: ["](.*)["].*/\1/g'))
    for (( m = 0; m < ${#temp_cal_images[@]}; m++ )); do
        cal_image="${temp_cal_images[m]}"
        # fix relative path
        if [[ "$cal_image" == "../"* ]]; then
            #cal_image=$(dirname "$temp_mosaic_asn")/"${temp_cal_images[m]}" # paths in asn are relative to asn dir.
            cal_image=$(dirname $(dirname "$temp_mosaic_asn"))/$(echo "${temp_cal_images[m]}" | perl -p -e 's%^../%%g') # paths in asn are relative to asn dir.
        fi
        rate_image=$(echo "$cal_image" | perl -p -e 's%/calibrated2_cals/%/calibrated1_rates/%g' | perl -p -e 's%_cal.fits$%_rate.fits%g')
        echo "rate_image = \"$rate_image\""
        masked_rate=$(echo "$rate_image" | perl -p -e 's%_rate.fits$%%g')"_masked_source_emission_rate.fits"
        dataset_name=$(basename "$rate_image" | perl -p -e 's%_rate.fits$%%g')
        multiobs_cal_images+=("$cal_image")
        multiobs_rate_images+=("$rate_image")
        multiobs_masked_rate_images+=("$masked_rate")
        multiobs_dataset_names+=("$dataset_name")
    done
done

echo "multiobs_masked_rate_images = ${multiobs_masked_rate_images[@]} (${#multiobs_masked_rate_images[@]})"


# loop each cal/rate image dir, merge all other source-emission-masked rates
for (( i = 0; i < ${#multiobs_rate_images[@]}; i++ )); do
    cal_image="${multiobs_cal_images[i]}"
    rate_image="${multiobs_rate_images[i]}"
    masked_rate_image="${multiobs_masked_rate_images[i]}"
    merged_masked_rate=$(dirname "$masked_rate_image")"/merged_other_visits_masked_source_emission_rate.fits"
    merged_masked_rate_list_file=$(dirname "$masked_rate_image")"/merged_other_visits_masked_source_emission_rate.list.txt"
    merged_masked_rate_updated=0
    output_cal_image=$(echo "$cal_image" | perl -p -e 's/_cal.fits$/_bkgsub_masked_source_emission_cal.fits/g')
    
    if [[ ! -f "$merged_masked_rate" ]] || [[ $overwrite -gt 0 ]]; then
        applicable_masked_rate_images=()
        for (( m = 0; m < ${#multiobs_masked_rate_images[@]}; m++ )); do
            if [[ "${multiobs_masked_rate_images[m]}" != "$masked_rate_image" ]]; then
                applicable_masked_rate_images+=("${multiobs_masked_rate_images[m]}")
            fi
        done
        if [[ ${#applicable_masked_rate_images[@]} -eq 0 ]]; then
            echo "Error! No applicable_masked_rate_images (no other dither/visit)?"
            exit 255
        fi
        
        echo $script_dir/util_merge_source_emission_masked_rate_data.py \
            ${applicable_masked_rate_images[@]} \
            "$merged_masked_rate"
        $script_dir/util_merge_source_emission_masked_rate_data.py \
            ${applicable_masked_rate_images[@]} \
            "$merged_masked_rate"
        echo "# $script_dir/util_merge_source_emission_masked_rate_data.py" > "$merged_masked_rate_list_file"
        for (( m = 0; m < ${#applicable_masked_rate_images[@]}; m++ )); do
            echo "${applicable_masked_rate_images[m]}" >> "$merged_masked_rate_list_file"
        done
        
        merged_masked_rate_updated=1
        
        if [[ -f "$output_cal_image" ]]; then
            mv "$output_cal_image" "$output_cal_image.backup"
        fi
    elif [[ -f "$merged_masked_rate" ]]; then
        echo "Found existing file \"$merged_masked_rate\" and overwrite is False. Skipping."
    fi
    # check result
    if [[ ! -f "$merged_masked_rate" ]]; then
        echo "Error! Failed to produce the output files: \"$merged_masked_rate\""
        exit 255
    else
        echo "Successfully produced the output files: \"$merged_masked_rate\""
    fi
    
    
    # then also re-run stage2
    if [[ ! -f "$output_cal_image" ]] || [[ $overwrite -gt 0 ]] || [[ $merged_masked_rate_updated -gt 0 ]]; then
        proc_args=(--darkobs "$merged_masked_rate")
        if [[ $overwrite -gt 0 ]] || [[ $merged_masked_rate_updated -gt 0 ]]; then
            proc_args+=(--overwrite)
        fi
        echo "*** Running ***" $script_dir/go-jwst-imaging-stage-2-step-3-redo-bkgsub.py \
            "$rate_image" \
            "$output_cal_image" \
            "${proc_args[@]}"
        $script_dir/go-jwst-imaging-stage-2-step-3-redo-bkgsub.py \
            "$rate_image" \
            "$output_cal_image" \
            "${proc_args[@]}"
        if [[ ! -f "$output_cal_image" ]] && [[ ! -L "$output_cal_image" ]]; then
            echo "Error! Failed to produce the output file: $output_cal_image"
            exit 255
        fi
    elif [[ -f "$output_cal_image" ]]; then
        echo "Found existing file \"$output_cal_image\" and overwrite is False and merged_masked_rate_updated is False. Skipping."
    fi
    # check result
    if [[ ! -f "$output_cal_image" ]]; then
        echo "Error! Failed to produce the output files: \"$output_cal_image\""
        exit 255
    else
        echo "Successfully produced the output files: \"$output_cal_image\""
    fi
    
done






