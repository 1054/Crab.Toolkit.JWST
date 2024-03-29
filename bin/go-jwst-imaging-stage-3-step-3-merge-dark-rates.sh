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
date_diff=7 # days, 
subtract_as_wisps=0 # subtract as wisps, i.e., a scaleable background, rather than as master dark rate.
mosaic_asn_files=()
while [[ $iarg -le $# ]]; do
    argstr="${!iarg}"
    if [[ "$argstr" == "--overwrite" ]]; then
        overwrite=1
    elif [[ "$argstr" == "--date-diff" ]]; then
        iarg=$((iarg+1))
        if [[ $iarg -le $# ]]; then
            date_diff="${!iarg}"
            echo "date_diff = \"$date_diff\""
        fi
    elif [[ "$argstr" == "--subtract-as-wisps" ]]; then
        subtract_as_wisps=0
    else
        mosaic_asn_files+=("${!iarg}")
        echo "mosaic_asn_files += \"${!iarg}\""
    fi
    iarg=$((iarg+1))
done


# define function to get detector name
function get_detector_name() {
    if [[ $# -gt 0 ]]; then
        file_basename=$(basename "$1")
        detector_name=$(echo "$file_basename" | perl -p -e 's/.*_((nrc(a|b)(1|2|3|4|long))|mirimage)_.*/\1/g')
        if [[ "$detector_name" != "$1" ]]; then
            echo "$detector_name"
        else
            echo "Unknown"
        fi
    else
        echo "Unknown"
    fi
}


# get all cal images in the asn files
multiobs_cal_images=()
multiobs_rate_images=()
multiobs_masked_rate_images=()
multiobs_dataset_names=()
multiobs_detector_names=()
for (( i = 0; i < ${#mosaic_asn_files[@]}; i++ )); do
    temp_mosaic_asn="${mosaic_asn_files[i]}"
    temp_mosaic_dir=$(dirname "$temp_mosaic_asn")
    temp_cal_images=($(cat "$temp_mosaic_asn" | grep "expname" | perl -p -e 's/^.*: ["](.*)["].*/\1/g'))
    for (( m = 0; m < ${#temp_cal_images[@]}; m++ )); do
        cal_image="${temp_cal_images[m]}"
        dataset_name=$(basename "$cal_image" | perl -p -e 's/_cal.fits$//g')
        # get cal_image full absolute path if it is only a filename, i.e., does not contain a path separator
        if [[ "$cal_image" != *"/"* ]]; then
            cal_image="$temp_mosaic_dir/$cal_image"
        fi
        # fix linked path
        if [[ -L "$cal_image" ]]; then
            cal_image=$(realpath "$cal_image")
        fi
        # fix relative path
        if [[ "$cal_image" == "../"* ]]; then
            cal_image=$(dirname "$temp_mosaic_dir")/$(echo "${cal_image}" | perl -p -e 's%^../%%g') # paths in asn are relative to asn dir.
        fi
        # check duplicates
        if [[ x" ${multiobs_cal_images[@]} "x == x*" $cal_image "*x ]]; then
            continue
        fi
        rate_image=$(echo "$cal_image" | perl -p -e 's%/calibrated2_cals/%/calibrated1_rates/%g' | perl -p -e 's%_cal.fits$%_rate.fits%g')
        echo "rate_image = \"$rate_image\""
        masked_rate=$(echo "$rate_image" | perl -p -e 's%_rate.fits$%%g')"_masked_source_emission_rate.fits"
        dataset_name=$(basename "$rate_image" | perl -p -e 's%_rate.fits$%%g')
        detector_name=$(get_detector_name "$rate_image")
        multiobs_cal_images+=("$cal_image")
        multiobs_rate_images+=("$rate_image")
        multiobs_masked_rate_images+=("$masked_rate")
        multiobs_dataset_names+=("$dataset_name")
        multiobs_detector_names+=("$detector_name")
    done
done

echo "multiobs_masked_rate_images = ${multiobs_masked_rate_images[@]} (${#multiobs_masked_rate_images[@]})"


# loop each cal/rate image dir, merge all other source-emission-masked rates
for (( i = 0; i < ${#multiobs_rate_images[@]}; i++ )); do
    cal_image="${multiobs_cal_images[i]}"
    rate_image="${multiobs_rate_images[i]}"
    detector_name="${multiobs_detector_names[i]}"
    masked_rate_image="${multiobs_masked_rate_images[i]}"
    merged_masked_rate=$(dirname "$masked_rate_image")"/merged_other_visits_masked_source_emission_rate.fits"
    merged_masked_rate_list_file=$(dirname "$masked_rate_image")"/merged_other_visits_masked_source_emission_rate.list.txt"
    merged_masked_rate_script_file=$(dirname "$masked_rate_image")"/merged_other_visits_masked_source_emission_rate.script.sh"
    merged_masked_rate_updated=0
    output_cal_image=$(echo "$cal_image" | perl -p -e 's/_cal.fits$/_bkgsub_masked_source_emission_cal.fits/g')
    
    if [[ ! -f "$merged_masked_rate" ]] || [[ $overwrite -gt 0 ]]; then
        applicable_masked_rate_images=()
        for (( m = 0; m < ${#multiobs_masked_rate_images[@]}; m++ )); do
            if [[ "${multiobs_masked_rate_images[m]}" != "$masked_rate_image" ]] && \
               [[ "${multiobs_detector_names[m]}" == "$detector_name" ]]; then
                applicable_masked_rate_images+=("${multiobs_masked_rate_images[m]}")
            fi
        done
        if [[ ${#applicable_masked_rate_images[@]} -eq 0 ]]; then
            echo "Error! No applicable_masked_rate_images (no other dither/visit)?"
            exit 255
        fi
        
        echo $script_dir/util_merge_source_emission_masked_rate_data.py \
            --check-date "$rate_image" \
            --date-diff $date_diff \
            ${applicable_masked_rate_images[@]} \
            "$merged_masked_rate"
        $script_dir/util_merge_source_emission_masked_rate_data.py \
            --check-date "$rate_image" \
            --date-diff $date_diff \
            ${applicable_masked_rate_images[@]} \
            "$merged_masked_rate"
        echo "# $script_dir/util_merge_source_emission_masked_rate_data.py" > "$merged_masked_rate_list_file"
        echo "#!/bin/bash" > "$merged_masked_rate_script_file"
        echo "cd \$(dirname \${BASH_SOURCE[0]})" >> "$merged_masked_rate_script_file"
        echo "$script_dir/util_merge_source_emission_masked_rate_data.py \\" >> "$merged_masked_rate_script_file"
        echo "  --check-date \"$rate_image\" \\" >> "$merged_masked_rate_script_file"
        echo "  --date-diff $date_diff \\" >> "$merged_masked_rate_script_file"
        for (( m = 0; m < ${#applicable_masked_rate_images[@]}; m++ )); do
            echo "${applicable_masked_rate_images[m]}" >> "$merged_masked_rate_list_file"
            echo "  \"${applicable_masked_rate_images[m]}\" \\" >> "$merged_masked_rate_script_file"
        done
        echo "  \"$merged_masked_rate\"" >> "$merged_masked_rate_script_file"
        chmod +x "$merged_masked_rate_script_file"
        
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
    
    
    # then re-run stage2 with an associated background so that `calwebb_image2.Image2Pipeline.bkg_subtract` will be used.
    if [[ ! -f "$output_cal_image" ]] || [[ $overwrite -gt 0 ]] || [[ $merged_masked_rate_updated -gt 0 ]]; then
        
        if [[ $subtract_as_wisps -eq 0 ]]; then
            
            # subtract the merged dark as master dark using jwst pipeline
            
            proc_args=(--darkobs "$merged_masked_rate" --skymatch)
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
        
        else
            
            # subtract the merged dark as wisps, i.e., a scaleable background, using 'util_remove_wisps_with_templates.py'
            # this is updating "$rate_image" in-place
            proc_args=(--template-file "$merged_masked_rate")
            echo "*** Running ***" $script_dir/util_remove_wisps_with_templates.py \
                "$rate_image" \
                "${proc_args[@]}"
            $script_dir/util_remove_wisps_with_templates.py \
                "$rate_image" \
                "${proc_args[@]}"
            
            # redo the rate->cal stage2
            proc_args=()
            if [[ $overwrite -gt 0 ]] || [[ $merged_masked_rate_updated -gt 0 ]]; then
                proc_args+=(--overwrite)
            fi
            echo "*** Running ***" $script_dir/go-jwst-imaging-stage-2-step-1.py \
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






