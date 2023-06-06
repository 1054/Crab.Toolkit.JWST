#!/bin/bash
# 
set -e


# Usage
usage () {
    echo "Usage: "
    echo "    go-run-random-aperture-statistics.sh \\"
    echo "        -group f115w \\"
    echo "        -images /path/to/*f115w*.fits \\"
    echo "        [-apercorr 1.0] \\"
    echo "        [-outdir output_random_aperture_statistics] \\"
    echo "        [-aperdiam 0.30] \\"
    echo "Notes: "
    echo "    Arguments in brackets are optional. Their default values are shown as above."
}


# Add sys path
script_dir=$(dirname $(realpath "${BASH_SOURCE[0]}"))
if [[ "$PATH" != *":$script_dir"* ]]; then
    export PATH="$PATH:$script_dir"
fi


# Read User input
iarg=1
argstr=""
argstr2=""
argmode="none"
argtype="none"
igroup=-1
group_names=()
aper_corrs=()
n_images=()
input_images=()
total_n_images=0
aperture_diameter=0.30 # arcsec
output_dir="output_random_aperture_statistics"
while [[ $iarg -le $# ]]; do
    argstr="${!iarg}"
    argstr2=$(echo "$argstr" | tr '[:upper:]' '[:lower:]' | perl -p -e 's/^[-]+/-/g')
    if [[ "$argstr2" == "-"* ]]; then
        # arg is an option
        if [[ "$argstr2" == "-group" ]] || [[ "$argstr2" == "-filter" ]]; then
            argmode="group"
        elif [[ "$argstr2" == "-apercorr" ]]; then
            argmode="apercorr"
        elif [[ "$argstr2" == "-aperdiam" ]]; then
            argmode="aperdiam"
        elif [[ "$argstr2" == "-outdir" ]]; then
            argmode="outdir"
        elif [[ "$argstr2" == "-images" ]]; then
            argmode="images"
        else
            argmode="none"
            echo "Unrecognized option: $argstr"
        fi
    else
        # arg is an input
        if [[ "$argmode" == "group" ]]; then
            group_names+=("$argstr")
            aper_corrs+=(1.0)
            n_images+=(0)
            igroup=$((igroup+1))
        elif [[ "$argmode" == "apercorr" ]]; then
            aper_corrs[igroup]=$argstr
        elif [[ "$argmode" == "aperdiam" ]]; then
            aperture_diameter=$argstr
        elif [[ "$argmode" == "outdir" ]]; then
            output_dir="$argstr"
        elif [[ "$argmode" == "images" ]]; then
            input_images+=("$argstr")
            #echo "n_images (${#n_images[@]}): ${n_images[@]}"
            n_images[igroup]=$((n_images[igroup]+1))
            total_n_images=$((total_n_images+1))
        else
            echo "Unrecognized input: $argstr"
        fi
    fi
    iarg=$((iarg+1))
done

igroup=$((igroup+1))

if [[ $igroup -eq 0 ]] || [[ $total_n_images -eq 0 ]]; then
    usage
    exit
fi

echo "Read ${igroup} group(s) of images."


# Make output directory
if [[ ! -d "$output_dir" ]]; then
    echo mkdir -p "$output_dir"
    mkdir -p "$output_dir"
fi


# Loop input groups
# Each group has N images
i_image=0
for (( igroup=0; igroup<${#group_names[@]}; igroup++ )); do
    group_name=${group_names[igroup]}
    group_name2=$(echo "$group_name" | perl -p -e 's/[^a-zA-Z0-9_-]/_/g' | perl -p -e 's/[_]+$//g')
    #group_uppercase=$(echo "${group_name}" | tr '[:lower:]' '[:upper:]')
    n_image=${n_images[igroup]}
    aper_corr=${aper_corrs[igroup]}
    image_files=()
    i2_image=$((i_image+n_image))
    while [[ $i_image -lt $i2_image ]]; do
        image_files+=(${input_images[i_image]})
        i_image=$((i_image+1))
    done
    
    out_dir="${output_dir}/${group_name2}"
    out_pdf="${output_dir}_${group_name2}.pdf"
    
    do_individual_image=0
    
    if [[ $do_individual_image -gt 0 ]]; then
        for (( k=0; k<${#image_files[@]}; k++ )); do
            image_file="${image_files[k]}"
            image_name=$(basename "$image_file" | perl -p -e 's/\.(fits|fits.gz)$//g' | perl -p -e 's/[^a-zA-Z0-9_-]/_/g' | perl -p -e 's/[_]+$//g')
            if [[ ! -f "${image_file}" ]]; then
                echo "Warning! Image file is not found: \"${image_file}\""
                continue
            fi
            k_image=$((k+1))
            out_path="${output_dir}/${group_name2}/${k_image}"
            echo "Running: util_compare_image_reduction.py \\"
            echo "    \"${image_file}\" \\"
            echo "    --aperture-size $aperture_diameter \\"
            echo "    --aperture-number 10000 \\"
            echo "    --aper-corr $aper_corr \\"
            echo "    --input-labels \"$k_image: $image_name\" \\"
            echo "    --no-legend \\"
            echo "    \"$out_path\""
            util_compare_image_reduction.py \
                "${image_file}" \
                --aperture-size $aperture_diameter \
                --aperture-number 10000 \
                --aper-corr $aper_corr \
                --input-labels "$k_image: $image_name" \
                --no-legend \
                "$out_path"
        done

        echo "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${out_pdf}  ${out_dir}/*.pdf"
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${out_pdf}  ${out_dir}/*.pdf
    
    else
        
        labels=""
        title=""
        for (( k=0; k<${#image_files[@]}; k++ )); do
            image_file="${image_files[k]}"
            image_name=$(basename "$image_file" | perl -p -e 's/\.(fits|fits.gz)$//g' | perl -p -e 's/[^a-zA-Z0-9_-]/_/g' | perl -p -e 's/[_]+$//g')
            k_image=$((k+1))
            if [[ $k -eq 0 ]]; then
                labels="$k_image"
                title="$k_image: $image_name"
            else
                labels="$labels, $k_image"
                title="$title, $k_image: $image_name"
            fi
        done
        out_path="${output_dir}/${group_name2}/all"
        echo "Running: util_compare_image_reduction.py \\"
        echo "    ${image_files[@]} \\"
        echo "    --aperture-size $aperture_diameter \\"
        echo "    --aperture-number 10000 \\"
        echo "    --aper-corr $aper_corr \\"
        echo "    --input-labels \"$labels\" \\"
        echo "    --title \"$title\" \\"
        echo "    --no-legend \\"
        echo "    \"$out_path\""
        util_compare_image_reduction.py \
            "${image_files[@]}" \
            --aperture-size $aperture_diameter \
            --aperture-number 10000 \
            --aper-corr $aper_corr \
            --input-labels "$labels" \
            --title "$title" \
            --no-legend \
            "$out_path"
        
    fi

done


echo "Running: $script_dir/go-run-random-aperture-statistics-collect-result-table.py \\"
echo "    $output_dir/*/*_stats.json \\"
echo "    $output_dir/output_rand_aper_stats.csv"
go-run-random-aperture-statistics-collect-result-table.py \
    $output_dir/*/*_stats.json \
    $output_dir/output_rand_aper_stats.csv



