#!/bin/bash
# 
set -e

script_dir=$(dirname "${BASH_SOURCE[0]}")

if [[ $# -eq 0 ]]; then
    echo "Usage: "
    echo "    Please input a \"calibrated3_mosaics/<dataset_name>\" directory."
    exit 255
fi

echo "script_dir: $script_dir"

input_dir=$(echo "$1" | perl -p -e 's%/$%%g')
echo "input_dir: $input_dir"

dataset_name=$(basename "$input_dir")
echo "dataset_name: $dataset_name"

if [[ "$input_dir" == *"_i2d.fits" ]]; then
    input_dir=$(dirname "$input_dir")
    echo "input_dir: $input_dir"
    dataset_name=$(echo "$dataset_name" | perl -p -e 's/_i2d.fits$//g')
    echo "dataset_name: $dataset_name"
fi

# find "catfile.txt"
if [[ ! -f "${input_dir}/catfile.txt" ]]; then
    echo "Error! File not found \"${input_dir}/catfile.txt\""
    exit 255
fi

# calfiles
cal_files=($(cat "${input_dir}/catfile.txt" | grep -v '^#' | awk '{print $1}'))
csv_files=($(cat "${input_dir}/catfile.txt" | grep -v '^#' | awk '{print $2}'))
reg_files=()

for (( i = 0; i < ${#csv_files[@]}; i++ )); do
    csv_file="${csv_files[i]}"
    reg_file=$(echo "${csv_file}" | perl -p -e 's/\.csv$/.ds9.reg/g')
    if [[ ! -f "$input_dir/$reg_file" ]]; then
        echo util_convert_catalog_x_y_to_ds9_region.py \
            "$input_dir/$csv_file" \
            "$input_dir/$reg_file"
        util_convert_catalog_x_y_to_ds9_region.py \
            "$input_dir/$csv_file" \
            "$input_dir/$reg_file"
    fi
    if [[ ! -f "$input_dir/$reg_file" ]]; then
        echo "Error! Failed to produce \"$input_dir/$reg_file\"!"
        exit 255
    fi
    reg_files+=("$reg_file")
done

# script
script_file="${input_dir}/examine_image_alignments_with_ds9.sh"
if [[ -f "$script_file" ]]; then
    mv "$script_file" "$script_file.backup"
fi
echo "#!/bin/bash" > "$script_file"
echo "set -e" >> "$script_file"
echo "cd \$(dirname \${BASH_SOURCE[0]})" >> "$script_file"
echo "ds9 -scale zscale \\" >> "$script_file"
echo "    \"${dataset_name}_i2d.fits\" \\" >> "$script_file"
echo "        -zoom to fit \\" >> "$script_file"
for (( i = 0; i < ${#cal_files[@]}; i++ )); do
    cal_filename=$(basename "${cal_files[i]}")
    reg_filename=$(basename "${reg_files[i]}")
    echo "    \"${cal_filename}\" \\" >> "$script_file"
    echo "        -regions load \"${reg_filename}\" \\" >> "$script_file"
    echo "        -zoom to fit \\" >> "$script_file"
done
echo "    " >> "$script_file"
echo "" >> "$script_file"
chmod +x "$script_file"

echo "Prepared script: \"$script_file\""
echo "Executing script: \"$script_file\""
$script_file


