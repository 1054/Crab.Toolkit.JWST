#!/bin/bash
# 
set -e

script_dir=$(dirname "${BASH_SOURCE[0]}")

if [[ $# -eq 0 ]]; then
    echo "Usage: "
    echo "    Please input fits files."
    exit 255
fi

echo "script_dir: $script_dir"

for (( i = 1; i <= $#; i++ )); do
    filepath="${!i}"
    if [[ "$filepath" == *".fits" ]]; then
        echo $script_dir/util_make_quicklook_image.py "$filepath" "($i/$#)"
        $script_dir/util_make_quicklook_image.py "$filepath"
    fi
done

