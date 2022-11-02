#!/bin/bash
# 
set -e

script_dir=$(dirname "${BASH_SOURCE[0]}")

for (( i = 1; i <= $#; i++ )); do
    filepath="${!i}"
    if [[ "$filepath" == *".fits" ]]; then
        echo $script_dir/go-make-quicklook-image.py "$filepath" "($i/$#)"
        $script_dir/go-make-quicklook-image.py "$filepath"
    fi
done

