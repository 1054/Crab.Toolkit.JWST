#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

if [[ ! -f $CONDA_PREFIX/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py.ori ]]; then
    cp -i $CONDA_PREFIX/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py \
        $CONDA_PREFIX/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py.ori
fi

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_outlier_detection_numpy_median_two_elements_bug.diff

