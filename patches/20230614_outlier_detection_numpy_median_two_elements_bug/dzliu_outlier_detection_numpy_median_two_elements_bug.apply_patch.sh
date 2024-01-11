#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

script_dir=$(dirname ${BASH_SOURCE[0]})
site_packages_dir=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')

if [[ ! -f ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py.ori ]]; then
    cp -i ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py \
          ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py.ori
fi

patch \
    ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py \
    < ${script_dir}/dzliu_outlier_detection_numpy_median_two_elements_bug.diff

