#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

script_dir=$(dirname ${BASH_SOURCE[0]})
site_packages_dir=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')

if [[ ! -f ${site_packages_dir}/mirage/ramp_generator/obs_generator.py.ori ]]; then
    cp -i ${site_packages_dir}/mirage/ramp_generator/obs_generator.py \
          ${site_packages_dir}/mirage/ramp_generator/obs_generator.py.ori
fi

patch \
    ${site_packages_dir}/mirage/ramp_generator/obs_generator.py \
    < ${script_dir}/dzliu_mirage_obs_generator_signalgain_nan_lam_too_large_fix_patch.txt


# 20240619: still needed for jwst==1.14.0 and mirage==2.4.0
