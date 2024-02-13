#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

site_packages_dir=$(python -c 'import site; print(site.getsitepackages()[0])')
if [[ x"$site_packages_dir" == x"" ]]; then
    echo "Error! Could not get Python site-packages directory! Please check: "
    echo "    python -c 'import site; print(site.getsitepackages()[0])'"
    exit
fi

if [[ ! -f $site_packages_dir/jwst/tweakreg/tweakreg_step.py ]]; then
    echo "Error! File not found: "
    echo "    $site_packages_dir/jwst/tweakreg/tweakreg_step.py"
    exit
fi

if [[ ! -f $site_packages_dir/jwst/tweakreg/tweakreg_step.py.backup ]]; then
    cp -i $site_packages_dir/jwst/tweakreg/tweakreg_step.py \
        $site_packages_dir/jwst/tweakreg/tweakreg_step.py.backup
fi

patch \
$site_packages_dir/jwst/tweakreg/tweakreg_step.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_tweakreg_step_model_container_from_asn_bug.diff

