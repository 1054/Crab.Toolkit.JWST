#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

if [[ ! -f $CONDA_PREFIX/lib/python3.9/site-packages/jwst/tweakreg/tweakreg_step.py.backup ]]; then
    cp -i $CONDA_PREFIX/lib/python3.9/site-packages/jwst/tweakreg/tweakreg_step.py \
        $CONDA_PREFIX/lib/python3.9/site-packages/jwst/tweakreg/tweakreg_step.py.backup
fi

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/jwst/tweakreg/tweakreg_step.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_tweakreg_step_model_container_from_asn_bug.diff

