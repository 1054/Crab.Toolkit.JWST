#!/bin/bash
# 
set -e

if [[ -z $CONDA_PREFIX ]]; then
    echo "Please set your conda environment \$CONDA_PREFIX first!"
    exit 255
fi

script_dir=$(dirname ${BASH_SOURCE[0]})
site_packages_dir=$(python -c 'from distutils.sysconfig import get_python_lib; print(get_python_lib())')

if [[ ! -f ${site_packages_dir}/mirage/seed_image/fits_seed_image.py.ori ]]; then
    cp -i ${site_packages_dir}/mirage/seed_image/fits_seed_image.py \
          ${site_packages_dir}/mirage/seed_image/fits_seed_image.py.ori
fi

patch \
    ${site_packages_dir}/mirage/seed_image/fits_seed_image.py \
    < ${script_dir}/dzliu_mirage_input_mosaic_image_fix_patch1.txt


if [[ ! -f ${site_packages_dir}/mirage/seed_image/blot_image.py.ori ]]; then
    cp -i ${site_packages_dir}/mirage/seed_image/blot_image.py \
          ${site_packages_dir}/mirage/seed_image/blot_image.py.ori
fi

patch \
    ${site_packages_dir}/mirage/seed_image/blot_image.py \
    < ${script_dir}/dzliu_mirage_input_mosaic_image_fix_patch2.txt


if [[ ! -f ${site_packages_dir}/mirage/seed_image/save_seed.py.ori ]]; then
    cp -i ${site_packages_dir}/mirage/seed_image/save_seed.py \
          ${site_packages_dir}/mirage/seed_image/save_seed.py.ori
fi

patch \
    ${site_packages_dir}/mirage/seed_image/save_seed.py \
    < ${script_dir}/dzliu_mirage_input_mosaic_image_fix_patch3.txt


#<20240111><NoNeed># if [[ ! -f ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py.ori ]]; then
#<20240111><NoNeed>#     cp -i ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py \
#<20240111><NoNeed>#           ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py.ori
#<20240111><NoNeed># fi
#<20240111><NoNeed># 
#<20240111><NoNeed># patch \
#<20240111><NoNeed>#     ${site_packages_dir}/jwst/outlier_detection/outlier_detection.py \
#<20240111><NoNeed>#     < ${script_dir}/dzliu_mirage_input_mosaic_image_fix_patch4.txt


if [[ ! -f ${site_packages_dir}/mirage/seed_image/crop_mosaic.py.ori ]]; then
    cp -i ${site_packages_dir}/mirage/seed_image/crop_mosaic.py \
          ${site_packages_dir}/mirage/seed_image/crop_mosaic.py.ori
fi

patch \
    ${site_packages_dir}/mirage/seed_image/crop_mosaic.py \
    < ${script_dir}/dzliu_mirage_input_mosaic_image_fix_patch5.txt


# 20230502: these patches are still needed for jwst==1.10.0 and mirage==2.4.0
# 20240111: these patches are still needed for jwst==1.12.5 and mirage==2.4.0 (except for patch4 now)
