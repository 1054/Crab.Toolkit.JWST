patch \
$CONDA_PREFIX/lib/python3.9/site-packages/mirage/seed_image/fits_seed_image.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_mirage_input_mosaic_image_fix_patch1.txt

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/mirage/seed_image/blot_image.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_mirage_input_mosaic_image_fix_patch2.txt

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/mirage/seed_image/save_seed.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_mirage_input_mosaic_image_fix_patch3.txt

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_mirage_input_mosaic_image_fix_patch4.txt

patch \
$CONDA_PREFIX/lib/python3.9/site-packages/mirage/seed_image/crop_mosaic.py \
< $(dirname ${BASH_SOURCE[0]})/dzliu_mirage_input_mosaic_image_fix_patch5.txt


# 20230502: these patches are still needed for jwst==1.10.0 and mirage==2.4.0
