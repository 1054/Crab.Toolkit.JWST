
## Problematic datasets: 01180 F090W 007+010+018 nrca3

Copy rate files:

```
mkdir making_custom_wisps_jw01180_F090W_nrca3
cd making_custom_wisps_jw01180_F090W_nrca3
find ../processing_NIRCam_Imaging -maxdepth 3 -mindepth 3 -name "jw*_nrca3_rate.fits" -print0 | xargs -0 -I % cp % ./
```

Draw a custom DS9 region:

```
ds9 jw01180007001_02101_00001_nrca3_rate.fits

cat << EOF > ds9_nrca3_avoiding_source_detection_region.reg
image
polygon(0.21155622,1097.5957,34.047462,1097.2471,66.296089,1079.3312,109.29426,1000.5012,173.79151,846.42446,198.87378,674.43178,191.70742,527.52137,173.79151,402.11004,145.12607,330.44643,84.211993,287.44826,70.786398,258.81889,86.338396,213.8909,122.62639,170.6909,158.91439,115.39491,186.56239,49.730913,191.94514,0.10802455,0.31116566,-0.072810003,-0.19642985,753.11521)
EOF
```

Mask out real sources:

```
util_detect_source_and_create_seed_image.py jw01180007001_02101_00001_nrca3_rate.fits \
    --exclude-region ds9_nrca3_avoiding_source_detection_region.reg \
    --overwrite --median-filter 3 --smooth-after 1 --smooth-cutoff 0.3

find . -maxdepth 1 -mindepth 1 -name "jw*_rate.fits" -print0 | \
    while IFS= read -r -d $'\0' rate_file; do \
    echo $rate_file; \
    util_detect_source_and_create_seed_image.py "$rate_file" \
        --exclude-region ds9_nrca3_avoiding_source_detection_region.reg \
        --overwrite --median-filter 3 --smooth-after 1 --smooth-cutoff 0.3; \
    done

ds9 -colorbar no -lock frame image -scale zscale *_unmasked.fits

```


Combine unmasked wisps:

```
util_merge_source_emission_masked_rate_data.py jw01180{007,010,018}*_nrca3_rate_galaxy_seed_image_unmasked.fits combined.fits
```


Try removing wisps:

```
util_remove_wisps_with_templates.py \
    jw01180007001_02101_00001_nrca3_rate.fits \
    jw01180007001_02101_00001_nrca3_rate_try_removing_wisps.fits \
    --template-file combined.fits

ds9 -colorbar no -lock frame image -lock scalelimits yes -scale zscale \
    jw01180007001_02101_00001_nrca3_rate.fits \
    jw01180007001_02101_00001_nrca3_rate_try_removing_wisps.fits
    

```


Deploy!

```
echo $NIRCAM_WISP_TEMPLATES

cp -i combined.fits $NIRCAM_WISP_TEMPLATES/wisps_nrca3_F090W.fits

cat << EOF > $NIRCAM_WISP_TEMPLATES/wisps_nrca3_F090W.custom.readme.txt
made by Daizhong Liu
see https://github.com/1054/Crab.Toolkit.JWST/notes/20230921_01180_custom_wisps_rate_files.md
EOF

```






## Problematic datasets: 01180 F090W nrcb4 007+010+018

Copy rate files:

```
mkdir making_custom_wisps_jw01180_F090W_nrcb4
cd making_custom_wisps_jw01180_F090W_nrcb4
find ../processing_NIRCam_Imaging -maxdepth 3 -mindepth 3 -name "jw*_nrcb4_rate.fits" -print0 | xargs -0 -I % cp % ./
```

Draw a custom DS9 region:

```
ds9 jw01180007001_02101_00001_nrcb4_rate.fits

cat << EOF > ds9_nrcb4_avoiding_source_detection_region.reg
image
polygon(0.78702591,171.28042,0.42868472,686.36638,123.11675,432.35025,231.9808,249.18216,435.88491,-1.3779713,173.63472,-1.5672784,0.78702591,-1.5672784)
EOF
```

Mask out real sources:

```
util_detect_source_and_create_seed_image.py jw01180007001_02101_00001_nrcb4_rate.fits \
    --exclude-region ds9_nrcb4_avoiding_source_detection_region.reg \
    --overwrite --median-filter 4 --smooth-after 1 --smooth-cutoff 0.35 --box-frac 0.02

ds9 -colorbar no -lock frame image -scale zscale jw01180007001_02101_00001_nrcb4_rate.fits jw01180007001_02101_00001_nrcb4_rate_galaxy_seed_image_unmasked.fits jw01180007001_02101_00001_nrcb4_rate_galaxy_seed_image_zeroonemask.fits jw01180007001_02101_00001_nrcb4_rate_galaxy_seed_image_bkg2d.fits  jw01180007001_02101_00001_nrcb4_rate_galaxy_seed_image_rms2d.fits

find . -maxdepth 1 -mindepth 1 -name "jw*_rate.fits" -print0 | \
    while IFS= read -r -d $'\0' rate_file; do \
    echo $rate_file; \
    util_detect_source_and_create_seed_image.py "$rate_file" \
        --exclude-region ds9_nrcb4_avoiding_source_detection_region.reg \
        --overwrite --median-filter 4 --smooth-after 1 --smooth-cutoff 0.35 \
        --box-frac 0.02; \
    done

ds9 -colorbar no -lock frame image -scale zscale *_unmasked.fits

```


Combine unmasked wisps:

```
util_merge_source_emission_masked_rate_data.py jw01180{007,010,018}*_nrcb4_rate_galaxy_seed_image_unmasked.fits combined.fits
```


Try removing wisps:

```
util_remove_wisps_with_templates.py \
    jw01180007001_02101_00001_nrcb4_rate.fits \
    jw01180007001_02101_00001_nrcb4_rate_try_removing_wisps.fits \
    --template-file combined.fits

ds9 -colorbar no -lock frame image -lock scalelimits yes -scale zscale \
    jw01180007001_02101_00001_nrcb4_rate.fits \
    jw01180007001_02101_00001_nrcb4_rate_try_removing_wisps.fits

```


Deploy!

```
echo $NIRCAM_WISP_TEMPLATES

cp -i combined.fits $NIRCAM_WISP_TEMPLATES/wisps_nrcb4_F090W.fits

cat << EOF > $NIRCAM_WISP_TEMPLATES/wisps_nrcb4_F090W.custom.readme.txt
made by Daizhong Liu
see https://github.com/1054/Crab.Toolkit.JWST/notes/20230921_01180_custom_wisps_rate_files.md
EOF

```





