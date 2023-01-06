# Example of removing claws

## Step 1

First open DS9,

```
ds9 -lock frame image -scale zscale jw01727044001_04101_0000*_cal.fits
```

Draw ds9 regions and save them in pixel (image) coordinates as:

```
jw01727044001_04101_00001_nrcb1_cal_claws.pix.reg
jw01727044001_04101_00002_nrcb1_cal_claws.pix.reg
jw01727044001_04101_00003_nrcb1_cal_claws.pix.reg
jw01727044001_04101_00004_nrcb1_cal_claws.pix.reg
```


## Step 2 

Then we run the following code: 

(assuming you have the `/path/to/Crab.Toolkit.JWST/bin` path in your system variable `$PATH` so that commands below can be executed in the command line)

```
for expo in 00001 00002 00003 00004; do

    util_detect_source_and_create_seed_image.py \
        jw01727044001_04101_${expo}_nrcb1_cal.fits \
        --exclude-region jw01727044001_04101_${expo}_nrcb1_cal_claws.pix.reg \
        --minpixarea 4 --smooth 2 \
        --overwrite

    util_detect_source_and_create_seed_image.py \
        jw01727044001_04101_${expo}_nrcb1_cal.fits \
        --include-region jw01727044001_04101_${expo}_nrcb1_cal_claws.pix.reg \
        --output-suffix "_claws_seed_image" \
        --minpixarea 1 --ignore-background --sigma 2 --smooth 1 --overwrite

    ds9 -tile grid layout 3 2 -width 1200 -height 550 -lock frame image -scale zscale \
        jw01727044001_04101_${expo}_nrcb1_cal.fits \
        jw01727044001_04101_${expo}_nrcb1_cal_galaxy_seed_image.fits \
        jw01727044001_04101_${expo}_nrcb1_cal_galaxy_seed_image_zeroonemask.fits \
        jw01727044001_04101_${expo}_nrcb1_cal_claws_seed_image.fits \
        jw01727044001_04101_${expo}_nrcb1_cal_claws_seed_image_zeroonemask.fits \
        -zoom to fit \
        -saveimage eps jw01727044001_04101_${expo}_nrcb1_cal_claws.eps

    gzip -k jw01727044001_04101_${expo}_nrcb1_cal_claws_seed_image_zeroonemask.fits

    echo "Output mask file: jw01727044001_04101_${expo}_nrcb1_cal_claws_seed_image_zeroonemask.fits.gz"

done
```

## Afterwards

Then we can apply the claw mask to "rate.fits" files. For each rate file with a corresponding claw mask file, we can open the rate file with `jwst.datamodels.open()` then set `model.dq = np.bitwise_or(model.dq, claw_mask.astype(int))`. 


