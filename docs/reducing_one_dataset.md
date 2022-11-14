# Reducing one dataset

## Aim

This is an exmaple of reducing one dataset, taking the dataset `jw01345002001_14201_00001_nrcb4` from program 01345 as an example. 

## How to

### Stage 0: Preparation

Download the toolkit into a local directory, for example `$HOME/Cloud/Github/Crab.Toolkit.JWST`.

Add the toolkit bin path into system PATH

```
export PATH=$PATH:$HOME/Cloud/Github/Crab.Toolkit.JWST/bin
```

Setup CRDS and NIRCam wisp templates. The latter should be obtained from [https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-features-and-caveats/nircam-claws-and-wisps](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-features-and-caveats/nircam-claws-and-wisps), and stored in a local directory, for example `$HOME/jwst_nircam_wisp_templates`. 

```
export CRDS_PATH=$HOME/jwst_crds_cache
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_CONTEXT=jwst_0995.pmap
export NIRCAM_WISP_TEMPLATES=$HOME/jwst_nircam_wisp_templates
```

Download JWST data into a local directory, for example `$HOME/mast_download_dir`

```
go-jwst-download-by-proposal-id.py 01345 --dataset jw01345002001_14201_00001_nrcb4 --calib-level 1 --download-dir $HOME/mast_download_dir
```

Setup per-dataset subdirectories in your working directory

```
go-jwst-imaging-precopy-uncal-files $HOME/mast_download_dir/mastDownload/JWST/jw*/jw*_uncal.fits
```

this will create a subdirectory `jw01345002001_14201_00001_nrcb4 ` under current directory, and the uncal file will be copied to `jw01345002001_14201_00001_nrcb4/uncals/jw01345002001_14201_00001_nrcb4_uncal.fits`. 

Precache CRDS files

```
go-jwst-imaging-precache-crds-reference-files
```

Presort visit gropus (optional, for parallel run on cluster)

```
go-jwst-imaging-presort-visits.py
```

Create qsub script (optional, for parallel run on cluster)

```
go-qsub-jwst-imaging-jobarray-group.bash
```


### Stage 1: uncal -> rate

```
go-jwst-imaging-stage-1
```
This will find all dataset directories under the current directory, then process the `{dataset}/uncals/{dataset}_uncal.fits` into `{dataset}/calibrated1_rates/{dataset}_rate.fits`.

This stage includes following steps:

1. Running `calwebb_detector1` pipeline. 
2. Following [Bagley+2022](https://arxiv.org/abs/2211.02495), applying flatfield correction using the jwst pipeline's `do_correction` function. 
3. Following [Bagley+2022](https://arxiv.org/abs/2211.02495), removing NIRCam F150W/F150W2/F200W/F210M wisps using NIRCam team's wisp templats (see preparation, need to set the system variable `NIRCAM_WISP_TEMPLATES` for our script).
4. Following [Bagley+2022](https://arxiv.org/abs/2211.02495), removing stripes (1/f noise), for NIRCam only. This step does more than [Bagley+2022](https://arxiv.org/abs/2211.02495) by destriping along two additional angles corresponding to the diffraction spikes of the PSF (of course with sufficient masking of source emission). This efficiently removes the tilted stripes heavily affecting PRIMER's NIRCam data.


### Stage 2: rate -> cal

```
go-jwst-imaging-stage-2
```
This will find all dataset directories under the current directory, then process the `{dataset}/calibrated1_rates/{dataset}_rate.fits` into `{dataset}/calibrated2_cals/{dataset}_cal.fits`.

This stage includes following steps:

1. Running `calwebb_image2` pipeline. 
2. Following [Bagley+2022](https://arxiv.org/abs/2211.02495), performing a constant background subtraction using the `skymatch` module of the jwst pipeline. 
3. _(The script for this step is used for later MIRI reprocessing.)_


### Stage 3: cal -> i2d (drizzling mosaic image)

```
go-jwst-imaging-stage-3
```
This will find all dataset directories under the current directory, group all `*/calibrated1_rates/*_cal.fits` by visit/instrument/filter, then generate a asn file for each group and run stage3 pipeline and output `calibrated3_mosaic/{obsnum}_{instrument}_{filter}/{dataset}_i2d.fits`.

This stage includes following steps:

1. Grouping cal.fits files, generating asn file, and running `calwebb_image3` pipeline to make drizzled mosaic image for each group. The output directory is `calibrated3_mosiac` under the current directory. A `list_obs.csv` table will be created therein. A subdirectory for each visit/intrument/filter group will be created therein, named like `{obsnum}_{instrument}_{filter}`. Images and temporary files will be in the subdirectories. 
2. _For MIRI Image-for-Image background subtraction only (if `--no-reprocess-MIRI` is not set)._ This step detects source emission in the final mosaic image, and uses that as the mask to mask out source emission in rate.fits files. The produced `*_masked_source_emission_rate.fits` files are like dark observations which will be used as the background. 
3. _For MIRI Image-for-Image background subtraction only (if `--no-reprocess-MIRI` is not set)._ This step merges the source-emission-masked rate files (`*_masked_source_emission_rate.fits`) made from other visits/exposures datasets that have the same filter/detector for each dataset. 
4. _For MIRI Image-for-Image background subtraction only (if `--no-reprocess-MIRI` is not set)._ Run `calwebb_Image2` pipeline with an asn file that includes both `expname` and `background` to do the Image-for-Image background subtraction. The output is `{dataset}/calibrated2_cals/`


### Stage 4: cal -> i2d (drizzling large mosaic image)

We can rerun the stage3 code with some options to make large mosaic image. Type `--help` to see the help content (run the stage-3-step-1 python code, not the stage-3 bash script):

```
go-jwst-imaging-stage-3-step-1.py --help
```

For example, to combine three visits and output to the directory `calibrated_mosaic_multiobs`: 

```
go-jwst-imaging-stage-3-step-1.py \
    jw01345001*/calibrated2_cals/*_cal.fits \
    jw01345002*/calibrated2_cals/*_cal.fits \
    jw01345003*/calibrated2_cals/*_cal.fits \
    --combine-obsnum \
    --pixel-scale 0.030 \
    --kernel square \
    --pixfrac 1.0 \
    calibrated_mosaic_multiobs
```


## Last updates: 

- 2022-11-14 Daizhong Liu








