# Making mirage simulation

## Aim

This is an exmaple of making mirage simulation for COSMOS-Web six December pointings. 

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

Download MIRAGE data from [https://mirage-data-simulator.readthedocs.io/en/latest/reference_files.html#reference-files](https://mirage-data-simulator.readthedocs.io/en/latest/reference_files.html#reference-files) to a local directory, for example, `$HOME/jwst_mirage_data` and set the system variable `MIRAGE_DATA`. 

```
export MIRAGE_DATA=$HOME/jwst_mirage_data
```

Prepare APT xml and pointing files. The need to be extracted from the *.aptx file using the JWST APT software [https://www.stsci.edu/scientific-community/software/astronomers-proposal-tool-apt](https://www.stsci.edu/scientific-community/software/astronomers-proposal-tool-apt). 

Put the APT xml and pointing files into an `apt_files` directory under the current directory

```
ls apt_files/*.xml
ls apt_files/*.xml
```

Prepare a star catalog and put that into a `input_catalogs` directory.

Okay then we can proceed to the next stage. 

## Stage 1: run the script

### To run a simulation using an input mosaic image plus a star catalog:

```
go-mirage-sim-mosaic.py \
    --xml-file 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.xml' \
    --pointing-file 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.pointing' \
    --mosaic-file input_mosaic_images/mosaic_image.fits \
    --star-catalog 'input_catalogs/ptsrc_pointings_BEST_sw_tot.cat' \
    --filter F200W \
    --dates '2023-01-01' \
    --pa-v3 293.09730273
```

The input mosaic image can be an intrisically no-PSF image, in which case a single-pixel-PSF will be used for the PSF-matching, or an image with a finite PSF, in which case please set `--mosaic-psf` with the PSF FWHM in arcsec or a PSF file (in the latter case, a Gaussian2D model will be fit to the PSF to estimate the FWHM, see `mirage/seed_image/fits_seed_image.py`). 

The output of the script are yaml files under the output directory `yaml_files`, and simulated data under the `sim_data` directory. Those `*_uncal.fits` files are simulated raw data files that will be processed by the JWST pipeline. 


## To run a simulation using an input galaxy catalog plus a star catalog:

```
go-mirage-sim-mosaic.py \
    --xml-file 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.xml' \
    --pointing-file 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.pointing' \
    --galaxy-catalog 'input_catalogs/galaxies.cat' \
    --star-catalog 'input_catalogs/ptsrc_pointings_BEST_sw_tot.cat' \
    --filter F200W \
    --dates '2023-01-01' \
    --pa-v3 293.09730273
```


## How MIRAGE handles flat field correction? 

It's already using the CRDS flat reference file for simulation! See its source code `mirage/ramp_generator/obs_generator.py`:

```
class Observation():
    ...
    def create(self, override_refs=None):
            ...
            # Multiply flat fields
            simexp = self.add_flatfield_effects(simexp)
            simzero = self.add_flatfield_effects(np.expand_dims(simzero, axis=1))[:, 0, :, :]
            ...
```

```
class Observation():
    ...
    def add_flatfield_effects(self, ramp):
        ...
        # PIXEL FLAT
        if self.runStep['pixelflat']:
            pixelflat, pixelflatheader = self.read_cal_file(self.params['Reffiles']['pixelflat'])
            ramp *= pixelflat
        ...
```


## How to make custom PSFs with MIRAGE? 

We need an old versoin of MIRAGE which has `psf_library.py` (can't find it now?), or webbpsf `gridded_library.py` [https://github.com/spacetelescope/webbpsf/blob/develop/webbpsf/gridded_library.py](https://github.com/spacetelescope/webbpsf/blob/develop/webbpsf/gridded_library.py). 

TBD



## Notes: 

- Require JWST pipeline version `>= 1.8.2`.
- How MIRAGE handles flat field correction? It's already using the CRDS flat reference file for simulation!



## Last updates: 

- 2022-11-14 Daizhong Liu








