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




## Last updates: 

- 2022-11-14 Daizhong Liu








