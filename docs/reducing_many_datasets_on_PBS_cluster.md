# Reducing many datasets on PBS cluster

## Aim

This is an exmaple of reducing many datasets on a PBS cluster, taking the JWST program 01345 as an example. A qsub script will be created to run parallel processing using JOBARRAY. 

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
go-jwst-download-by-proposal-id.py 01345 --calib-level 1 --download-dir $HOME/mast_download_dir
```

Launch your conda environment

```
conda activate YOUR_JWST_CONDA_ENV
```

Setup things using the script

```
go-jwst-imaging-stage-0
```

it will almost do every preparation as written in the [Reducing one dataset](reducing_one_dataset.md) document, then create a qsub script, then ask you if you want to submit that script to the cluster. You can check the source code of these scripts to see how it works or get inspired. 


## Last updates: 

- 2022-11-14 Daizhong Liu
