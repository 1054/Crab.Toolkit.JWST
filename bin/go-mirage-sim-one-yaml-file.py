#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = '/n17data/dzliu/Data/jwst_crds_cache'
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = '/n23data1/hjmcc/jwst/mirage/mirage_data'
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

import click
import astropy.units as u
import numpy as np
#import pysiaf
import photutils # used by mirage.seed_image.fits_seed_image
#from photutils.psf.matching import TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow # see "mirage/seed_image/fits_seed_image.py", this is a bug in mirage!
for key in 'TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow'.split(','): # workaround for the above bug
    setattr(photutils, key.strip(), getattr(photutils.psf.matching, key.strip()))
import jwst
import mirage # pip install --upgrade git+https://github.com/spacetelescope/mirage.git
from astropy.coordinates import SkyCoord, FK5
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from functools import partial
from jwst.stpipe import Step
from jwst.outlier_detection import outlier_detection, OutlierDetectionStep
#outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, outlier_detection.OutlierDetection)
#outlier_detection.OutlierDetection.output_ext = 'fits'
#outlier_detection.OutlierDetection.make_output_path = OutlierDetectionStep.make_output_path
from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.catalogs import catalog_generator
from mirage.dark import dark_prep
from mirage.imaging_simulator import ImgSim
from mirage.ramp_generator import obs_generator
from mirage.seed_image import catalog_seed_image, blot_image
from mirage.seed_image.fits_seed_image import ImgSeed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirage-sim-mosaic')


DEFAULT_NEW_PSF_LIBRARY = '/automnt/n17data/dzliu/Work/JWST-MIRAGE-Simulation/20221008_make_mirage_psf_library/mirage_data/nircam/test_webbpsf_library'




####################
### MAIN PROGRAM ###
####################

@click.command()
@click.argument(yaml_file, type=click.Path(exists=True))
@click.option('--output-dir', type=click.Path(exists=False), default=None, help='Can specifiy a new output directory.')
@click.option('--output-file', type=str, default=None, help='Can specifiy a new output filename.')
@click.option('--new-psf-library', type=click.Path(exists=True), default=DEFAULT_NEW_PSF_LIBRARY, help='New psf library directory.')
@click.option('--overwrite', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        yaml_file, 
        output_dir, 
        output_file, 
        new_psf_library,
        overwrite, 
        verbose, 
    ):
    
    if verbose:
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    # get sim_data_filepath
    with open(yamlfile, 'r') as fp:
        yamldict = yaml.safe_load(fp)
        if output_dir is None:
            output_dir = yamldict['Output']['directory']
        if output_file is None:
            output_file = yamldict['Output']['file']
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    sim_data_filepath = os.path.join(output_dir, output_file)

    # Check sim data output file
    if os.path.isfile(sim_data_filepath):
        if overwrite_simdata:
            shutil.move(sim_data_filepath, sim_data_filepath+'.backup')
        else:
            if verbose:
                logger.info('Found existing data file {!r} and overwrite is set to False. Skipping.'.format(sim_data_filepath))
                return
    
    # Update yaml file in the output directory, we will use that instead of the input yaml file for the simulation
    original_yaml_file = yaml_file
    new_yaml_file = os.path.join(output_dir, os.path.basename(yaml_file))
    
    # Update yaml dict
    yamldict['Output']['directory'] = output_dir
    yamldict['Output']['file'] = output_file
    
    if new_psf_library is not None:
        yamldict['simSignals']['psfpath'] = new_psf_library
    
    # Write yaml dict to disk
    if os.path.isfile(new_yaml_file):
        shutil.move(new_yaml_file, new_yaml_file+'.backup')
    with open(new_yaml_file, 'w') as fp:
        yaml.dump(yamldict, fp)
    if verbose:
        logger.info('Updated {!r} with new PSF library {!r}'.format(output_yamlfile, new_psf_library))
    
    # Run simulation to generate data file
    if verbose:
        logger.info('*'*100)
        logger.info('*** Running Actual Simulation for {!r}'.format(
            new_yaml_file).ljust(96) + ' ***')
        logger.info('*'*100)
    m = ImgSim() # override_dark=dark
    m.paramfile = new_yaml_file
    m.create()




if __name__ == '__main__':
    main()



