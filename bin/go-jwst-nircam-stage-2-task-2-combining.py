#!/usr/bin/env python
#
"""
Process JWST rate data under "calibrated1_rates/" and named as "jw*_rate.fits"
into calibrated FITS data under "calibrated2_cals/" and named as "jw*_cal.fits". 

The calwebb_image2 pipeline: Calibrated Slope Images

    The Stage 2 calwebb_image2 pipeline applies instrumental corrections and calibrations 
    to the slope images output from Stage 1. This includes background subtraction, the 
    creation of a full World Coordinate System (WCS) for the data, application of the 
    flat field, and flux calibration. In most cases the final output is an image in 
    units of surface brightness. Whereas the input files had suffixes of *_rate.fits*, 
    the output files have suffixes of *_cal.fits*.

    In addition to the steps above, by default the Stage 2 pipeline will also run the 
    Resample step on the calibrated images, in order to remove the effects of instrument 
    distortion. This step outputs files with the suffix *_i2d.fits* that contain "rectified" 
    images. However, these files are meant only for user examination of the data. It is 
    the *_cal.fits* files that are passed on to Stage 3 of the pipeline.

    All JWST imaging mode data, regardless of instrument, are processed through the 
    calwebb_image2 pipeline. The steps and the order in which they are performed is the 
    same for all data. See Figure 1 on the calwebb_image2 algorithm page for a map of the 
    steps are performed on the input data.

Inputs

    A 2D countrate image (*_rate.fits) in units of DN/sec. The user can input a single 
    image file or an association file listing several files, in which case the processing 
    steps will be applied to each input exposure, one at a time.
Outputs

    A 2D calibrated, but unrectified, exposure (*_cal.fits) in units of MJy/sr
    A 2D resampled, or rectified, image (*_i2d.fits) in units of MJy/sr
    Note: At this stage, the resampled *_i2d.fits images are intended for quick-look use 
    only, while the *_cal.fits files are passed through for Stage 3 processing. We have 
    chosen to skip ResampleStep of calwebb_image2 to save on both processing time and disk 
    space. If you wish to perform this step to inspect the outputs, change the skip: true 
    in the jwst.resample.resample_step.ResampleStep dictionary of the image2_edited.asdf 
    parameter file (line 136) to skip: false. Alternatively comment out the line 
    image2.resample.skip = True in the cell using the run() method.
    
From ceers_nircam_reduction.ipynb

"""

# Packages that allow us to get information about objects:
import os, sys, re, json, copy, datetime, time, glob, shutil
import asdf
from collections import OrderedDict

# Numpy library:
import numpy as np

# For downloading data
import requests

# Astropy tools:
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch

# matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300

# Import JWST pipeline-related modules

# List of possible data quality flags
from jwst.datamodels import dqflags

# The entire calwebb_image2 pipeline
from jwst.pipeline import calwebb_image2

# Individual steps that make up calwebb_image2
from jwst.background import BackgroundStep
from jwst.assign_wcs import AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleStep
from jwst import datamodels

# importing an individual pipeline step
from jwst.skymatch import SkyMatchStep

# Import jwst package itself
import jwst

# Setup logging
import logging

# Define utility functions
def get_script_dir():
    """Get current script file's directory path."""
    return os.path.abspath(os.path.dirname(__file__))

def get_script_name():
    """Get current script file name without the suffix and replaced some characters to underscores."""
    return re.sub(r'[^a-zA-Z0-9_]', r'_', os.path.splitext(os.path.basename(__file__))[0])

def setup_logger():
    logger_streamhandler = logging.StreamHandler()
    logger_streamhandler_formatter = logging.Formatter("[%(asctime)-8s] %(message)s", "%H:%M:%S")
    logger_streamhandler.setFormatter(logger_streamhandler_formatter)
    logger_streamhandler.setLevel(logging.DEBUG)

    log_file = get_script_name()
    log_time = datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    logger_filehandler = logging.FileHandler(f"log_{log_file}_{log_time}.txt", mode='a')
    logger_filehandler_formatter = logging.Formatter("[%(asctime)-15s] %(message)s", "%Y-%m-%d %H:%M:%S")
    logger_filehandler.setFormatter(logger_filehandler_formatter)

    logger = logging.getLogger()
    while len(logger.handlers) > 0:
        del logger.handlers[0]
    logger.addHandler(logger_streamhandler)
    logger.addHandler(logger_filehandler)
    logger.setLevel(logging.DEBUG)
    
    return logger


# Main 
if __name__ == '__main__':

    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Print JWST pipeline version
    logger.info('JWST pipeline version: {}'.format(jwst.__version__))

    # check CRDS 
    try:
        logger.info("CRDS_PATH: {}".format(os.environ['CRDS_PATH']))
    except KeyError:
        logger.error("Error! CRDS_PATH environment variable not set!")
        sys.exit(-1)
        
    try:
        logger.info("CRDS_SERVER_URL: {}".format(os.environ['CRDS_SERVER_URL']))
    except KeyError:
        logger.error("Error! CRDS_SERVER_URL environment variable not set!")
        sys.exit(-1)
    
    
    # Find all NIRCam "jw*/calibrated1_rates"
    input_files = []
    input_files.extend(glob.glob("jw*nrca1/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrca2/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrca3/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrca4/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcb1/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcb2/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcb3/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcb4/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcalong/calibrated1_rates/jw*_rate.fits"))
    input_files.extend(glob.glob("jw*nrcblong/calibrated1_rates/jw*_rate.fits"))
    
    
    # find files to process
    if (len(input_files) == 0):
        logger.error("Error! No input file \"jw*nrc*/calibrated1_rates/jw*_rate.fits\" is found!")
        sys.exit(-1)
    input_files = sorted(input_files)
    
    
    # loop individual files, find unique 'obs_id' and 'target_name'
    list_of_obs_id = []
    list_of_target_name = []
    list_of_target_RA = []
    list_of_target_Dec = []
    list_of_instrument = []
    list_of_ifilter = []
    for input_filepath in input_files:
        logger.info("Processing {} -> {}".format(input_filepath, output_filepath))
        
        # read fits header
        header = fits.getheader(input_filepath, 0)
        obs_id = header['OBSERVTN'].strip()
        target_name = header['TARGPROP'].strip()
        target_RA = header['TARG_RA']
        target_Dec = header['TARG_DEC']
        instrument = header['INSTRUME']
        ifilter = header['FILTER']
        list_of_obs_id.extend(obs_id)
        list_of_target_name.extend(target_name)
        list_of_target_RA.extend(target_RA)
        list_of_target_Dec.extend(target_Dec)
        list_of_instrument = 
        list_of_ifilter = 
        
        # log
        logger.info("Processed {} -> {}".format(input_filepath, output_filepath))
    
    
    # get unique obs_id
    list_unique_obs_id = list(set(sorted(list_of_obs_id)))
    
    
    # loop unique obs_id
    for unique_obs_id in list_unique_obs_id:
        
    
    # prepare association file to process all rate files into one single output file
    if len(list_unique_obs_id) > 1:
        # 
        # TODO: Need to find groups
        #     "jw<sci>_<group>_<scan>_<instr>_(uncal|rate|cal).fits"
        # 
        # The Stage 2 pipeline can be called on a single fits file, or a collection of fits files. 
        # When calling on multiple files, the input is a json-formatted file called an "association" 
        # file that lists each of the fits files to be processed.
        asn_dict = OrderedDict()
        asn_dict['asn_type'] = 'None'
        asn_dict['asn_rule'] = 'DMS_Level2_Base'
        asn_dict['version_id'] = None
        asn_dict['code_version'] = jwst.__version__
        asn_dict['degraded_status'] = 'No known degraded exposures in association.'
        asn_dict['program'] = 'noprogram'
        asn_dict['constraints'] = 'No constraints'
        asn_dict['asn_id'] = obs_id
        asn_dict['asn_pool'] = 'none'
        asn_dict['products'] = []
        for input_file, output_file in list(zip(input_files, output_files)):
            input_filepath = os.path.join(input_dir, input_file)
            product_dict = OrderedDict()
            product_dict['name'] = input_file.replace(f"{input_suffix}.fits", "")
            product_dict['members'] = [
                    {'expname': input_filepath,
                     'exptype': 'science'}
                ]
            asn_dict['products'].append(product_dict)
        
        combined_name = 'level2_combined'
        if re.match(r'(jw[0-9]+)_([0-9]+)_([0-9]+)_([a-zA-Z0-9]+)_rate.fits', input_files[0]):
            combined_name = re.sub(r'(jw[0-9]+)_([0-9]+)_([0-9]+)_([a-zA-Z0-9]+)_rate.fits', r'\1_\2_\4_combined', input_files[0])
        
        asn_file = os.path.join(output_dir, f'{combined_name}_asn.json')
        
        if os.path.isfile(asn_file):
            shutil.move(asn_file, asn_file+'.backup')
        
        with open(asn_file, 'w') as fp:
            json.dump(asn_dict, fp, indent=4)
        
        # prepare a single output file 
        output_file = f"{combined_name}_cal.fits"
        output_filepath = os.path.join(output_dir, output_file)
        
        logger.info("Processing {} -> {}".format(input_files, output_filepath))
        
        # prepare to run
        pipeline_object = calwebb_image2.Image2Pipeline()
        pipeline_object.output_dir = output_dir
        pipeline_object.output_file = os.path.splitext(output_file)[0]
        pipeline_object.output_ext = ".fits"
        pipeline_object.save_results = True
        #pipeline_object.resample.skip = True # turn off the ResampleStep, comment out to produce the individual rectified *_i2d.fits for quick-look checks
        pipeline_object.resample.pixfrac = 1.0 # default
        pipeline_object.bkg_subtract.sigma = 3.0 # default
        
        # run
        run_output = pipeline_object.run(input_filepath)
        
        # check
        assert os.path.isfile(output_filepath)
        
        # log
        logger.info("Processed {} -> {}".format(input_files, output_filepath))
    
    
    
    
    



