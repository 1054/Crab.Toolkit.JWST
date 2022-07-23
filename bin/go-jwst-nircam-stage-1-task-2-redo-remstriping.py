#!/usr/bin/env python
#
"""
Redo remstriping

"""

# Packages that allow us to get information about objects:
import os, sys, re, json, copy, datetime, time, glob, shutil
import asdf

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

# Set CRDS
if not ('CRDS_PATH' in os.environ):
    os.environ['CRDS_PATH'] = os.path.expanduser('~/jwst_crds_cache')
if not ('CRDS_SERVER_URL' in os.environ):
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
#if not ('CRDS_CONTEXT' in os.environ):
#    os.environ['CRDS_CONTEXT'] = 

# Import JWST pipeline-related modules

# List of possible data quality flags
from jwst.datamodels import dqflags

# The entire calwebb_detector1 pipeline
from jwst.pipeline import calwebb_detector1

# Individual steps that make up calwebb_detector1
from jwst.dq_init import DQInitStep
from jwst.saturation import SaturationStep
from jwst.superbias import SuperBiasStep
from jwst.ipc import IPCStep                                                                                    
from jwst.refpix import RefPixStep                                                                
from jwst.linearity import LinearityStep
from jwst.persistence import PersistenceStep
from jwst.dark_current import DarkCurrentStep
from jwst.jump import JumpStep
from jwst.ramp_fitting import RampFitStep
from jwst import datamodels

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
    
    # read user input
    input_dirs = []
    input_files = []
    input_filepaths = []
    iarg = 1 
    while iarg < len(sys.argv):
        argstr = sys.argv[iarg]
        print("re.match(r'^(.*/|)jw.*_[0-9].*_[0-9].*_.*_.*$', '{}')".format(argstr))
        if re.match(r'^(.*/|)jw.*_[0-9].*_[0-9].*_.*_.*$', argstr):
            if os.path.isdir(argstr):
                list_files = os.path.listdir(argstr + os.sep + 'calibrated1_rates')
                for list_file in list_files:
                    if re.match(r'^jw.*_[0-9].*_[0-9].*_.*_rate.fits$', list_file):
                        input_dirs.append(argstr)
                        input_files.append(list_file)
                        input_filepaths.append(os.path.join(argstr, list_file))
            elif os.path.isfile(argstr):
                input_dirs.append(os.path.dirname(argstr))
                input_files.append(os.path.basename(argstr))
                input_filepaths.append(argstr)
        iarg+=1
    if len(input_dirs) == 0:
        print('Please input directories like "jw.*_[0-9].*_[0-9].*_[0-9].*_.*"  or files like "jw.*_[0-9].*_[0-9].*_[0-9].*_rate.fits"')
        sys.exit(255)
    

    # check CRDS, needed by remstriping.py
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
    
    
    # loop individual files
    for input_file, input_filepath in zip(input_files, input_filepaths):
        
        output_filepath = input_filepath
        output_dir = os.path.dirname(output_filepath)
        logger.info("Processing {} -> {}".format(input_filepath, output_filepath))
        
        # get fits header, check if NIRCam image, do 'remstriping'
        header = fits.getheader(output_filepath, 0)
        if header['INSTRUME'].strip().upper() == 'NIRCAM':
            
            # set CRDS_CONTEXT
            if ('CRDS_CONTEXT' not in os.environ) or (os.environ['CRDS_CONTEXT'] == ''):
                os.environ['CRDS_CONTEXT'] = header['CRDS_CTX']
            
            # additionally, following CEERS, 
            # measure and remove the horizontal and vertical striping from the two countrate images
            # -- see ceers_nircam_reduction.ipynb
            
            # prepare a temporary file
            temp_filepath = os.path.splitext(output_filepath)[0] + '_remstriping_rate.fits'
            backup_filepath = os.path.splitext(output_filepath)[0] + '_orig.fits'
            print('shutil.copy2({!r}, {!r})'.format(output_filepath, temp_filepath))
            shutil.copy2(output_filepath, temp_filepath)
            
            # run dzliu tool to make seed image, i.e., initial source masking
            from util_make_seed_image_for_rate_image import make_seed_image_for_rate_image
            print('make_seed_image_for_rate_image({!r}, overwrite=False)'.format(temp_filepath))
            make_seed_image_for_rate_image(temp_filepath, overwrite=False)
            
            # run CEERS team's script to remove striping
            from remstriping import measure_striping
            print('measure_striping({!r}, apply_flat=True, mask_sources=True, seedim_directory={!r}, threshold=0.0)'.format(temp_filepath, output_dir))
            measure_striping(temp_filepath, apply_flat=True, mask_sources=True, seedim_directory=output_dir, threshold=0.0)
            
            # make sure output fits file has the correct size -- NOT REALLY USEFUL BECAUSE THE CODE STOPPED AT `measure_striping`
            assert os.path.isfile(temp_filepath)
            # header_out = fits.getheader(temp_filepath, 1)
            # if os.path.getsize(temp_filepath) < len(header_out.tostring()) + header_out['NAXIS1']*header_out['NAXIS2']*np.abs(header_out['BITPIX']//8):
            #     shutil.move(temp_filepath, temp_filepath+'.corrupted')
            #     raise Exception('Error! The output file "{0}" has a too small size! Corrupted? Renaming it as "{0}.corrupted"!'.format(temp_filepath))
            
            # backup output_filepath and copy temporary file to output_filepath # do not backup it earlier as `measure_striping` may fail.
            shutil.copy2(output_filepath, backup_filepath)
            shutil.copy2(temp_filepath, output_filepath)
        
        
        # log
        logger.info("Processed {} -> {}".format(input_filepath, output_filepath))
    



