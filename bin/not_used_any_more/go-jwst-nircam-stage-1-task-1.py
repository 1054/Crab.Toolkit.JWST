#!/usr/bin/env python
#
"""
Process JWST uncalibrated FITS data under "uncals/" and named as "jw*_uncal.fits"
into calibrated FITS data under "calibrated1_rates/" and named as "jw*_rate.fits". 

The calwebb_detector1 pipeline: Ramps to Slopes

Inputs

    A raw exposure (*_uncal.fits) containing the 4-dimensional raw data from all detector readouts: (ncols x nrows x ngroups x nintegrations).

Outputs

    A 2D countrate image (*_rate.fits) resulting from averaging over the exposure's integrations.
    A 3D countrate image (*_rateints.fits) containing the results of each integration in separate extensions.
    
From ceers_nircam_reduction.ipynb

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
    
    # Read user input
    iarg = 1
    arg_str = ''
    do_remstripping = True
    while iarg < len(sys.argv):
        arg_str = sys.argv[iarg].lower().replace('--', '-')
        if arg_str == '-no-remstripping':
            print('User set no remstriping.')
            do_remstripping = False

    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Print JWST pipeline version
    logger.info('JWST pipeline version: {}'.format(jwst.__version__))

    # Define input and output dirs for this stage
    input_dir = "uncals"
    input_suffix = "_uncal"
    output_dir = "calibrated1_rates"
    output_suffix = "_rate"
    if not os.path.isdir(input_dir):
        logger.error("Error! Input directory \"{}\" does not exist!".format(input_dir))
        sys.exit(-1)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    # check override_gain_file
    override_gain_file = ""
    #override_gain_file = "gains_v2.1.0/jwst_nircam_gain_nrca1.fits" # TODO

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

    #try:
    #    print(os.environ['CRDS_CONTEXT'])
    #except KeyError:
    #    print('CRDS_CONTEXT environment variable not set!')

    # check stage asdf file
    #asdf_file = get_script_name() + ".asdf"
    #if not os.path.exists(asdf_file):
    #    logger.error("Error! The asdf file is not found \"{}\"!".format(asdf_file))
    #    sys.exit(-1)
        
    # print asdf file tree
    #asdf_obj = asdf.open(asdf_file)
    #logger.info("asdf file tree: \n" + str(asdf_obj.tree))
    #asdf_obj.close()
    
    
    # find files to process
    input_files = [t for t in os.listdir(input_dir) if t.endswith(f"{input_suffix}.fits")]
    if (len(input_files) == 0):
        logger.error("Error! No input file \"{}/*{}\" is found!".format(input_dir, f"{input_suffix}.fits"))
        sys.exit(-1)
    input_files = sorted(input_files)
    output_files = [t.replace(f"{input_suffix}.fits", f"{output_suffix}.fits") for t in input_files]
    
    
    # loop individual files
    for input_file, output_file in list(zip(input_files, output_files)):
        input_filepath = os.path.join(input_dir, input_file)
        output_filepath = os.path.join(output_dir, output_file)
        logger.info("Processing {} -> {}".format(input_filepath, output_filepath))
        
        # prepare to run
        pipeline_object = calwebb_detector1.Detector1Pipeline()
        pipeline_object.output_dir = output_dir
        pipeline_object.save_results = True
        pipeline_object.ipc.skip = False # turn on IPCStep (default is already on)
        pipeline_object.persistence.skip = False # turn on PersistenceStep (default is already on)
        pipeline_object.persistence.maximum_cores = 'all' # 
        pipeline_object.jump.maximum_cores = 'all' # 
        pipeline_object.ramp_fit.maximum_cores = 'all' # 
        pipeline_object.save_calibrated_ramp = True #<dzliu>#
        # Specify the name of the gain file that will override 
        # the existing gain reference file used for the jump and ramp_fit steps
        if override_gain_file is not None and override_gain_file != "":
            pipeline_object.jump.override_gain = override_gain_file
            pipeline_object.ramp_fit.override_gain = override_gain_file
        
        # run
        run_output = pipeline_object.run(input_filepath)
        
        # check
        assert os.path.isfile(output_filepath)
    
        
        # get fits header, check if NIRCam image, do 'remstriping'
        if do_remstripping:
            header = fits.getheader(output_filepath, 0)
            if header['INSTRUME'].strip().upper() != 'NIRCAM':
                do_remstripping = False
        
        if do_remstripping:
            # set CRDS_CONTEXT
            if ('CRDS_CONTEXT' not in os.environ) or (os.environ['CRDS_CONTEXT'] == ''):
                os.environ['CRDS_CONTEXT'] = header['CRDS_CTX']
            
            # additionally, following CEERS, 
            # measure and remove the horizontal and vertical striping from the two countrate images
            # -- see ceers_nircam_reduction.ipynb
            
            # prepare a temporary file
            temp_filepath = os.path.splitext(output_filepath)[0] + '_remstriping_rate.fits'
            backup_filepath = os.path.splitext(output_filepath)[0] + '_orig.fits'
            shutil.copy2(output_filepath, temp_filepath)
            
            # run dzliu tool to make seed image, i.e., initial source masking
            from util_make_seed_image_for_rate_image import make_seed_image_for_rate_image
            make_seed_image_for_rate_image(temp_filepath)
            
            # run CEERS team's script to remove striping
            from remstriping import measure_striping
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
    



