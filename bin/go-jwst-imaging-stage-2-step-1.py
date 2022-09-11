#!/usr/bin/env python
#
"""
Process a JWST rate data, output cal data. 

This script runs the JWST calwebb_image2 pipeline.

Inputs
    
    jw*_rate.fits

Outputs

    jw*_cal.fits
    jw*_i2d.fits
    
Last update:
    
    2022-09-10 DZLIU.


More notes from CEERS in "ceers_nircam_reduction.ipynb":

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

Notes:
    
    The Stage 2 pipeline can be called on a single fits file, or a collection of fits files. 
    When calling on multiple files, the input is a json-formatted file called an "association" 
    file that lists each of the fits files to be processed. 
    
    In this script we process each file separately without using an "association" file
    except for the additional sky subtraction step. 

"""

# Packages that allow us to get information about objects:
import os, sys, re, json, copy, datetime, time, glob, shutil
import click

# Set CRDS
if not ('CRDS_PATH' in os.environ):
    os.environ['CRDS_PATH'] = os.path.expanduser('~/jwst_crds_cache')
if not ('CRDS_SERVER_URL' in os.environ):
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
#if not ('CRDS_CONTEXT' in os.environ):
#    os.environ['CRDS_CONTEXT'] = 

# Import JWST pipeline-related modules

# The entire calwebb_image2 pipeline
from jwst.pipeline import calwebb_image2

# Individual steps that make up calwebb_image2
from jwst.background import BackgroundStep
from jwst.assign_wcs import AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleStep
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

def parse_jwst_data_name(input_str, raise_exception=True):
    regex_format = r'^jw([0-9]{5})([0-9]{3})([0-9]{3})_([0-9]{2})([0-9]{1})([0-9]{2})_([0-9]{5})_([a-z0-9]+)(.*)$'
    regex_match = re.match(regex_format, input_str)
    if regex_match is not None:
        return dict(zip(['program', 'obs_num', 'visit_num', 
                         'visit_group', 'parallel', 'activity', 
                         'exposure', 'detector', 'extra'], 
                        regex_match.groups()
               ))
    else:
        if raise_exception:
            raise Exception('Error! The input prefix does not seem to have the right format: {}'.format(regex_format))
        return None




# Main
@click.command()
@click.argument('input_rate_file', type=click.Path(exists=True))
@click.argument('output_cal_file', type=click.Path(exists=False))
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_rate_file, 
        output_cal_file, 
        overwrite = False, 
    ):
    
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
    
    # Set input_filepath, output_filepath
    input_filepath = input_rate_file
    output_filepath = output_cal_file
    output_filename = os.path.splitext(os.path.basename(output_filepath))[0] # no .fits suffix, has _cal suffix.
    output_dir = os.path.dirname(output_filepath)
    
    # Check output_filepath existence
    if os.path.isfile(output_filepath):
        if overwrite:
            shutil.move(output_filepath, output_filepath+'.backup')
        else:
            logger.info('Found existing output file {!r} and overwrite is False. Do nothing.'.format(
                output_filepath))
            return
    
    # Check output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    
    
    # Print progress
    logger.info("Processing {} -> {}".format(input_filepath, output_filepath))
    
    
    
    # prepare to run
    pipeline_object = calwebb_image2.Image2Pipeline()
    pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    pipeline_object.resample.skip = False # turn on the ResampleStep to produce the individual rectified *_i2d.fits for quick-look checks
    #pipeline_object.resample.pixfrac = 1.0 # default
    #pipeline_object.bkg_subtract.sigma = 3.0 # default
    
    
    
    # run
    run_output = pipeline_object.run(input_filepath)
    
    
    
    # Check
    assert os.path.isfile(output_filepath)
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()



