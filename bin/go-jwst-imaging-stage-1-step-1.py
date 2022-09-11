#!/usr/bin/env python
#
"""
Process a JWST uncalibrated data, output rate data. 

This script runs the JWST calwebb_detector1 pipeline: Ramps to Slopes. 

Inputs
    
    jw*_uncal.fits

Outputs

    jw*_rate.fits
    jw*_rateints.fits
    
Last update:
    
    2022-09-10 DZLIU.

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
@click.argument('input_uncal_file', type=click.Path(exists=True))
@click.argument('output_rate_file', type=click.Path(exists=False))
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_uncal_file, 
        output_rate_file, 
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
    input_filepath = input_uncal_file
    output_filepath = output_rate_file
    if os.path.isfile(output_rate_file):
        if overwrite:
            shutil.move(output_rate_file, output_rate_file+'.backup')
        else:
            logger.info('Found existing output file {!r} and overwrite is False. Do nothing.'.format(
                output_filepath))
            return
    
    # Check output_dir
    if output_filepath.find(os.sep) >= 0:
        output_dir = os.path.dirname(output_filepath)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = './'
    
    # Print progress
    logger.info("Processing {} -> {}".format(input_filepath, output_filepath))
    
    # Prepare parameters to run the JWST pipeline
    pipeline_object = calwebb_detector1.Detector1Pipeline()
    pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    pipeline_object.ipc.skip = False # turn on IPCStep (default is already on)
    pipeline_object.persistence.skip = False # turn on PersistenceStep (default is already on)
    pipeline_object.persistence.maximum_cores = 'all' # 
    pipeline_object.jump.maximum_cores = 'all' # 
    pipeline_object.ramp_fit.maximum_cores = 'all' # 
    pipeline_object.save_calibrated_ramp = True #<dzliu>#
    
    # Run
    run_output = pipeline_object.run(input_filepath)
    
    # Check
    assert os.path.isfile(output_filepath)
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()


