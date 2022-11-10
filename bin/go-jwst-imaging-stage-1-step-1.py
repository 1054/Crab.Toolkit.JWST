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
import jwst
import stcal
import stdatamodels
from asdf import AsdfFile

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

# Setup logging
import logging

# Version control
from packaging import version

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
@click.command()
@click.argument('input_uncal_file', type=click.Path(exists=True))
@click.argument('output_rate_file', type=click.Path(exists=False))
@click.option('--maximum-cores', type=str, default='all')
@click.option('--remove-snowballs/--no-remove-snowballs', is_flag=True, default=True)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_uncal_file, 
        output_rate_file, 
        maximum_cores, 
        remove_snowballs, 
        overwrite, 
    ):
    
    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Print JWST pipeline version
    logger.info('JWST pipeline version: {}'.format(jwst.__version__))
    
    if version.parse(jwst.__version__) < version.parse('1.8.0'):
        logger.error("Error! Requiring JWST pipeline version newer than 1.8.0!") # for the jump step to deal with snowballs
        sys.exit(-1)
    
    # Print STCAL pipeline version
    logger.info('STCAL pipeline version: {}'.format(stcal.__version__))
    
    if version.parse(stcal.__version__) < version.parse('1.2.0'):
        logger.error("Error! Requiring STCAL pipeline version newer than 1.2.0!") # for the jump step to deal with snowballs
        sys.exit(-1)
    
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
        
    try:
        logger.info("CRDS_CONTEXT: {}".format(os.environ['CRDS_CONTEXT']))
    except KeyError:
        logger.error("Error! CRDS_CONTEXT environment variable not set!")
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
    
    # Get instrument
    instrument_name = ''
    with datamodels.open(input_filepath) as model:
        instrument_name = model.meta.instrument.name
    
    # Prepare parameters to run the JWST pipeline
    pipeline_object = calwebb_detector1.Detector1Pipeline()
    pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    pipeline_object.ipc.skip = False # turn on IPCStep (default is already on)
    pipeline_object.persistence.skip = False # turn on PersistenceStep (default is already on)
    pipeline_object.persistence.maximum_cores = maximum_cores
    pipeline_object.jump.maximum_cores = maximum_cores
    if remove_snowballs:
        pipeline_object.jump.min_jump_area = 100.0 # 5.0 is the default, minimum area to trigger large events processing
        pipeline_object.jump.expand_factor = 1.5 # 2.0 is the default, The expansion factor for the enclosing circles or ellipses
        pipeline_object.jump.sat_required_snowball = False
        pipeline_object.jump.expand_large_events = True
        if instrument_name.upper() == 'MIRI':
            pipeline_object.jump.use_ellipses = True # "stcal/jump/jump.py" says that "Use ellipses rather than circles (better for MIRI)"
    pipeline_object.ramp_fit.maximum_cores = maximum_cores
    pipeline_object.save_calibrated_ramp = True
    
    # Run
    run_output = pipeline_object.run(input_filepath)
    
    # Check
    assert os.path.isfile(output_filepath)
    
    # Save pars
    asdf_filepath = os.path.splitext(output_filepath)[0] + '_calwebb_detector1.asdf'
    if os.path.isfile(asdf_filepath):
        shutil.move(asdf_filepath, asdf_filepath+'.backup')
    asdf_object = AsdfFile(pipeline_object.get_pars())
    asdf_object.write_to(asdf_filepath)
    logger.info('Parameters are saved into {}'.format(asdf_filepath))
    
    # Write history about removed snowballs
    if remove_snowballs:
        with datamodels.open(output_filepath) as model:
            history_entry = stdatamodels.util.create_history_entry(
                'Removed snowballs with JWST pipeline {} JumpStep module (min_jump_area = {}, expand_factor = {}, sat_required_snowball = {}, expand_large_events = {})'.format(
                    jwst.__version__, pipeline_object.jump.min_jump_area, pipeline_object.jump.expand_factor, pipeline_object.jump.sat_required_snowball, pipeline_object.jump.expand_large_events
                )
            )
            model.history.append(history_entry)
            model.save(output_filepath, overwrite=True)
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()



