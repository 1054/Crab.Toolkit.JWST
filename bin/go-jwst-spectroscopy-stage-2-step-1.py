#!/usr/bin/env python
#
"""
Process a JWST rate data, output cal data. 

This script runs the JWST calwebb_spec2 pipeline.

Inputs
    
    jw*_rate.fits

Outputs

    jw*_cal.fits
    jw*_s2d.fits
    
Last update:
    
    2022-12-24 DZLIU.

Notes:
    

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
from jwst.pipeline import calwebb_spec2
from jwst import datamodels
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

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
@click.command()
@click.argument('input_rate_file', type=click.Path(exists=True))
@click.argument('output_cal_file', type=click.Path(exists=False))
@click.option('--sourcecat', 'input_sourcecat_file', type=click.Path(exists=True))
@click.option('--segmap', 'input_segmap_file', type=click.Path(exists=True))
@click.option('--direct-image', 'input_direct_image_file', type=click.Path(exists=True))
@click.option('--user-flat-file', type=click.Path(exists=True), default=None, help='A user-specified flat file, FITS format.')
@click.option('--user-flat-dir', type=click.Path(exists=True), default=None, help='A directory to find user-specified flat file e.g. "*{filter}*.fits".')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_rate_file, 
        output_cal_file, 
        input_sourcecat_file, 
        input_segmap_file, 
        input_direct_image_file, 
        user_flat_file, 
        user_flat_dir, 
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
        
    try:
        logger.info("CRDS_CONTEXT: {}".format(os.environ['CRDS_CONTEXT']))
    except KeyError:
        logger.error("Error! CRDS_CONTEXT environment variable not set!")
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
    
    
    # prepare asn file
    asn_file = os.path.join(os.path.dirname(output_filepath), 'asn_for_calwebb_spec2.json')
    asn_items = [(os.path.abspath(input_rate_file), 'science'),
                 (os.path.abspath(input_sourcecat_file), 'sourcecat'),
                 (os.path.abspath(input_segmap_file), 'segmap'),
                 (os.path.abspath(input_direct_image_file), 'direct_image'),
                ]
    asn_dict = asn_from_list(items, rule=DMSLevel2bBase)
    
    
    # prepare to run
    pipeline_object = calwebb_spec2.Spec2Pipeline()
    pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    
    
    # run
    run_output = pipeline_object.run(asn_file)
    
    
    # Check
    assert os.path.isfile(output_filepath)
    
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()



