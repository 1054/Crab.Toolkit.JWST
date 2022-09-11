#!/usr/bin/env python
#
"""
Subtract background for a cal image using jwst.calwebb_image2.bkg_subtract. 

This requires some --dark-obs rate data.

"""

# Packages that allow us to get information about objects:
import os, sys, re, json, copy, datetime, time, glob, shutil
import click
from astropy.io import fits

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

# Import stdatamodels package
import stdatamodels

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
@click.option('--dark-obs', '--darkobs', 'input_dark_files', type=click.Path(exists=True), multiple=True, required=True)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_rate_file, 
        output_cal_file, 
        input_dark_files, 
        overwrite, 
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
    
    # Check output_filepath existence or history
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
    
    
    # read program info
    header = fits.getheader(input_filepath, 0)
    program = header['PROGRAM'].strip()
    obs_id = header['OBSERVTN'].strip()
    target_group = header['TARGPROP'].strip()
    target_name = header['TARGNAME'].strip()
    if target_name == '':
        target_name = target_group
    
    
    # prepare association content
    asn_dict = OrderedDict()
    asn_dict['asn_type'] = 'None'
    asn_dict['asn_rule'] = 'DMS_Level3_Base'
    asn_dict['version_id'] = None
    asn_dict['code_version'] = jwst.__version__
    asn_dict['degraded_status'] = 'No known degraded exposures in association.'
    asn_dict['program'] = program # 'noprogram'
    asn_dict['constraints'] = 'No constraints'
    asn_dict['asn_id'] = obs_id
    asn_dict['target'] = target_name
    asn_dict['asn_pool'] = 'none'
    asn_dict['products'] = []
    product_dict = OrderedDict()
    product_dict['name'] = output_filename
    product_dict['members'] = [
            {'expname': input_filepath, 
             'exptype': 'science'
            }
        ]
    for input_dark_file in input_dark_files:
        product_dict['members'].append(
            {'expname': input_dark_file, 
             'exptype': 'background'
            }
        )
    asn_dict['products'].append(product_dict)
    
    asn_filepath = os.path.join(output_dir, output_filename + '_asn.json')
    
    if os.path.isfile(asn_filepath):
        shutil.move(asn_filepath, asn_filepath+'.backup')
    
    with open(asn_filepath, 'w') as fp:
        json.dump(asn_dict, fp, indent=4)
    
    
    # prepare to run
    pipeline_object = calwebb_image2.Image2Pipeline()
    pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    pipeline_object.save_calibrated_ramp = True #<dzliu>#

    pipeline_object.bkg_subtract.save_combined_background = True
    
    
    # run
    run_output = pipeline_object.run(asn_filepath)
    
    
    # check output file existence
    assert os.path.isfile(output_filepath)
    
    
    # add history entry
    with datamodels.open(output_filepath) as model:
        timeobj = datetime.datetime.now()
        historystr = 'Done background subtraction; go-jwst-imaging-stage-2-step-3-redo-bkgsub.py' + \
                     '; time = ' + timeobj.strftime('%Y-%m-%d %H:%M:%S')
        historyobj = stdatamodels.util.create_history_entry(historystr)
        model.history.append(historyobj)
        
        model.save(
            output_filepath,
            overwrite=True,
        )
    
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()



