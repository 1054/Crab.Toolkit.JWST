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
from jwst.pipeline import calwebb_image2
from jwst import datamodels
from jwst.associations.asn_from_list import asn_from_list
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
@click.option('--sourcecat', 'input_sourcecat_file', type=click.Path(exists=True), default=None)
@click.option('--segmap', 'input_segmap_file', type=click.Path(exists=True), default=None)
@click.option('--direct-image', 'input_direct_image_file', type=click.Path(exists=True), default=None)
@click.option('--msa-meta-file', 'input_msa_meta_file', type=click.Path(exists=True), default=None, help='The auxiliary data shipped with the NRS_MSASPEC data.')
@click.option('--user-flat-file', type=click.Path(exists=True), default=None, help='A user-specified flat file, FITS format.')
@click.option('--user-flat-dir', type=click.Path(exists=True), default=None, help='A directory to find user-specified flat file e.g. "*{filter}*.fits".')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_rate_file, 
        output_cal_file, 
        input_sourcecat_file, 
        input_segmap_file, 
        input_direct_image_file, 
        input_msa_meta_file, 
        user_flat_file, #<TODO><20230915># not used
        user_flat_dir, #<TODO><20230915># not used
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
    
    
    # check WFSS
    with datamodels.open(input_rate_file) as model:
        exp_type = model.meta.exposure.type
    if exp_type in calwebb_spec2.WFSS_TYPES:
        # see calwebb_spec2.py `if exp_type in WFSS_TYPES:`
        asn_items = [(os.path.basename(input_rate_file), 'science'),
                     (os.path.basename(input_sourcecat_file), 'sourcecat'),
                     (os.path.basename(input_segmap_file), 'segmap'),
                     (os.path.basename(input_direct_image_file), 'direct_image'),
                    ]
        linked_rate_file = os.path.join(output_dir, os.path.basename(input_rate_file))
        linked_sourcecat_file = os.path.join(output_dir, os.path.basename(input_sourcecat_file))
        linked_segmap_file = os.path.join(output_dir, os.path.basename(input_segmap_file))
        linked_direct_image_file = os.path.join(output_dir, os.path.basename(input_direct_image_file))
        if os.path.isfile(linked_rate_file) or os.path.islink(linked_rate_file):
            os.remove(linked_rate_file)
        if os.path.isfile(linked_sourcecat_file) or os.path.islink(linked_sourcecat_file):
            os.remove(linked_sourcecat_file)
        if os.path.isfile(linked_segmap_file) or os.path.islink(linked_segmap_file):
            os.remove(linked_segmap_file)
        if os.path.isfile(linked_direct_image_file) or os.path.islink(linked_direct_image_file):
            os.remove(linked_direct_image_file)
        os.symlink(os.path.relpath(input_rate_file, output_dir), linked_rate_file)
        os.symlink(os.path.relpath(input_sourcecat_file, output_dir), linked_sourcecat_file)
        os.symlink(os.path.relpath(input_segmap_file, output_dir), linked_segmap_file)
        os.symlink(os.path.relpath(input_direct_image_file, output_dir), linked_direct_image_file)
    elif exp_type in ['NRS_MSATA', 'NRS_TACONFIRM']:
        asn_items = [(os.path.basename(input_rate_file), 'science'),
                    ]
        linked_rate_file = os.path.join(output_dir, os.path.basename(input_rate_file))
        if os.path.isfile(linked_rate_file) or os.path.islink(linked_rate_file):
            os.remove(linked_rate_file)
        os.symlink(os.path.relpath(input_rate_file, output_dir), linked_rate_file)
    elif exp_type == 'NRS_MSASPEC':
        if input_msa_meta_file is None:
            logger.info('Error! We need a --msa-meta-file!')
            raise Exception('Error! We need a --msa-meta-file!')
        asn_items = [(os.path.basename(input_rate_file), 'science'),
                    ]
        linked_rate_file = os.path.join(output_dir, os.path.basename(input_rate_file))
        linked_msa_meta_file = os.path.join(output_dir, os.path.basename(input_msa_meta_file))
        if os.path.isfile(linked_rate_file) or os.path.islink(linked_rate_file):
            os.remove(linked_rate_file)
        if os.path.isfile(linked_msa_meta_file) or os.path.islink(linked_msa_meta_file):
            os.remove(linked_msa_meta_file)
        os.symlink(os.path.relpath(input_rate_file, output_dir), linked_rate_file)
        os.symlink(os.path.relpath(input_msa_meta_file, output_dir), linked_msa_meta_file)
    else:
        raise NotImplementedError('Not implemented for the exp_type {!r} of data file {!r}!'.format(exp_type, input_rate_file))
    
    
    # prepare asn file
    asn_file = os.path.join(output_dir, 'asn_for_calwebb_spec2.json')
    asn_obj = asn_from_list(asn_items, 
        product_name=output_filename,
        with_exptype=True,
        #rule=DMSLevel2bBase,
    )
    #asn_obj.filename = asn_file
    _file_name, serialized = asn_obj.dump()
    if os.path.isfile(asn_file):
        shutil.move(asn_file, asn_file+'.backup')
    with open(asn_file, 'w') as fp:
        fp.write(serialized)
    
    
    # prepare to run
    pipeline_object = calwebb_spec2.Spec2Pipeline()
    #pipeline_object.output_dir = output_dir
    pipeline_object.save_results = True
    
    
    # run
    current_dir = os.getcwd()
    os.chdir(output_dir)
    run_output = pipeline_object.run(os.path.basename(asn_file))
    os.chdir(current_dir)
    
    
    # Check
    assert os.path.isfile(output_filepath)
    
    
    # Print progress
    logger.info("Processed {} -> {}".format(input_filepath, output_filepath))




# Entry
if __name__ == '__main__':
    main()



