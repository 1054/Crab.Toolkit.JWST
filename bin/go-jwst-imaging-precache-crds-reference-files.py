#!/usr/bin/env python
#
"""
Download CRDS cache files.

By Daizhong Liu.

"""

import os, sys, re, datetime
assert os.environ["CRDS_PATH"] != ''
assert os.environ["CRDS_SERVER_URL"] != ''
#os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
#os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Import CRDS
import crds

# Import JWST pipeline
from jwst import datamodels
from jwst.pipeline import calwebb_detector1

# Import click
import click

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
@click.argument('jwst_uncal_file', type=click.Path(exists=True))
def main(jwst_uncal_file):

    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Check input
    
    # Print CRDS pipeline version
    logger.info('CRDS pipeline version: {}'.format(crds.__version__))
    
    # 
    pipeline_context = crds.client.get_default_context('jwst')
    logger.info('pipeline_context: {}'.format(pipeline_context))
    
    #crds.get_cached_mapping(pipeline_context) # need local file to exist
    
    #crds.rmap.load_mapping(pipeline_context)
    
    # with datamodels.open(jwst_uncal_file) as model:
    #     params = {
    #         'INSTRUME': model.meta.instrument.name, 
    #         'DATE': model.meta.date, 
    #         'TIME': model.meta.observation.time,
    #     }
    #     crds.getreferences(
    #         parameters = parameters, 
    #         reftypes = ['DARK'], 
    #         context = pipeline_context,
    #         ignore_cache = False,
    #         observatory = 'jwst',
    #     )
    
    logger.info('Detector1Pipeline._precache_references: {!r}'.format(jwst_uncal_file))
    pipeline_object = calwebb_detector1.Detector1Pipeline()
    pipeline_object._precache_references(jwst_uncal_file)




# Main
if __name__ == '__main__':
    
    main()



