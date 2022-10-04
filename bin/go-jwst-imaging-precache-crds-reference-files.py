#!/usr/bin/env python
#
"""
Download CRDS cache files.

By Daizhong Liu.

"""

import os, sys, re, datetime, glob
assert os.environ["CRDS_PATH"] != ''
assert os.environ["CRDS_SERVER_URL"] != ''
#os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
#os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Import CRDS
import crds

# Import JWST pipeline
from jwst import datamodels
from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2

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
@click.argument('jwst_uncal_files', nargs=-1, type=click.Path(exists=True))
def main(jwst_uncal_files):

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
    
    all_jwst_uncal_files = []
    
    for jwst_uncal_file in jwst_uncal_files:
        if jwst_uncal_file.find('*') >= 0:
            all_jwst_uncal_files.extend(glob.glob(jwst_uncal_file))
        else:
            all_jwst_uncal_files.append(jwst_uncal_file)
    
    all_jwst_uncal_files = list(set(sorted(all_jwst_uncal_files)))
    
    for jwst_uncal_file in all_jwst_uncal_files:
        
        jwst_data_base_name = os.path.basename(jwst_uncal_file)
        regex_match = re.match(r'^(jw[0-9]+_[0-9]+_[0-9]+)_([a-z0-9]+)_uncal.fits$', jwst_data_base_name)
        if regex_match:
            jwst_data_base_str = regex_match.group(1)
            jwst_data_detector_str = regex_match.group(2)
            payload = crds.client.get_aui_best_references(pipeline_context, [jwst_data_base_str+'.'+jwst_data_detector_str])
            fcache = crds.api.FileCacher(pipeline_context, ignore_cache=False, raise_exceptions=False)
            fcache.get_local_files([pipeline_context])
        
        for key in payload.keys():
            status = payload[key][0]
            if status is not False:
                bestrefs = payload[key][1]
                fcache.get_local_files(bestrefs)
                #for bestref in bestrefs:
                #    crds.client.get_flex_uri(bestref)
        
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
        
        
        if regex_match:
            try_rate_name = jwst_data_base_str + '_' + jwst_data_detector_str + '_rate.fits'
            try_rate_file = os.path.join(os.path.dirname(os.path.dirname(jwst_uncal_file)), 
                'calibrated1_rates', try_rate_name)
            if os.path.isfile(try_rate_file):
                logger.info('Image2Pipeline._precache_references: {!r}'.format(jwst_uncal_file))
                pipeline_object = calwebb_image2.Image2Pipeline()
                pipeline_object._precache_references(jwst_uncal_file)
                




# Main
if __name__ == '__main__':
    
    main()



