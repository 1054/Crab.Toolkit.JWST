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
@click.argument('jwst_data_dir_name', type=str)
def main(jwst_data_dir_name):

    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Check input
    jwst_data_dir_path = jwst_data_dir_name
    if jwst_data_dir_name.find(os.sep) >= 0:
        jwst_data_dir_name = os.path.basename(jwst_data_dir_name)
    
    regex_pattern = r'(jw[0-9]{5}[0-9]{3}[0-9]{3}_[0-9]{5}_[0-9]{5})_([a-zA-Z0-9]+)' # e.g., {jw01345004001_04201_00003}_{nrcb4}
    regex_match = re.match(regex_pattern, jwst_data_dir_name)
    if not regex_match:
        raise Exception('Error! The input jwst_data_dir_name "{}" could not match the regular expression "{}"'.format(jwst_data_dir_name, regex_pattern))
    
    jwst_data_base_str = regex_match.group(1)
    jwst_data_detector_str = regex_match.group(2)
    
    # Print CRDS pipeline version
    logger.info('CRDS pipeline version: {}'.format(crds.__version__))
    
    # 
    pipeline_context = crds.client.get_default_context('jwst')
    logger.info('pipeline_context: {}'.format(pipeline_context))
    
    #crds.get_cached_mapping(pipeline_context) # need local file to exist
    
    #crds.rmap.load_mapping(pipeline_context)
    
    #crds.getreferences(parameters={'INSTRUME':'NIRCAM', }, 
    #                   reftypes=['DARK'], 
    #                   context=pipeline_context,
    #                   ignore_cache=False,
    #                   observatory='jwst')
    
    payload = crds.client.get_aui_best_references(pipeline_context, [jwst_data_base_str+'.'+jwst_data_detector_str])
    # {'JW01837014001_02101_00001_MIRIMAGE': [True, [
    #     'jwst_miri_abvegaoffset_0001.asdf', 'jwst_miri_apcorr_0008.fits', 
    #     'jwst_miri_area_0004.fits', 'jwst_miri_dark_0082.fits', 
    #     'jwst_miri_distortion_0047.asdf', 'jwst_miri_drizpars_0001.fits', 
    #     'jwst_miri_filteroffset_0006.asdf', 'jwst_miri_flat_0785.fits', 
    #     'jwst_miri_gain_0008.fits', 'jwst_miri_ipc_0009.fits', 
    #     'jwst_miri_linearity_0032.fits', 'jwst_miri_mask_0030.fits', 
    #     'jwst_miri_pars-darkpipeline_0001.asdf', 'jwst_miri_pars-detector1pipeline_0001.asdf', 
    #     'jwst_miri_pars-outlierdetectionstep_0037.asdf', 'jwst_miri_pars-sourcecatalogstep_0018.asdf', 
    #     'jwst_miri_pars-tweakregstep_0020.asdf', 'jwst_miri_photom_0074.fits', 
    #     'jwst_miri_readnoise_0085.fits', 'jwst_miri_reset_0070.fits', 
    #     'jwst_miri_rscd_0017.fits', 'jwst_miri_saturation_0027.fits'
    # ]]}
    
    fcache = crds.api.FileCacher(pipeline_context, ignore_cache=False, raise_exceptions=False)
    
    fcache.get_local_files([pipeline_context])
    
    for key in payload.keys():
        status = payload[key][0]
        if status is not False:
            bestrefs = payload[key][1]
            fcache.get_local_files(bestrefs)
            #for bestref in bestrefs:
            #    crds.client.get_flex_uri(bestref)
    
    
    # 
    fitsfile = f'{jwst_data_dir_path}/uncals/{jwst_data_dir_name}_uncal.fits'
    if os.path.exists(fitsfile):
        pipeline_object = calwebb_detector1.Detector1Pipeline()
        pipeline_object._precache_references(fitsfile)




# Main
if __name__ == '__main__':
    
    main()



