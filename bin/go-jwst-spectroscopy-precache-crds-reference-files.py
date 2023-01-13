#!/usr/bin/env python
#
"""
Download CRDS cache files.

By Daizhong Liu.

"""

import os, sys, re, datetime, glob, shutil
assert os.environ["CRDS_PATH"] != ''
assert os.environ["CRDS_SERVER_URL"] != ''
#os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
#os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Import CRDS
import crds
from stpipe import crds_client

# Import JWST pipeline
from jwst import datamodels
from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_spec2
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_spec3

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
        
        # jwst_data_base_name = os.path.basename(jwst_uncal_file)
        # regex_match = re.match(r'^(jw[0-9]+_[0-9]+_[0-9]+)_([a-z0-9]+)_uncal.fits$', jwst_data_base_name)
        # if regex_match:
        #     jwst_data_base_str = regex_match.group(1)
        #     jwst_data_detector_str = regex_match.group(2)
        #     payload = crds.client.get_aui_best_references(pipeline_context, [jwst_data_base_str+'.'+jwst_data_detector_str])
        #     fcache = crds.api.FileCacher(pipeline_context, ignore_cache=False, raise_exceptions=False)
        #     fcache.get_local_files([pipeline_context])
        
        # for key in payload.keys():
        #     status = payload[key][0]
        #     if status is not False:
        #         bestrefs = payload[key][1]
        #         fcache.get_local_files(bestrefs)
        #         #for bestref in bestrefs:
        #         #    crds.client.get_flex_uri(bestref)
        
        # # with datamodels.open(jwst_uncal_file) as model:
        # #     params = {
        # #         'INSTRUME': model.meta.instrument.name, 
        # #         'DATE': model.meta.date, 
        # #         'TIME': model.meta.observation.time,
        # #     }
        # #     crds.getreferences(
        # #         parameters = parameters, 
        # #         reftypes = ['DARK'], 
        # #         context = pipeline_context,
        # #         ignore_cache = False,
        # #         observatory = 'jwst',
        # #     )
        
        
        with datamodels.open(jwst_uncal_file) as model:
            exp_type = model.meta.exposure.type
        
        logger.info('Detector1Pipeline._precache_references: {!r}'.format(jwst_uncal_file))
        pipeline_object = calwebb_detector1.Detector1Pipeline()
        pipeline_object._precache_references(jwst_uncal_file)
        
        
        logger.info('Spec2Pipeline._precache_references: {!r}'.format(jwst_uncal_file))
        if exp_type in ['NRS_IMAGE', 'NRS_WATA', 'NRS_MSATA', 'NRS_TACONFIRM', 'NRS_CONFIRM', 'NRS_FOCUS', 'NRS_MIMF']:
            pipeline_object = calwebb_image2.Image2Pipeline()
            pipeline_object._precache_references(jwst_uncal_file)
            
            # No need stage 3
            # see https://jwst-pipeline.readthedocs.io/_/downloads/en/stable/pdf/
            # Table 4
            
        else:
            pipeline_object = calwebb_spec2.Spec2Pipeline()
            pipeline_object._precache_references(jwst_uncal_file)
            # 
            #shutil.copy2(jwst_uncal_file, jwst_uncal_file+'.tmp')
            #with datamodels.open(jwst_uncal_file+'.tmp') as model:
            #    if model.meta.exposure.type in ['NRS_MSATA', 'NRS_TACONFIRM']:
            #        model.meta.exposure.type = 'NRS_MSASPEC'
            #    #pipeline_object._precache_references(model) # this will also fail
            #    ovr_refs = {reftype: pipeline_object.get_ref_override(reftype) 
            #                    for reftype in pipeline_object.reference_file_types 
            #                    if pipeline_object.get_ref_override(reftype) is not None}
            #    fetch_types = sorted(set(pipeline_object.reference_file_types) - set(ovr_refs.keys()))
            #    for key in ['sflat', 'area']:
            #        if key in fetch_types:
            #            fetch_types.remove(key)
            #    logger.info("Prefetching reference files for dataset: " + repr(model.meta.filename) +
            #                " reftypes = " + repr(fetch_types)) # following "stpipe/pipeline.py"
            #    crds_refs = crds_client.get_multiple_reference_paths(
            #        model.get_crds_parameters(), 
            #        fetch_types, 
            #        model.crds_observatory
            #    )
            #os.remove(jwst_uncal_file+'.tmp')
            
            
            logger.info('Spec3Pipeline._precache_references: {!r}'.format(jwst_uncal_file))
            pipeline_object = calwebb_spec3.Spec3Pipeline()
            pipeline_object._precache_references(jwst_uncal_file) # this will fail
            # 
            #shutil.copy2(jwst_uncal_file, jwst_uncal_file+'.tmp')
            # with datamodels.open(jwst_uncal_file+'.tmp') as model:
            #     if model.meta.exposure.type in ['NRS_MSATA', 'NRS_TACONFIRM']:
            #         model.meta.exposure.type = 'NRS_MSASPEC'
            #     #pipeline_object._precache_references(model) # this will also fail
            #     ovr_refs = {reftype: pipeline_object.get_ref_override(reftype) 
            #                     for reftype in pipeline_object.reference_file_types 
            #                     if pipeline_object.get_ref_override(reftype) is not None}
            #     fetch_types = sorted(set(pipeline_object.reference_file_types) - set(ovr_refs.keys()))
            #     for key in ['area']:
            #         if key in fetch_types:
            #             fetch_types.remove(key)
            #     logger.info("Prefetching reference files for dataset: " + repr(model.meta.filename) +
            #                 " reftypes = " + repr(fetch_types)) # following "stpipe/pipeline.py"
            #     crds_refs = crds_client.get_multiple_reference_paths(
            #         model.get_crds_parameters(), 
            #         fetch_types, 
            #         model.crds_observatory
            #     )
            #     ref_path_map = dict(list(crds_refs.items()) + list(ovr_refs.items()))
            #     for (reftype, refpath) in sorted(ref_path_map.items()):
            #         how = "Override" if reftype in ovr_refs else "Prefetch"
            #         logger.info(f"{how} for {reftype.upper()} reference file is '{refpath}'.")
            #         crds_client.check_reference_open(refpath)
            # os.remove(jwst_uncal_file+'.tmp')
        
        
        
        
        # create a timestamp file
        if os.path.dirname(jwst_uncal_file) == 'uncals':
            timestamp_file = os.path.dirname(os.path.dirname(jwst_uncal_file))+os.sep+'crds_cached'
            timestamp_str = datetime.datetime.now().strftime('%Y-%m-%d %Hh%Mm%Ss')
            with open(timestamp_file, 'w') as fp:
                fp.write(timestamp_str+'\n')
        




# Main
if __name__ == '__main__':
    
    main()



