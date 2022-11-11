#!/usr/bin/env python 
# 

"""
Apply flat field, write history. 

"""

import os, sys, re, copy, shutil, glob, datetime
import click
import numpy as np
import tqdm
from astropy.io import fits
#from astropy.nddata.bitmask import (
#    bitfield_to_boolean_mask,
#    interpret_bit_flags,
#)
#from astropy.stats import sigma_clipped_stats
#from scipy.optimize import curve_fit
#from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
#from tqdm import tqdm

# CRDS
if 'CRDS_PATH' not in os.environ:
    search_paths = [os.path.expanduser('~/jwst_crds_cache')]
    for search_path in search_paths:
        if os.path.isdir(search_path):
            os.environ['CRDS_PATH'] = search_path
            break
if 'CRDS_PATH' not in os.environ:
    raise Exception('Error! CRDS_PATH is not set!')
if 'CRDS_SERVER_URL' not in os.environ:
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
if 'CRDS_CONTEXT' not in os.environ:
    raise Exception('Error! CRDS_CONTEXT is not set!')

# jwst
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.datamodels.dqflags import pixel as dqflags_pixel
from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
from stdatamodels.util import create_history_entry
import crds

# code name and version
CODE_NAME = 'util_apply_flatfield.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221109'
CODE_HOMEPAGE = ''

# logging
import logging
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)


# Apply flat
def apply_flat_to_model(model):
    logger.info('Applying flat to model {}'.format(model))
    try:
        crds_context = os.environ['CRDS_CONTEXT']
    except KeyError:
        crds_context = crds.get_default_context()
    # pull flat from CRDS using the current context
    #<DZLIU># enabling both NIRCAM and MIRI, 
    #<DZLIU># following https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/reference_files.html
    instrument_name = model.meta.instrument.name
    if instrument_name.upper() in ['NIRCAM', 'NIRISS']: #<DZLIU>#
        crds_dict = {'INSTRUME':instrument_name, #<DZLIU>#
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'PUPIL':model.meta.instrument.pupil, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    elif instrument_name.upper() == 'MIRI': #<DZLIU>#
        crds_dict = {'INSTRUME':instrument_name, #<DZLIU>#
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'BAND':model.meta.instrument.band, 
                     'READPATT':model.meta.exposure.readpatt, 
                     'SUBARRAY':model.meta.subarray.name, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    else:
        raise ValueError("Error! The `instrument_name` (model.meta.instrument.name) can only be NIRCAM and MIRI!")
    logger.info('crds_dict: {}'.format(crds_dict)) #<DZLIU>#
    flats = crds.getreferences(crds_dict, reftypes=['flat'], 
                               context=crds_context)
    # if the CRDS loopup fails, should return a CrdsLookupError, but 
    # just in case:
    try:
        flatfile = flats['flat']
    except KeyError:
        logger.error('Flat was not found in CRDS with the parameters: {}'.format(crds_dict))
        raise Exception('Flat was not found in CRDS with the parameters: {}'.format(crds_dict))

    logger.info('Using flat: %s'%(os.path.basename(flatfile)))
    with FlatModel(flatfile) as flat:
        # use the JWST Calibration Pipeline flat fielding Step 
        output_model, applied_flat = do_correction(model, flat)
        output_model.meta.ref_file.flat.name = 'crds://'+os.path.basename(flatfile)
    
    return output_model, applied_flat




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--apply-flat/--no-apply-flat', is_flag=True, default=True, help='Apply CRDS flat correction before destriping?')
def main(
        input_file, 
        output_file, 
        apply_flat, 
    ):
    """
    Input rate image. 
    """
    
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
    
    
    # open input model
    with ImageModel(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Applied flatfield'):
                    logger.info('{!r} already had flatfield applied. Skipping!'.format(input_file))
                    return
        
        # Apply flat (and dq)
        apply_flat_status = None
        if apply_flat:
            output_model, applied_flat = apply_flat_to_model(model)
            apply_flat_status = output_model.meta.cal_step.flat_field
        else:
            logger.warning('The apply_flat is False. Flatfield is not applied!')
            return
    
    # prepare output file
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_applying_flatfield.fits'
        if os.path.isfile(output_orig):
            logger.info('Found backup file {}, skipping backup.'.format(output_orig))
        else:
            logger.info('Backing up input file as {}'.format(output_orig))
            shutil.copy2(output_file, output_orig) # do not overwrite
    else:
        if os.path.isfile(output_file):
            shutil.move(output_file, output_file+'.backup')
        elif os.path.dirname(output_file) != '' and not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        shutil.copy2(input_file, output_file) # copy input to output then update the model.data
    
    
    # add history entry following CEERS
    out_timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    stepdescription = f'Applied flatfield with {CODE_NAME} ({out_timestamp})'
    software_dict = {'name':CODE_NAME,
                     'author':CODE_AUTHOR,
                     'version':CODE_VERSION,
                     'homepage':CODE_HOMEPAGE}
    entry = create_history_entry(stepdescription, software=software_dict)
    # TODO directly calling output_model.save() does not work!
    with ImageModel(output_file) as out_model:
        out_model.update(output_model)
        out_model.data = output_model.data
        out_model.err = output_model.err
        out_model.dq = output_model.dq
        out_model.history.append('')
        out_model.history.append(entry)
        out_model.history.append('')
        logger.info('Adding model history: {}'.format(entry))
        out_model.save(output_file)
        logger.info('Saved flatfield-applied data into {!r}'.format(output_file))
    
    
    # save additional file
    output_flatfile = re.sub(r'\.fits$', r'', output_file) + '_applying_flatfield_flatfile.fits'
    with ImageModel(output_flatfile) as out_model:
        out_model = applied_flat
        out_model.save(output_flatfile)
    logger.info('Saved the applied flat into {!r}'.format(output_flatfile))





if __name__ == '__main__':
    main()



