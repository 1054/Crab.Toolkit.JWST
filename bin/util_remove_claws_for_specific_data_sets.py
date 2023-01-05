#!/usr/bin/env python 
# 

"""
Remove claws for specific data sets. 

We need some region files under the `CODE_DATADIR` directory with the 
same name as the data set file, e.g., 'jw01837022001_02201_00002_nrca1_rate.reg'. 

"""

import os, sys, re, copy, shutil, glob, datetime
import click
import numpy as np
import tqdm
from astropy.io import fits
from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import sigma_clipped_stats
#from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
#from scipy.optimize import curve_fit
#from tqdm import tqdm
#import lmfit
#from scipy.optimize import minimize, LinearConstraint
#from lmfit import Minimizer, Parameters, report_fit

# jwst
from jwst.datamodels import ImageModel, FlatModel, dqflags
#from jwst.datamodels.dqflags import pixel as dqflags_pixel
#from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
#from stdatamodels import util
#import crds

# code name and version
CODE_NAME = 'util_remove_claws_for_specific_data_sets.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230105'
CODE_HOMEPAGE = 'https://github.com/1054/Crab.Toolkit.JWST'
CODE_DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'util_remove_claws_for_specific_data_sets_region_files')

# logging
import logging
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--template-dir', type=click.Path(exists=True), default=CODE_DATADIR, help='Where to find the custom claws region files.')
def main(
        input_file, 
        output_file, 
        template_dir, 
    ):
    """
    Input rate image. 
    """
    
    # find region file(s)
    dataset_name = os.path.splitext(os.path.basename(input_file))[0]
    region_files = glob.glob(os.path.join(template_dir, dataset_name + '*.reg'))
    if len(region_files) == 0:
        region_files = glob.glob(os.path.join(template_dir, '*', dataset_name + '*.reg'))
    if len(region_files) == 0:
        logger.warning('No region file is found associating with the data set {!r} in the directory {!r}'.format(
            dataset_name, template_dir, 
        ))
        return
    
    # check input image
    with ImageModel(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Removed claws'):
                    logger.info('{!r} already had claws removed. Skipping!'.format(input_file))
                    return
    
        # get image
        image = model.data.copy()
        imerr = model.err.copy()
        
        # check instrument_name
        instrument_name = model.meta.instrument.name
        if instrument_name.upper() != 'NIRCAM':
            logger.warning('The input is not NIRCam data, cannot remove claws.')
            return
    
    
    # 
    
    
    
    # output file
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_removing_claws.fits'
        if not os.path.isfile(output_orig):
            shutil.copy2(output_file, output_orig) # do not overwrite
    else:
        if os.path.isfile(output_file):
            shutil.move(output_file, output_file+'.backup')
        elif os.path.dirname(output_file) != '' and not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        shutil.copy2(input_file, output_file) # copy input to output then update the model.data
    
    out_timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with ImageModel(output_file) as out_model:
        out_model.data = subtracted_image
        # add history entry following CEERS
        stepdescription = f'Removed claws with {CODE_NAME} ({out_timestamp})'
        software_dict = {'name':CODE_NAME,
                         'author':CODE_AUTHOR,
                         'version':CODE_VERSION,
                         'homepage':CODE_HOMEPAGE}
        substr = util.create_history_entry(stepdescription, software=software_dict)
        out_model.history.append(substr)
        print('out_model.history', out_model.history)
        out_model.save(output_file)
        logger.info('Saved claws-removed data into {!r}'.format(output_file))
    
    
    



if __name__ == '__main__':
    main()



