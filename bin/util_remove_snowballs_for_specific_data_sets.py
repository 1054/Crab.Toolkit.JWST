#!/usr/bin/env python 
# 

"""
Remove snowballs for specific data sets. 

We need some region files under the `CODE_DATADIR` directory with the 
same name as the data set file, e.g., 'jw01837022001_02201_00002_nrca1_rate.reg'. 

20230604: It can also read region files and directly add those regions into DQ mask. 
The regions need to be pixel regions.

"""

import os, sys, re, copy, shutil, glob, datetime
import click
import numpy as np
#import tqdm
from astropy.io import fits
from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS, FITSFixedWarning
#from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
#from scipy.optimize import curve_fit
#from tqdm import tqdm
#import lmfit
#from scipy.optimize import minimize, LinearConstraint
#from lmfit import Minimizer, Parameters, report_fit
from regions import Regions, PixelRegion, SkyRegion
from reproject import reproject_interp

# jwst
from jwst import datamodels

from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)

from jwst.datamodels.dqflags import pixel as dqflags_pixel
#{'GOOD': 0, 
#'DO_NOT_USE': 1, 'SATURATED': 2, 'JUMP_DET': 4, 'DROPOUT': 8, 
#'OUTLIER': 16, 'PERSISTENCE': 32, 'AD_FLOOR': 64, 'RESERVED_4': 128, 
#'UNRELIABLE_ERROR': 256, 'NON_SCIENCE': 512, 'DEAD': 1024, 'HOT': 2048, 
#'WARM': 4096, 'LOW_QE': 8192, 'RC': 16384, 'TELEGRAPH': 32768, 
#'NONLINEAR': 65536, 'BAD_REF_PIXEL': 131072, 'NO_FLAT_FIELD': 262144, 'NO_GAIN_VALUE': 524288, 
#'NO_LIN_CORR': 1048576, 'NO_SAT_CHECK': 2097152, 'UNRELIABLE_BIAS': 4194304, 'UNRELIABLE_DARK': 8388608, 
#'UNRELIABLE_SLOPE': 16777216, 'UNRELIABLE_FLAT': 33554432, 'OPEN': 67108864, 'ADJ_OPEN': 134217728, 
#'UNRELIABLE_RESET': 268435456, 'MSA_FAILED_OPEN': 536870912, 'OTHER_BAD_PIXEL': 1073741824, 'REFERENCE_PIXEL': 2147483648
#}

from jwst.datamodels import ImageModel, FlatModel, dqflags
#from jwst.datamodels.dqflags import pixel as dqflags_pixel
#from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
from stdatamodels import util
#import crds

import warnings
warnings.filterwarnings('ignore', category=FITSFixedWarning)

# code name and version
CODE_NAME = 'util_remove_snowballs_for_specific_data_sets.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230105'
CODE_HOMEPAGE = 'https://github.com/1054/Crab.Toolkit.JWST'
CODE_DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                            'util_remove_snowballs_for_specific_data_sets_region_and_mask_files')

# logging
import logging
logging.basicConfig(level='INFO')
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--template-dir', type=click.Path(exists=True), default=CODE_DATADIR, help='Where to find the custom snowballs region files.')
@click.option('--fill-value', type=float, default=None, help='Fill the data pixel value in the mask area. Default is None. Their DQ will be updated anyway.')
def main(
        input_file, 
        output_file, 
        template_dir, 
        fill_value, 
    ):
    """
    Input rate image. 
    """
    
    # find region and mask file(s)
    # TODO: region files
    dataset_name = os.path.splitext(os.path.basename(input_file))[0]
    dataset_name = re.sub(r'^(.*)(_rate|_cal)$', r'\1', dataset_name)
    region_files = glob.glob(os.path.join(template_dir, dataset_name + '*.reg'))
    mask_files = glob.glob(os.path.join(template_dir, dataset_name + '*.fits.gz')) + \
                 glob.glob(os.path.join(template_dir, dataset_name + '*.fits'))
    if len(region_files) == 0: # search subsubdirectories if needed
        region_files = glob.glob(os.path.join(template_dir, '*', dataset_name + '*.reg'))
    if len(mask_files) == 0: # search subsubdirectories if needed
        mask_files = glob.glob(os.path.join(template_dir, '*', dataset_name + '*.fits.gz')) + \
                     glob.glob(os.path.join(template_dir, '*', dataset_name + '*.fits'))
    if len(region_files) == 0 and len(mask_files) == 0:
        logger.warning('No region or mask file is found associating with the data set {!r} in the directory {!r}. Skipping.'.format(
            dataset_name, template_dir, 
        ))
        return
    elif len(mask_files) > 0:
        logger.info('Found mask file(s): {}'.format(mask_files))
        if len(region_files) > 0:
            region_files = []
    elif len(region_files) > 0:
        logger.info('Found region file(s): {}'.format(region_files))
    
    
    # check input image
    logger.info('Reading input data: {!r}'.format(input_file))
    with datamodels.open(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith(f'Removed snowballs with {CODE_NAME}'):
                    logger.info(f'{input_file!r} already had snowballs removed with {CODE_NAME}. Skipping!')
                    return
        
        # get image
        rate_image = model.data.copy()
        
        # check instrument_name
        instrument_name = model.meta.instrument.name
        if instrument_name.upper() != 'NIRCAM':
            logger.warning('The input is not NIRCam data, cannot remove snowballs.')
            return
    
    # read fits header (I don't know how to get the header from model.meta)
    rate_image_header = fits.getheader(input_file, 'SCI')
    
    
    # for loop
    seed_masks = None
    for seed_image_file in mask_files:
        
        # read seed image
        seed_image, seed_image_header = fits.getdata(seed_image_file, header=True)
        
        # get seed emission mask
        # reproject if needed
        if seed_image.shape == rate_image.shape:
            seed_mask = (seed_image>1e-6) # 1e-6 is just a threshold value, maybe use 0 is also okay.
        else:
            logger.info('Reprojecting seed image to rate image')
            seed_image_reproj = reproject_interp(
                (seed_image, seed_image_header), 
                rate_image_header,
                return_footprint=False,
            )
            seed_mask = (seed_image_reproj>1e-6) # 1e-6 is just a threshold value, maybe use 0 is also okay.
        
        if seed_masks is None:
            seed_masks = seed_mask
        else:
            seed_masks = np.bitwise_or(seed_masks, seed_mask)
    
    
    # 
    wcs = None
    for region_file in region_files:
        
        # read pixel region
        region_list = Regions.read(region_file, format='ds9')
        for iregion, region_object in enumerate(region_list):
            if not isinstance(region_object, PixelRegion):
                ith = '{}th'.format(iregion) if iregion >= 3 else (['1st', '2nd', '3rd'][iregion])
                logger.info('Warning! The {} region in {!r} is not a pixel region! The WCS may be inaccurate due to astrometry.'.format(ith, region_file))
                if wcs is None:
                    wcs = WCS(rate_image_header, naxis=2)
                pixel_region = region_object.to_pixel(wcs)
            else:
                pixel_region = region_object
            pixel_mask = pixel_region.to_mask().to_image(rate_image.shape).astype(bool)
            
            if seed_masks is None:
                seed_masks = pixel_mask
            else:
                seed_masks = np.bitwise_or(seed_masks, pixel_mask)
    
    
    # prepare mask_seed_emission_dq -- a specific combination of flags for the masked bad pixels
    mask_seed_emission_dq = interpret_bit_flags(
        'JUMP_DET',
        flag_name_map=dqflags_pixel
    )
    
    
    # output file
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_removing_snowballs.fits'
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
        
        # set DQ to the pixels with seed emission mask
        out_model.dq[seed_masks] = np.bitwise_or(
            out_model.dq[seed_masks], 
            mask_seed_emission_dq
        )
        
        if fill_value is not None:
            out_model.data[seed_masks] = fill_value
        
        # add history entry following CEERS
        stepdescription = f'Removed snowballs with {CODE_NAME} ({out_timestamp})'
        software_dict = {'name':CODE_NAME,
                         'author':CODE_AUTHOR,
                         'version':CODE_VERSION,
                         'homepage':CODE_HOMEPAGE}
        substr = util.create_history_entry(stepdescription, software=software_dict)
        out_model.history.append(substr)
        #print('out_model.history', out_model.history)
        out_model.save(output_file)
        logger.info('Saved snowballs-removed data into {!r}'.format(output_file))
    
    
    



if __name__ == '__main__':
    main()



