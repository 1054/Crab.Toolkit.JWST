#!/usr/bin/env python
# 
"""
Use a source emission seed image as a mask and set data quality flag in the rate.fits file.

Usage: 
    ./util_mask_rate_data_with_seed_image.py \
        jw01837022001_02101_00002_mirimage/calibrated1_rates/jw01837022001_02101_00002_mirimage_rate.fits \
        calibrated3_mosaics/jw01837_obs022_MIRI_F770W/jw01837_obs022_MIRI_F770W_destripping_i2d_remstriping_rate_MIRI_F770W_galaxy_seed_image.fits \
        jw01837022001_02101_00002_mirimage/calibrated1_rates/jw01837022001_02101_00002_mirimage_rate_masked_source_emission.fits

By Daizhong Liu @MPE. 

Last update: 2022-09-10.
"""
import os, sys, re, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.modeling import models as apy_models
from astropy.modeling import fitting as apy_fitting
from astropy.convolution import convolve as apy_convolve
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from photutils.background import Background2D
from reproject import reproject_interp
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300

if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

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

from jwst.datamodels import ImageModel, FlatModel


def mask_rate_data_with_seed_image(
        rate_image_file, 
        seed_image_file, 
        out_image_file, 
        check_rate_image = False, # whether to check file name ends with "_rate.fits"
        miri_flat_file = None, # for MIRI only, will read DQ from the flat file.
        overwrite = True, 
    ):

    # 
    print('Running mask_rate_data_with_seed_image')
    print('rate_image_file: {!r}'.format(rate_image_file))
    print('seed_image_file: {!r}'.format(seed_image_file))
    print('out_image_file: {!r}'.format(out_image_file))
    
    # check input rate image
    if check_rate_image:
        if not rate_image_file.endswith('_rate.fits'):
            raise Exception('Error! The input FITS image file name should ends with \"_rate.fits\"!')
        data_name = re.sub(r'^(.*)_rate\.fits$', r'\1', rate_image_file)
    else:
        data_name = os.path.splitext(rate_image_file)[0]
        
    # check out_image_file
    if os.path.isfile(out_image_file) and not overwrite:
        print('Found existing output image {!r} and overwrite is set to False. Will not do anything.'.format(out_image_file))
        return
    
    # read input rate image header
    rate_image, rate_image_header = fits.getdata(rate_image_file, extname='SCI', header=True)
    
    # read seed image
    source_image, source_image_header = fits.getdata(seed_image_file, header=True)
    
    # get source emission mask
    # reproject if needed
    if source_image.shape == rate_image.shape:
        source_mask = (source_image>1e-6)
    else:
        print('Reprojecting seed image to rate image')
        source_image_reproj = reproject_interp(
            (source_image, source_image_header), 
            rate_image_header,
            return_footprint=False,
        )
        source_mask = (source_image_reproj>1e-6)
    
    # prepare mask_source_emission_dq
    mask_source_emission_dq = interpret_bit_flags(
                                    'DO_NOT_USE+SATURATED+NO_FLAT_FIELD+UNRELIABLE_DARK+UNRELIABLE_FLAT+OTHER_BAD_PIXEL',
                                    flag_name_map=dqflags_pixel
                                )
    
    # read input rate image as ImageModel
    with datamodels.open(rate_image_file) as image_model:
        
        print(type(image_model.dq))
        
        print('image_model.dq[884,145]', image_model.dq[884,145]) # Lyot CoronGraph
        print('image_model.dq[688,666]', image_model.dq[688,666]) # source emission
        
        # set DQ to the pixels with source emission mask
        image_model.dq[source_mask] = np.bitwise_or(
            image_model.dq[source_mask], 
            mask_source_emission_dq
        )
        
        image_model.data[source_mask] = 0.0
        
        print('image_model.dq[884,145]', image_model.dq[884,145]) # Lyot CoronGraph
        print('image_model.dq[688,666]', image_model.dq[688,666]) # source emission
        
        # write out_image_file
        if os.path.isfile(out_image_file):
            shutil.move(out_image_file, out_image_file+'.backup')
        elif not os.path.isdir(os.path.dirname(out_image_file)):
            os.makedirs(os.path.dirname(out_image_file))
        image_model.write(out_image_file)
        
        print('Output to {!r}'.format(out_image_file))




@click.command()
@click.argument('rate_image_file', type=click.Path(exists=True))
@click.argument('seed_image_file', type=click.Path(exists=True))
@click.argument('out_image_file', type=click.Path(exists=False))
@click.option('--miri-flat-file', type=click.Path(exists=True), default=None)
@click.option('--check-rate-image/--no-check-rate-image', is_flag=True, default=False)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=True)
def main(rate_image_file, 
         seed_image_file, 
         out_image_file, 
         miri_flat_file, 
         check_rate_image, 
         overwrite):
    
    mask_rate_data_with_seed_image(
        rate_image_file, 
        seed_image_file, 
        out_image_file, 
        miri_flat_file = miri_flat_file, 
        check_rate_image = check_rate_image, 
        overwrite = overwrite, 
    )



if __name__ == '__main__':
    main()

    
