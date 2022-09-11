#!/usr/bin/env python
# 
"""
Merge a list of rate data which have source emission masked out.

Usage: 
    ./util_merge_source_emission_masked_rate_data.py \
        jw01837014001_02101_00002_mirimage/calibrated1_rates/jw01837014001_02101_00002_mirimage_rate_masked_source_emission.fits \
        jw01837023001_02101_00002_mirimage/calibrated1_rates/jw01837023001_02101_00002_mirimage_rate_masked_source_emission.fits \
        jw01837022001_02101_00002_mirimage/calibrated1_rates/merged_other_visit_rates_with_source_emission_mask.fits

By Daizhong Liu @MPE. 

Last update: 2022-09-10.
"""
import os, sys, re, copy, shutil
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

# logging
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def merge_source_emission_masked_rate_data(
        rate_image_files, 
        out_image_file, 
        do_sigma_clip = False, 
        sigma = 3.0, 
        overwrite = True, 
    ):

    # 
    logger.info('Running merge_source_emission_masked_rate_data')
    logger.info('rate_image_files: {!r} ({})'.format(rate_image_files, len(rate_image_files)))
    logger.info('out_image_file: {!r}'.format(out_image_file))
        
    # check out_image_file
    if os.path.isfile(out_image_file) and not overwrite:
        logger.info('Found existing output image {!r} and overwrite is set to False. Will not do anything.'.format(out_image_file))
        return
    
    # prepare mask_source_emission_dq
    # same as in "util_mask_rate_data_with_seed_image.py"
    mask_source_emission_dq = interpret_bit_flags(
                                    'DO_NOT_USE+SATURATED+NO_FLAT_FIELD+UNRELIABLE_DARK+UNRELIABLE_FLAT+OTHER_BAD_PIXEL',
                                    flag_name_map=dqflags_pixel
                                )
    
    # read input rate image header
    rate_images_sci = None
    rate_images_err = None
    rate_images_count = None
    rate_images_dq = None
    for i, rate_image_file in enumerate(rate_image_files):
        rate_image_sci, rate_image_header = fits.getdata(rate_image_file, extname='SCI', header=True)
        rate_image_err = fits.getdata(rate_image_file, extname='ERR', header=False)
        rate_image_dq = fits.getdata(rate_image_file, extname='DQ', header=False)
        if rate_images_sci is None:
            rate_images_sci = np.full(rate_image_sci.shape, fill_value=0.0)
            rate_images_err = np.full(rate_image_sci.shape, fill_value=0.0)
            rate_images_count = np.full(rate_image_sci.shape, fill_value=0.0)
            rate_images_dq = copy.copy(rate_image_dq)
        # 
        #mask_good_pixels = bitfield_to_boolean_mask(
        #    rate_image_dq, 
        #    interpret_bit_flags('GOOD',
        #        flag_name_map=dqflags_pixel
        #    ),
        #    good_mask_value=True,
        #    dtype=bool
        #)
        mask_good_pixels = (rate_image_dq == 0)
        # 
        if do_sigma_clip:
            mask = np.logical_and(mask_good_pixels, rate_image_sci>0.0)
            sigma_clip_median = np.median(rate_image_sci[mask].ravel())
            sigma_clip_stddev = np.std(rate_image_sci[mask].ravel() - sigma_clip_median)
            sigma_clip_thresh = sigma_clip_median + sigma * sigma_clip_stddev
            logger.info('sigma clipping: median {}, stddev {}, median+sigma*stddev {}'.format(
                sigma_clip_median, sigma_clip_stddev, sigma_clip_thresh
                ))
            sigma_clip_mask = np.logical_and(mask, rate_image_sci>sigma_clip_thresh)
            rate_image_sci[sigma_clip_mask] = 0.0
            rate_image_dq[sigma_clip_mask] = mask_source_emission_dq
        # 
        mask = np.logical_and(np.isfinite(rate_image_sci), rate_image_sci>0.0)
        # 
        rate_images_sci[mask] += rate_image_sci[mask]
        rate_images_err[mask] += rate_image_err[mask]
        rate_images_count[mask] += 1.0
        rate_images_dq = np.bitwise_or(rate_images_dq, rate_image_dq)
    
    # compute output rate image
    out_image_sci = rate_images_sci * 0.0
    out_image_err = rate_images_err * 0.0
    mask = (rate_images_count > 0)
    out_image_sci[mask] = rate_images_sci[mask]/rate_images_count[mask]
    out_image_err[mask] = rate_images_err[mask]/rate_images_count[mask] #<TODO># error propagation?
    
    # restore dq to good pixel if the masked source emission dq has valid pixel values
    out_image_dq = copy.copy(rate_images_dq)
    mask = np.logical_and(out_image_dq == mask_source_emission_dq, out_image_sci > 0.0)
    out_image_dq[mask] = 0 # good
    
    # read the first input rate image as a ImageModel template and save output rate data
    with datamodels.open(rate_image_file) as image_model:
        
        image_model.data = out_image_sci
        image_model.err[out_image_err>0] = out_image_err[out_image_err>0]
        image_model.dq = out_image_dq
        
        # write out_image_file
        if os.path.isfile(out_image_file):
            shutil.move(out_image_file, out_image_file+'.backup')
        elif not os.path.isdir(os.path.dirname(out_image_file)):
            os.makedirs(os.path.dirname(out_image_file))
        image_model.write(out_image_file)
        
        logger.info('Output to {!r}'.format(out_image_file))




@click.command()
@click.argument('rate_image_files', nargs=-1, type=click.Path(exists=True))
@click.argument('out_image_file', nargs=1, type=click.Path(exists=False))
@click.option('--do-sigma-clip/--no-sigma-clip', is_flag=True, default=False)
@click.option('--sigma', type=float, default=3.0, help='Clipping sigma. Only clipping rate image pixels brighter than this sigma.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=True)
def main(rate_image_files, 
         out_image_file, 
         do_sigma_clip, 
         sigma, 
         overwrite):
    
    merge_source_emission_masked_rate_data(
        rate_image_files, 
        out_image_file, 
        do_sigma_clip = do_sigma_clip, 
        sigma = sigma, 
        overwrite = overwrite, 
    )



if __name__ == '__main__':
    main()

    
