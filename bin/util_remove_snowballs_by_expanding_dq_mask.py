#!/usr/bin/env python 
# 

"""
Remove snowballs by expanding large cosmic ray DQ masks.

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
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
from tqdm import tqdm

# jwst
from jwst import datamodels
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.datamodels.dqflags import pixel as dqflags_pixel
from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
from stdatamodels import util
import crds

# code name and version
CODE_NAME = 'util_remove_snowballs_by_expanding_dq_mask.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221109'
CODE_HOMEPAGE = ''

# logging
import logging
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--mask-image', type=click.Path(exists=True), default=None, help='A mask image containing pixels which will not be analyzed, e.g., due to bright emission.')
@click.option('--mask-threshold', type=float, default=None, help='Pixels with values above this threshold in the mask image will be not be analyzed.')
@click.option('--apply-flat', is_flag=True, default=None, help='Apply CRDS flat correction before destripping?')
def main(
        input_file, 
        output_file, 
        mask_image, 
        mask_threshold, 
        apply_flat, 
    ):
    """
    Input rate image. 
    """
    with datamodels.open(input_file) as model:
        
        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Removed snowballs'):
                    logger.info('{!r} already had snowballs removed. Skipping!'.format(input_file))
                    return
        
        # Get DQ mask
        dqmask = bitfield_to_boolean_mask(
                model.dq,
                interpret_bit_flags('~DO_NOT_USE+NON_SCIENCE', # SkyMatchStep().dqbits is '~DO_NOT_USE+NON_SCIENCE' in default
                    flag_name_map=dqflags_pixel), 
                good_mask_value=1,
                dtype=np.uint8
            ) * np.isfinite(model.data) # following "jwst/skymatch/skymatch_step.py", this is valid data
        
        jumpmask = bitfield_to_boolean_mask(
                model.dq,
                interpret_bit_flags('~JUMP_DET', 
                    flag_name_map=dqflags_pixel), 
                good_mask_value=0,
                dtype=np.uint8
            ) # jump_det pixels have a pixel value of 1
        
        satumask = bitfield_to_boolean_mask(
                model.dq,
                interpret_bit_flags('~SATURATED', 
                    flag_name_map=dqflags_pixel), 
                good_mask_value=0,
                dtype=np.uint8
            ) # saturated pixels have a pixel value of 1
        
        
        raise NotImplementedError()
        
        
        from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
        jumpmask_medfil = median_filter(jumpmask, size=8) # filters out small mask areas
        #jumpmask_medfil = gaussian_filter(jumpmask, sigma=1) # filters out small mask areas
        #jumpmask_gausfil = gaussian_filter(jumpmask_medfil, sigma=1) # filters out small areas but shrinks large ones
        jumpmask_gausfil = jumpmask_medfil
        #jumpmask_dilated = binary_dilation(jumpmask_gausfil, iterations=30) # expand areas
        jumpmask_dilated = jumpmask_gausfil
        from astropy.convolution import Gaussian2DKernel, convolve
        jumpmask_expanded = convolve(jumpmask_dilated, kernel=Gaussian2DKernel(x_stddev=0.5))
        jumpmask_expanded[jumpmask_expanded<0.1] = 0
        jumpmask_expanded[jumpmask_expanded>=0.1] = 1
        jumpmask_expanded = np.logical_or(jumpmask>0, jumpmask_expanded>0).astype(float)

        # # edge detection
        # from scipy.ndimage import sobel
        # sx = sobel(jumpmask_gausfil, axis=0, mode='constant')
        # sy = sobel(jumpmask_gausfil, axis=1, mode='constant')
        # sob = np.hypot(sx, sy)
        # jumpmask_edgedet = sob
        
        
        
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
            validmask = np.isfinite(model.data)
            out_model.data[validmask] = image[validmask] # only copy valid pixel values
            #out_model.data = image
            # check flat
            if apply_flat:
                out_model.meta.cal_step.flat_field = model.meta.cal_step.flat_field
            # add history entry following CEERS
            stepdescription = f'Removed snowballs with {CODE_NAME} ({out_timestamp})'
            software_dict = {'name':CODE_NAME,
                             'author':CODE_AUTHOR,
                             'version':CODE_VERSION,
                             'homepage':CODE_HOMEPAGE}
            substr = util.create_history_entry(stepdescription, software=software_dict)
            out_model.history.append(substr)
            print('out_model.history', out_model.history)
            out_model.save(output_file)
            logger.info('Saved destripped data into {!r}'.format(output_file))
        
        
        output_bkg = re.sub(r'\.fits$', r'', output_file) + '_bkg.fits'
        shutil.copy2(output_file, output_bkg)
        with ImageModel(output_bkg) as out_bkg_model:
            out_bkg_model.data = bkg
            out_bkg_model.history.append(substr)
            out_bkg_model.save(output_bkg)
            logger.info('Saved destripping background into {!r}'.format(output_bkg))
        
    



if __name__ == '__main__':
    main()



