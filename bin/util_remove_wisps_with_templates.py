#!/usr/bin/env python 
# 

"""
Measure and renormalize strips in JWST images. 

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
from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
from scipy.optimize import curve_fit
#from tqdm import tqdm
#import lmfit
#from scipy.optimize import minimize, LinearConstraint
from lmfit import Minimizer, Parameters, report_fit

# jwst
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.datamodels.dqflags import pixel as dqflags_pixel
from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
from stdatamodels import util
import crds

# code name and version
CODE_NAME = 'util_remove_wisps_with_templates.py'
CODE_AUTHOR = 'Daizhong Liu'
#CODE_VERSION = '20221109'
CODE_VERSION = '20230830' # adding --template-file
CODE_HOMEPAGE = ''

# logging
import logging
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)



# define template fitting function
def template_fitting_function(pars, x, data=None, error=None, weight=None, template=None):
    if weight is None:
        weight = 1.0
    if error is None:
        return np.abs((pars['A'] * template + pars['B'] - data)) * weight
    return np.abs((pars['A'] * template + pars['B'] - data) / error) * weight

# iter_cb
def iter_cb(params, iter, resid, *fcn_args, **fcn_kws):
    logger.info('iteration {}, params {}, residual sum {}'.format(iter, params, np.sum(resid**2)))
    #params.pretty_print()

# define template fitting 
def do_template_fitting(data_image, error_image, template_image):
    
    # first define a broad mask with all valid pixels in both data and template
    # so that we can compute a background difference, i.e. param B.
    valid_mask = np.logical_and.reduce((np.isfinite(template_image), np.isfinite(data_image), (error_image>0.0)))
    init_param_B = np.median(data_image[valid_mask]) - np.median(template_image[valid_mask])
    
    # then define a strict mask in the templates image so that we will focus on matching that part of image
    # there are some annoying small 3x3 areas in the templates image that have a high pixel value and affect later fitting
    # we need to do some median filter then dilation to remove them and get a better mask
    template_max = np.max(template_image[valid_mask])
    template_mean = np.mean(template_image[valid_mask])
    template_stddev = np.std(template_image[valid_mask]-template_mean)
    template_P95 = np.percentile(template_image[valid_mask], 95.) # use top 5% pixels
    logger.info('template_max: {}'.format(template_stddev))
    logger.info('template_mean: {}'.format(template_mean))
    logger.info('template_stddev: {}'.format(template_stddev))
    logger.info('template_P95: {}'.format(template_P95))
    if template_mean+2.0*template_stddev < template_P95:
        signal_mask = np.logical_and(valid_mask, template_image>=(template_mean+2.0*template_stddev)) # select template signal > 2-sigma
        logger.info('np.count_nonzero(signal_mask): {} (2-sigma)'.format(np.count_nonzero(signal_mask)))
    else:
        signal_mask = np.logical_and(valid_mask, template_image>=template_P95) # select template signal > 2-sigma
        logger.info('np.count_nonzero(signal_mask): {} (top 5%)'.format(np.count_nonzero(signal_mask)))
    
    # median filter the mask then smooth the mask
    #filtered_mask = gaussian_filter(signal_mask, sigma=0.5) # filters out small mask areas
    filtered_mask = median_filter(signal_mask, size=5) # filters out small mask areas
    #dilated_mask = binary_dilation(filtered_mask, iterations=6) # expand areas
    dilated_mask = convolve(filtered_mask, kernel=Gaussian2DKernel(x_stddev=5))
    dilated_mask[dilated_mask<0.1] = 0
    dilated_mask[dilated_mask>=0.1] = 1
    final_mask = np.logical_and(valid_mask, dilated_mask)
    init_param_A = np.median((data_image[final_mask]) / (template_image[final_mask] + init_param_B))
    
    # set parameters for the scaling equation:
    # image = A * template + B
    # we need to fix the background difference B because we only focus to the fitting within the strict mask
    # otherwise some weird result like a negative B value can happen
    fit_params = Parameters()
    fit_params.add('A', value=init_param_A, min=0.0)
    fit_params.add('B', value=init_param_B, vary=False) # vary=False
    #fit_params.pretty_print()
    
    # run the fitting of template within the mask
    #mask = valid_mask
    #mask = signal_mask
    mask = final_mask
    x = np.arange(np.count_nonzero(mask))
    array = template_image[mask].ravel()
    minval = np.min(array)
    weight = None # (array-minval) / (np.max(array)-minval)
    fitter = Minimizer(
        template_fitting_function, 
        fit_params, 
        fcn_args = (x, ), 
        fcn_kws = dict(
            data=data_image[mask].ravel(), 
            error=error_image[mask].ravel(), 
            weight=weight,
            template=array,
        ),
        iter_cb = iter_cb, 
    )
    fit_result = fitter.minimize()
    #fit_result = fitter.minimize(method='ampgo') # too slow to find global minima
    #fit_result.params.pretty_print()
    scaled_template = fit_result.params['A'] * template_image
    scaled_template[np.isnan(scaled_template)] = 0.0
    subtracted_image = data_image - scaled_template
    return subtracted_image, scaled_template, fit_result.params['A'], fit_result.params['B'], mask



# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--template-dir', type=click.Path(exists=True), default=None, help='Where to find the downloaded NIRCam team wisps template files. If not given we will search current and home directories.')
@click.option('--template-file', type=click.Path(exists=True), default=None, help='Which exact template FITS file to use. Use with caution! This will be applied to all input images.')
def main(
        input_file, 
        output_file, 
        template_dir, 
        template_file, 
    ):
    """
    Input rate image. 
    """
    with ImageModel(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Removed wisps'):
                    logger.info('{!r} already had wisps removed. Skipping!'.format(input_file))
                    return
    
        # get image
        image = model.data.copy()
        imerr = model.err.copy()
        
        # check instrument_name
        instrument_name = model.meta.instrument.name
        if instrument_name.upper() != 'NIRCAM':
            logger.warning('The input is not NIRCam data, cannot remove wisps.')
            return
        
        # check filter, only 'F150W', 'F150W2', 'F200W', 'F210M' allowed
        filter_name = model.meta.instrument.filter
        if filter_name.upper() not in ['F150W', 'F150W2', 'F200W', 'F210M']:
            logger.warning('The input NIRCam data filter {} is not one of the {} that can remove wisps.'.format(
                filter_name, repr(['F150W', 'F150W2', 'F200W', 'F210M'])
            ))
            return
        
        # check detector, only 'NRCA3', 'NRCA4', 'NRCB3', 'NRCB4' allowed
        detector_name = model.meta.instrument.detector
        if detector_name.upper() not in ['NRCA3', 'NRCA4', 'NRCB3', 'NRCB4']:
            logger.warning('The input NIRCam data detector {} is not one of the {} that can remove wisps.'.format(
                detector_name, repr(['NRCA3', 'NRCA4', 'NRCB3', 'NRCB4'])
            ))
            return
    
    
    # find wisps template file
    if template_file is None:
        search_dirs = []
        if template_dir is None:
            if 'NIRCAM_WISP_TEMPLATES' in os.environ:
                logger.info('template_dir is not given but found system variable NIRCAM_WISP_TEMPLATES')
                search_dirs.append(os.environ['NIRCAM_WISP_TEMPLATES'])
            else:
                logger.info('template_dir is not given and system variable NIRCAM_WISP_TEMPLATES is not found, will search current dir and home dir')
                search_dirs.extend(['.', 'wisp_templates', os.path.expanduser('~/wisp_templates')])
        else:
            if not os.path.isdir(template_dir):
                errmsg = 'Error! template_dir {} does not exist?'.format(template_dir)
                logger.error(errmsg)
                raise Exception(errmsg) # template_dir not found
            search_dirs.append(template_dir)
        template_filepath = None
        searched_paths = []
        for search_dir in search_dirs:
            template_filename = 'wisps_{}_{}.fits'.format(detector_name.lower(), filter_name.upper())
            search_path = os.path.join(search_dir, template_filename)
            searched_paths.append(search_path)
            if os.path.isfile(search_path):
                template_filepath = search_path
                break
        if template_filepath is None:
            errmsg = 'Error! Could not find template file! Searched paths: {}. Please specify --template-dir.'.format(repr(searched_paths))
            logger.error(errmsg)
            raise Exception(errmsg) # template file not found
    else:
        template_filepath = template_file
    
    
    # read wisps template image
    logger.info('Reading template file: {}'.format(template_filepath))
    template_image, template_header = fits.getdata(template_filepath, header=True)
    # fix nan or invalid pixels
    if np.isfinite(template_image[0,0]) and (template_image[0,0] == template_image[1,0]):
        invalid_mask = (template_image==template_image[0,0]) # TODO: mind float == operation
        template_image[invalid_mask] = np.nan
    
    
    # fit wisps scaling factor
    subtracted_image, scaled_template, scaling_factor, scaling_baseline, mask = do_template_fitting(
        image, imerr, template_image
    )
    
    
    # output file
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_removing_wisps.fits'
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
        stepdescription = f'Removed wisps with {CODE_NAME} ({out_timestamp})'
        software_dict = {'name':CODE_NAME,
                         'author':CODE_AUTHOR,
                         'version':CODE_VERSION,
                         'homepage':CODE_HOMEPAGE}
        substr = util.create_history_entry(stepdescription, software=software_dict)
        out_model.history.append(substr)
        print('out_model.history', out_model.history)
        out_model.save(output_file)
        logger.info('Saved wisps-removed data into {!r}'.format(output_file))
    
    
    output_filepath_additional = re.sub(r'\.fits$', r'', output_file) + '_removing_wisps_scaled_template.fits'
    output_header_0 = fits.getheader(output_file, 0)
    output_header_1 = fits.getheader(output_file, 1) # this contains wcs
    output_header_1['HISTORY'] = ''
    output_header_1['HISTORY'] = 'Scaled template = A * template, where A = {}, using the template {}'.format(
                                    scaling_factor, template_filepath)
    output_header_1['HISTORY'] = ''
    output_hdu_0 = fits.PrimaryHDU(data=None, header=output_header_0)
    output_hdu_1 = fits.ImageHDU(data=scaled_template, header=output_header_1)
    output_hdulist = fits.HDUList([output_hdu_0, output_hdu_1])
    output_hdulist.writeto(output_filepath_additional, overwrite=True)
    logger.info('Saved wisps scaled template into {!r}'.format(output_filepath_additional))
    
    
    output_filepath_additional = re.sub(r'\.fits$', r'', output_file) + '_removing_wisps_scaled_template_with_bkg.fits'
    output_header_0 = fits.getheader(output_file, 0)
    output_header_1 = fits.getheader(output_file, 1) # this contains wcs
    output_header_1['HISTORY'] = ''
    output_header_1['HISTORY'] = 'Scaled template = A * template + B, where A = {} and B = {}, using the template {}'.format(
                                    scaling_factor, scaling_baseline, template_filepath)
    output_header_1['HISTORY'] = ''
    output_hdu_0 = fits.PrimaryHDU(data=None, header=output_header_0)
    output_hdu_1 = fits.ImageHDU(data=scaled_template+scaling_baseline, header=output_header_1)
    output_hdulist = fits.HDUList([output_hdu_0, output_hdu_1])
    output_hdulist.writeto(output_filepath_additional, overwrite=True)
    logger.info('Saved wisps scaled template into {!r}'.format(output_filepath_additional))
    
    
    output_filepath_additional = re.sub(r'\.fits$', r'', output_file) + '_removing_wisps_mask.fits'
    output_header_0 = fits.getheader(output_file, 0)
    output_header_1 = fits.getheader(output_file, 1) # this contains wcs
    output_hdu_0 = fits.PrimaryHDU(data=None, header=output_header_0)
    output_hdu_1 = fits.ImageHDU(data=mask.astype(int), header=output_header_1)
    output_hdulist = fits.HDUList([output_hdu_0, output_hdu_1])
    output_hdulist.writeto(output_filepath_additional, overwrite=True)
    logger.info('Saved wisps scaled template into {!r}'.format(output_filepath_additional))
    
    



if __name__ == '__main__':
    main()



