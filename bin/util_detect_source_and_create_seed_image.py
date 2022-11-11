#!/usr/bin/env python
# 
"""
Detect source emission and 2D background in the input image, output seed image.

Usage: 
    ./util_detect_source_and_create_seed_image.py some_image.fits

By Daizhong Liu @MPE. 

Last updates: 
    2022-07-27.
    2022-10-03 rewritten based on "util_make_seed_image_for_rate_image.py".
    2022-11-11 

"""
import os, sys, re, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.modeling import models as apy_models
from astropy.modeling import fitting as apy_fitting
from astropy.convolution import convolve as apy_convolve
from astropy.convolution import Gaussian2DKernel
from photutils.background import Background2D
from scipy.ndimage import median_filter as scipy_median_filter
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# mpl.rcParams['savefig.dpi'] = 300
# mpl.rcParams['figure.dpi'] = 300

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

# code name and version
CODE_NAME = 'util_detect_source_and_create_seed_image.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221111'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)



# define functions
def prep_crds_getreferences_kwargs(model):
    crds_getreferences_kwargs = {}
    
    import crds
    crds_context = None
    if 'CRDS_CONTEXT' in os.environ: 
        if os.environ['CRDS_CONTEXT'] != '':
            crds_context = os.environ['CRDS_CONTEXT']
    if crds_context is None:
        crds_context = crds.get_default_context()
    
    # CRDS query parameters
    # https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/reference_files.html
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
        raise ValueError("Error! The input model (model.meta.instrument.name) is not NIRCAM or MIRI!")
    
    crds_getreferences_kwargs['parameters'] = crds_dict
    crds_getreferences_kwargs['reftypes'] = ['flat']
    crds_getreferences_kwargs['context'] = crds_context
    
    return crds_getreferences_kwargs



def get_flat_from_CRDS(model, raise_exception=True):
    logger.info('Querying flat from CRDS')
    import crds
    
    if isinstance(model, str):
        with datamodels.open(model) as model:
            crds_getreferences_kwargs = prep_crds_getreferences_kwargs(model)
    else:
        crds_getreferences_kwargs = prep_crds_getreferences_kwargs(model)
    logger.info('crds_getreferences_kwargs: {}'.format(crds_getreferences_kwargs)) #<DZLIU>#
    flats = crds.getreferences(**crds_getreferences_kwargs)
    
    # check if CRDS got the flat correctly
    try:
        flatfile = flats['flat']
    except KeyError:
        errmsg = 'Flat was not found in CRDS with the parameters: {}'.format(crds_dict)
        logger.error(errmsg)
        if raise_exception:
            raise Exception(errmsg)
        else:
            return None
    
    return flatfile
    


def apply_flat_to_model(model):
    
    flatfile = get_flat_from_CRDS(model)
    
    logger.info('Applying flat: %s'%(flatfile))
    
    with FlatModel(flatfile) as flat:
        # use the JWST Calibration Pipeline flat fielding Step 
        model, applied_flat = do_correction(model, flat)
    
    return model, applied_flat



def get_dqmask_from_dq(dq, image=None):
    dqmask = bitfield_to_boolean_mask(
        dq,
        interpret_bit_flags(
            '~(DO_NOT_USE,SATURATED,NON_SCIENCE,JUMP_DET,OUTLIER,PERSISTENCE,AD_FLOOR,RESERVED_4,UNRELIABLE_FLAT,OTHER_BAD_PIXEL)', # SkyMatchStep().dqbits is '~DO_NOT_USE+NON_SCIENCE' in default
            flag_name_map=dqflags_pixel), # -1073742468
        good_mask_value=1,
        dtype=int
        #dtype=np.uint8
    ) # >0 means good data, following "jwst/skymatch/skymatch_step.py"
    # -2147090176
    if image is not None:
        dqmask *= np.isfinite(image)
    return dqmask



def get_dqmask_for_image(
        fits_image, 
        combine_flat_dqmask = True, # only valid if the input is not i2d but a rate or cal image.
        flat_file = None, 
    ):
    """The input image can be a rate image or an i2d image. 
    """
    
    with fits.open(fits_image) as hdul:
        dq = None
        dqmask = None
        extension = ''
        
        if extension == '':
            try:
                dq = hdul['DQ'].data
                extension = 'DQ'
            except KeyError:
                dq = None
        
        if extension == '':
            try:
                context_image = hdul['CON'].data
                while len(context_image.shape) > 2:
                    context_image = context_image[0]
                dqmask = (context_image != 0).astype(int)
                extension = 'CON'
                # CON: 2-D context image, which encodes information about 
                # which input images contribute to a specific output pixel
                # -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/science_products.html
                # context 0 means no input image for the output mosaic pixel
            except KeyError:
                dq = None
        
        if dqmask is None and dq is not None:
            dqmask = get_dqmask_from_dq(dq)
    
    if extension == '' and combine_flat_dqmask:
        # 
        #DZLIU: 2022-07-27 the problem is that MIRI cal.fits DQ are not distinguishing 4QPM or coronograph support shadow areas!
        #                  the DQ flags are only correct in the 'flat' calibration file "jwst_miri_flat_0786.fits"
        # 
        #if dq is None:
        #    raise Exception('The option apply_flat should be set only when the input data has a DQ extension!')
        if flat_file is None:
            flat_file = get_flat_from_CRDS(fits_image)
        logger.info('Combining DQ flags from flat reference file: {}'.format(flat_file))
        with FlatModel(flat_file) as flat:
            if dq is None:
                dq = flat.dq
            else:
                dq = np.bitwise_or(dq, flat.dq)
        dqmask = get_dqmask_from_dq(dq)
    
    return dqmask



def detect_source_and_background_for_image(
        fits_image, 
        lsigma = None, 
        usigma = None, 
        sigma = 3.0, 
        box_size = None, 
        box_frac = 0.05, 
        median_filter = 1, 
        smooth = 0.0, 
        minpixarea = 1, # 
        flat_file = None, # for MIRI only, will combine DQ from the flat file.
        with_filter_in_output_name = False, # with a filter name in the output filename
        output_dir = None, 
        output_suffix = '_galaxy_seed_image', 
        overwrite = False, 
        verbose = True, 
    ):
    
    if verbose:
        logger.info('Detecting source and background for image: {!r}'.format(fits_image))
    
    # strip off file name suffix
    fits_name = os.path.splitext(os.path.basename(fits_image))[0]
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(fits_image))
    elif not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    # read image
    hdulist = fits.open(fits_image)
    pheader = hdulist[0].header
    hdulist.close()
    image, header = fits.getdata(fits_image, header=True)
    if len(image.shape) == 0 or image.shape == (0,):
        image, header1 = fits.getdata(fits_image, ext=1, header=True)
        for key in header1:
            header[key] = header1[key]
    
    # prepare output file name
    if with_filter_in_output_name:
        if 'PUPIL' in pheader:
            output_suffix = '_' + pheader['INSTRUME'].strip() + \
                            '_' + pheader['FILTER'].strip() + \
                            '_' + pheader['PUPIL'].strip() + \
                            output_suffix
        else:
            output_suffix = '_' + pheader['INSTRUME'].strip() + \
                            '_' + pheader['FILTER'].strip() + \
                            output_suffix
    output_fits_image = output_dir + os.sep + fits_name + output_suffix + '.fits' # this extension is needed by 'remstriping.py'
    
    # check output file
    if verbose:
        logger.info('Output image file: {!r}'.format(output_fits_image))
    if os.path.isfile(output_fits_image): 
        if not overwrite:
            if verbose:
                logger.info('Found existing output file: "{}" and overwrite is set to False. Will do nothing.'.format(output_fits_image))
            return
        else:
            shutil.move(output_fits_image, output_fits_image+'.backup')
    
    
    # get dqmask
    combine_flat_dqmask = False
    if pheader['INSTRUME'].strip().upper() == 'MIRI':
        combine_flat_dqmask = True
    elif flat_file is not None:
        combine_flat_dqmask = True
    dqmask = get_dqmask_for_image(fits_image, combine_flat_dqmask=combine_flat_dqmask, flat_file=flat_file)
    if dqmask is None:
        dqmask = np.isfinite(image).astype(int)
    else:
        dqmask *= np.isfinite(image).astype(int)
    
    
    # get valid pixels, write to disk
    valid_mask = (dqmask > 0)
    output_mask_fits_image = re.sub(r'\.fits$', '_valid_pixels.fits', output_fits_image) # 1 means valid pixel
    output_hdu = fits.PrimaryHDU(data=valid_mask.astype(int), header=header)
    output_hdu.writeto(output_mask_fits_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_mask_fits_image))
    
    valid_data = image[valid_mask].ravel()
    pixval_min = np.nanmin(valid_data)
    pixval_max = np.nanmax(valid_data)
    
    # analyze Background2D with a 1/10 box size
    box_size = min(header['NAXIS1'], header['NAXIS2']) // 10
    bkg = Background2D(image, box_size = box_size, 
                       mask = ~valid_mask, # mask means invalid here. dqmask 1 means valid.
                       exclude_percentile = 50.0, # default is 10.0 percent bad pixel per box
                      ) 
    pixval_median = bkg.background_median
    pixval_rms = np.nanmean(bkg.background_rms)
    bkg2d_map = bkg.background
    rms2d_map = bkg.background_rms
    if verbose:
        logger.info('Computed background median {} rms {}'.format(pixval_median, pixval_rms))
    
    
    # get lsigma, usigma
    if lsigma is None and sigma is not None:
        lsigma = sigma
    if usigma is None and sigma is not None:
        usigma = sigma
    
    if verbose:
        logger.info('Sigma-clipping with lsigma {} usigma {}'.format(lsigma, usigma))
    
    
    # improve the mask
    #nonbright_mask = np.logical_and(image > pixval_median - lsigma * pixval_rms, 
    #                                image < pixval_median + usigma * pixval_rms)
    nonbright_mask = np.logical_and(image > bkg2d_map - lsigma * pixval_rms, 
                                    image < bkg2d_map + usigma * pixval_rms)
    
    
    # analyze Background2D with a 1/30 box size
    box_size = min(header['NAXIS1'], header['NAXIS2']) // 30
    bkg2 = Background2D(image, box_size = box_size, 
                       mask = ~np.logical_and(valid_mask, nonbright_mask), # mask means invalid here. 
                       exclude_percentile = 50.0, # default is 10.0 percent bad pixel per box
                      ) 
    pixval_median = bkg2.background_median
    pixval_rms = np.nanmean(bkg2.background_rms)
    bkg2d_map = bkg2.background
    rms2d_map = bkg2.background_rms
    if verbose:
        logger.info('Computed background median {} rms {}'.format(pixval_median, pixval_rms))
    
    output_bkg2d_image = re.sub(r'\.fits$', '_bkg2d.fits', output_fits_image) # 1 means valid pixel
    output_hdu = fits.PrimaryHDU(data=bkg2d_map, header=header)
    output_hdu.writeto(output_bkg2d_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_bkg2d_image))
    
    output_rms2d_image = re.sub(r'\.fits$', '_rms2d.fits', output_fits_image) # 1 means valid pixel
    output_hdu = fits.PrimaryHDU(data=rms2d_map, header=header)
    output_hdu.writeto(output_rms2d_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_rms2d_image))
    
    
    # update the mask
    #bright_mask = (image > bkg2d_map + usigma * rms2d_map)
    bright_mask = (image > bkg2d_map + usigma * pixval_rms)
    
    
    # mask sources by 3-sigma and invalid data
    #mask = np.logical_or(bright_mask, dqmask == 0) # mask bright sources (>3sigma) and invalid pixels (dqmask == 0)
    mask = bright_mask # mask bright sources (>3sigma), IGNORING invalid pixels (dqmask == 0)
    arr = mask.astype(int)
    
    
    # median filter
    if median_filter > 0:
        if verbose:
            logger.info('Median-filtering pixels with size {} pixels'.format(median_filter))
        arr = scipy_median_filter(arr, median_filter)
    
    
    # minpixarea
    if minpixarea > 0:
        header['MPIXAREA'] = (minpixarea, 'Mininum pixel area.')
        if verbose:
            logger.info('Filtering pixels with at least {} surrounding valid pixels'.format(minpixarea))
        ny, nx = arr.shape
        arrbig = np.zeros((2+ny, 2+nx)).astype(int)
        arrbig[1:-1, 1:-1] = arr[:, :]
        arrsum = arrbig[0:-2, 0:-2] + \
                 arrbig[0:-2, 1:-1] + \
                 arrbig[0:-2, 2:  ] + \
                 arrbig[1:-1, 0:-2] + \
                 arrbig[1:-1, 2:  ] + \
                 arrbig[2:  , 0:-2] + \
                 arrbig[2:  , 1:-1] + \
                 arrbig[2:  , 2:  ]
        arrmsk = (arrsum >= minpixarea).astype(int)
        arr = arr * arrmsk
    
    
    arr = arr.astype(float)
    
    # smoothing/expanding the mask
    if smooth > 0.0:
        header['SMOOTH'] = (smooth, 'Gaussian2DKernel stddev in pixel')
        if verbose:
            logger.info('Smoothing with Gaussian2DKernel kernel stddev {}'.format(smooth))
        kernel = Gaussian2DKernel(x_stddev=smooth) # 1 pixel stddev kernel
        arr = apy_convolve(arr, kernel)
        arr[arr<1e-6] = 0.0
    
    # save fits file
    hdu = fits.PrimaryHDU(data=arr, header=header)
    hdu.writeto(output_fits_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_fits_image))
    
    # save masked image
    output_masked_image = re.sub(r'\.fits$', '_nonbrigthmask.fits', output_fits_image) # 1 means valid pixel
    output_hdu = fits.PrimaryHDU(data=nonbright_mask.astype(int), header=header)
    output_hdu.writeto(output_masked_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_masked_image))
    
    # save masked image
    masked_image = np.full(image.shape, fill_value=0.0)
    #mask1 = np.logical_and(arr>0, dqmask > 0)
    mask1 = (arr>0)
    masked_image[mask1] = image[mask1]
    output_masked_image = re.sub(r'\.fits$', '_masked.fits', output_fits_image) # 1 means valid pixel
    output_hdu = fits.PrimaryHDU(data=masked_image, header=header)
    output_hdu.writeto(output_masked_image, overwrite=True)
    if verbose:
        logger.info('Output to "{}"'.format(output_masked_image))




@click.command()
@click.argument('fits_image', type=click.Path(exists=True))
@click.option('--lsigma', type=float, default=None)
@click.option('--usigma', type=float, default=None)
@click.option('--sigma', type=float, default=3.0)
@click.option('--box-size', type=int, default=None)
@click.option('--box-frac', type=float, default=0.05)
@click.option('--median-filter', type=int, default=1) # median_filter
@click.option('--minpixarea', type=int, default=1)
@click.option('--smooth', type=float, default=0.0)
@click.option('--flat-file', type=click.Path(exists=True), default=None)
@click.option('--output-dir', type=click.Path(exists=False), default=None)
@click.option('--with-filter-in-output-name', is_flag=True, default=False)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
@click.option('--verbose/--no-verbose', is_flag=True, default=True)
def main(
        fits_image, 
        lsigma, 
        usigma, 
        sigma, 
        box_size,
        box_frac,
        median_filter, 
        minpixarea, 
        smooth, 
        flat_file, 
        output_dir,
        with_filter_in_output_name, 
        overwrite,
        verbose, 
    ):
    
    detect_source_and_background_for_image(
        fits_image, 
        lsigma = lsigma, 
        usigma = usigma, 
        sigma = sigma, 
        box_size = box_size, 
        box_frac = box_frac, 
        median_filter = median_filter, 
        minpixarea = minpixarea, 
        smooth = smooth, 
        flat_file = flat_file, 
        output_dir = output_dir, 
        with_filter_in_output_name = with_filter_in_output_name, 
        overwrite = overwrite, 
        verbose = verbose, 
    )




if __name__ == '__main__':
    
    main()



