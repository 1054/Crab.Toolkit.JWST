#!/usr/bin/env python
# 
"""
Finding source emission in the image and deriving pixel histogram and 1D fitting.

Usage: 
    ./util_make_seed_image_for_rate_image.py ngc0628_miri_lvl3_f770w_i2d_hackedbgr.fits --min=-1.0 --max=4.0

By Daizhong Liu @MPE. 

Last update: 2022-07-27.
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
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300

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


def make_seed_image_for_rate_image(
        fits_image, 
        sigma = 3.0, 
        smooth = 1.0, 
        nbin = 91, 
        bin_min = None, 
        bin_max = None, 
        fit_min = None, 
        fit_max = None, 
        fit_half = True, 
        fit_core = False, 
        dynamical_range = 2.0, 
        check_rate_image = False, 
        overwrite = True, 
    ):
    
    # check fits_image
    if check_rate_image:
        if not fits_image.endswith('_rate.fits') and not fits_image.endswith('_cal.fits'):
            raise Exception('Error! The input FITS image file name should ends with \"_rate.fits\"!')
        fits_name = re.sub(r'^(.*)(_rate|_cal)\.fits$', r'\1', fits_image)
    else:
        fits_name = os.path.splitext(fits_image)[0]
    
    # read image
    hdulist = fits.open(fits_image)
    pheader = hdulist[0].header
    hdulist.close()
    image, header = fits.getdata(fits_image, header=True)
    
    # check output file
    if 'PUPIL' in pheader:
        output_fits_image = fits_name + '_{}_{}_{}_galaxy_seed_image.fits'.format(
            pheader['INSTRUME'].strip(), pheader['FILTER'].strip(), pheader['PUPIL'].strip())
    else:
        output_fits_image = fits_name + '_{}_{}_galaxy_seed_image.fits'.format(
            pheader['INSTRUME'].strip(), pheader['FILTER'].strip())
    
    if os.path.isfile(output_fits_image): 
        if not overwrite:
            print('Found existing output file: "{}" and overwrite is set to False. Will do nothing.'.format(output_fits_image))
            return
        else:
            shutil.move(output_fits_image, output_fits_image+'.backup')
    
    output_histogram_figure = re.sub(r'^(.*)_galaxy_seed_image.fits$', r'\1_histogram.pdf', output_fits_image)
    if os.path.isfile(output_histogram_figure): 
        shutil.move(output_histogram_figure, output_histogram_figure+'.backup')
    
    
    # deal with masked data using the error map
    if fits_image.find('_i2d') >= 0:
        #error = fits.getdata(fits_image, header=False, extname='ERR')
        #dqmask = np.isfinite(error).astype(int)
        context_image = fits.getdata(fits_image, header=False, extname='CON')
        dqmask = context_image
    # deal with masked data using jwst.datamodels
    else:
        image_model = ImageModel(
                    fits_image
                )
        # 
        #DZLIU: 2022-07-27 the problem is that MIRI cal.fits DQ are not distinguishing 4QPM or coronograph support shadow areas!
        #                  the DQ flags are only correct in the 'flat' calibration file "jwst_miri_flat_0786.fits"
        # 
        instrument_name = image_model.meta.instrument.name
        if instrument_name.upper() == 'MIRI':
            crds_dict = {'INSTRUME':instrument_name, #<DZLIU>#
                         'DETECTOR':image_model.meta.instrument.detector, 
                         'FILTER':image_model.meta.instrument.filter, 
                         'BAND':image_model.meta.instrument.band, 
                         'READPATT':image_model.meta.exposure.readpatt, 
                         'SUBARRAY':image_model.meta.subarray.name, 
                         'DATE-OBS':image_model.meta.observation.date,
                         'TIME-OBS':image_model.meta.observation.time}
            print('crds_dict:', crds_dict) #<DZLIU>#
            import crds
            try:
                crds_context = os.environ['CRDS_CONTEXT']
            except KeyError:
                crds_context = crds.get_default_context()
            flats = crds.getreferences(crds_dict, reftypes=['flat'], 
                                       context=crds_context)
            # if the CRDS loopup fails, should return a CrdsLookupError, but 
            # just in case:
            try:
                flatfile = flats['flat']
            except KeyError:
                print('Flat reference file was not found in CRDS with the parameters: {}'.format(crds_dict))
                exit()
            
            log.info('Combining DQ flags from flat reference file: %s'%(os.path.basename(flatfile)))
            with FlatModel(flatfile) as flat:
                image_model.dq = np.bitwise_or(image_model.dq, flat.dq)
        # 
        dqmask = bitfield_to_boolean_mask(
                    image_model.dq,
                    interpret_bit_flags('~(DO_NOT_USE,NON_SCIENCE,SATURATED,RESERVED_4,OTHER_BAD_PIXEL)', # SkyMatchStep().dqbits is '~DO_NOT_USE+NON_SCIENCE' in default
                        flag_name_map=dqflags_pixel), 
                    good_mask_value=1,
                    dtype=np.uint8
                ) * np.isfinite(image_model.data) # following "jwst/skymatch/skymatch_step.py", this is valid data
        #print('dqmask', dqmask, np.count_nonzero(dqmask))
        #DZLIU: https://github.com/spacetelescope/jwst/wiki/DQ-Flags
        #DZLIU: https://jwst-reffiles.stsci.edu/source/data_quality.html
        #DZLIU: checked MIRI "jwst_crds_cache/references/jwst/miri/jwst_miri_flat_0786.fits" 
        # DQ flag: struct.pack('>i', -2147483648).hex() = 80 (10000000) 00 00 00 (00000000)  -- science data
        # DQ flag: struct.pack('>i', -2147483581).hex() = 80 (10000000) 00 00 43 (10000011)  -- miri blank area
        # DQ flag: struct.pack('>i', -2147483583).hex() = 80 (10000000) 00 00 41 (10000001)  -- miri Lyot coronagraph support rack
        # DQ flag: struct.pack('>i', -2147483645).hex() = 80 (10000000) 00 00 03 (00000011)  -- miri Lyot coronagraph support rack
        
    
    
    # get valid pixels
    valid_mask = (dqmask > 0)
    output_mask_fits_image = re.sub(r'\.fits$', '_mask.fits', output_fits_image)
    fits.PrimaryHDU(data=valid_mask.astype(int), header=header).writeto(output_mask_fits_image, overwrite=True)
    print('Output to "{}"'.format(output_mask_fits_image))
    
    valid_data = image[valid_mask].ravel()
    pixval_min = np.nanmin(valid_data)
    pixval_max = np.nanmax(valid_data)
    
    # analyze pixel histogram
    nbin = int(np.ceil(float(nbin)/2.)*2+1) # 50
    dynamical_range = dynamical_range # 2.0 #<TODO>#
    clip_sigma = sigma # 1.0 #<TODO>
    conv_kern = smooth # 0 #<TODO># do a convolution to the not-used-for-background mask
    
    which_method = 2
    if which_method == 1:
        # method 1, too slow
        pixval_median = np.nanpercentile(valid_data, 50)
        if bin_min is None:
            bin_min = max(1.0/dynamical_range*pixval_median, pixval_min)
        if bin_max is None:
            bin_max = min(dynamical_range*pixval_median, pixval_max)
        hist, bin_edges = np.histogram(valid_data, bins=nbin, range=(bin_min, bin_max))
        bin_centers = (bin_edges[1:]+bin_edges[0:-1])/2.0
        
        recheck_valid_hist = 0
        while np.count_nonzero(hist>0) < 5 and recheck_valid_hist < 10:
            # insufficient histograms, shrink the dynamical range
            shrink_factor = float(np.count_nonzero(hist>0)+1) / float(len(hist))
            dynamical_range = dynamical_range * shrink_factor
            bin_min = max(1.0/dynamical_range*pixval_median, pixval_min)
            bin_max = min(dynamical_range*pixval_median, pixval_max)
            hist, bin_edges = np.histogram(valid_data, bins=nbin, range=(bin_min, bin_max))
            bin_centers = (bin_edges[1:]+bin_edges[0:-1])/2.0
            recheck_valid_hist += 1
        
    elif which_method == 2:
        # method 2
        box_size = min(header['NAXIS1'], header['NAXIS2']) // 10
        bkg = Background2D(image, box_size = box_size, mask = (dqmask == 0)) # mask means invalid here. dqmask 1 means valid.
        pixval_median = bkg.background_median
        pixval_rms = np.nanmax(bkg.background_rms)
        if bin_min is None:
            bin_min = max((pixval_median - 7.*pixval_rms), pixval_min)
        if bin_max is None:
            bin_max = min((pixval_median + 7.*pixval_rms), pixval_max)
        print('pixval_median:', pixval_median)
        print('pixval_rms:', pixval_rms)
        print('bin_min:', bin_min)
        print('bin_max:', bin_max)
        hist, bin_edges = np.histogram(valid_data, bins=np.linspace(bin_min, bin_max, num=nbin+1, endpoint=True))
        print('bin_edges:', bin_edges)
        print('hist:', hist)
        bin_centers = (bin_edges[1:]+bin_edges[0:-1])/2.0
        dynamical_range_1 = pixval_median / bin_min
        dynamical_range_2 = bin_max / pixval_median
        dynamical_range = max(dynamical_range_1, dynamical_range_2)
        
    else:
        raise NotImplementedError('which_method {} is invalid.'.format(which_method))
    
    
    # plot pixel histogram
    fig, ax = plt.subplots(ncols=1, nrows=1)
    ax.bar(bin_edges[0:-1], hist, width=(bin_edges[1:]-bin_edges[0:-1]), alpha=0.8, color='C0', align='edge')
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim([ax.get_ylim()[0], ax.get_ylim()[1]*1.09]) # expanded a little bit the ylim
    
    #ax.plot([pixval_median]*2, ax.get_ylim(), ls='dashed', color='C1')
    
    # get histogram peak (mode)
    pixval_argmode = np.argwhere(hist==np.max(hist)).ravel()[0]
    pixval_mode = bin_centers[pixval_argmode]
    
    # prepare fitting range
    if fit_min is None:
        fit_min = np.min(bin_centers)
    if fit_max is None:
        fit_max = np.max(bin_centers)
    if fit_half:
        # fit only the lower half of the histogram
        if pixval_argmode < len(bin_centers)-1:
            # slightly extending to one bin over the mode (peak of the histogram)
            fit_max = bin_centers[pixval_argmode+1]
        else:
            fit_max = pixval_mode
    
    # model fitting 1D Gaussian histogram
    model_init = apy_models.Gaussian1D(amplitude=np.max(hist), mean=pixval_median, stddev=pixval_median*0.2)
    fitter = apy_fitting.LevMarLSQFitter()
    print('fitting with bins from min {} to max {}'.format(fit_min, fit_max))
    bin_mask = np.logical_and(bin_centers >= fit_min, bin_centers <= fit_max)
    model_fitted = fitter(model_init, bin_centers[bin_mask], hist[bin_mask])
    hist_fitted = model_fitted(bin_centers)
    pixval_mean_fitted = model_fitted.mean.value
    pixval_stddev_fitted = model_fitted.stddev.value
    pixval_1sigma = pixval_mean_fitted + 1.0 * pixval_stddev_fitted
    pixval_3sigma = pixval_mean_fitted + 3.0 * pixval_stddev_fitted
    
    # refit with the core part of the histogram
    if fit_core and np.count_nonzero(valid_data > pixval_1sigma) > 0:
       print('refitting with bins from min {} to 1-sigma {}'.format(fit_min, pixval_1sigma))
       bin_mask = np.logical_and(bin_mask, bin_centers < pixval_1sigma)
       model_fitted = fitter(model_init, bin_centers[bin_mask], hist[bin_mask])
       hist_fitted = model_fitted(bin_centers)
       pixval_mean_fitted = model_fitted.mean.value
       pixval_stddev_fitted = model_fitted.stddev.value
       pixval_1sigma = pixval_mean_fitted + 1.0 * pixval_stddev_fitted
       pixval_3sigma = pixval_mean_fitted + 3.0 * pixval_stddev_fitted
    
    # plot vertical lines
    ax.plot([pixval_mean_fitted]*2, ax.get_ylim(), ls='dotted', color='red')
    ax.text(ax.get_xlim()[0] + 0.02 * (ax.get_xlim()[1]-ax.get_xlim()[0]), 
            ax.get_ylim()[0] + 0.97 * (ax.get_ylim()[1]-ax.get_ylim()[0]), 
            'mean={:.8g}, \nstddev={:.8g}'.format(pixval_mean_fitted, pixval_stddev_fitted), 
            ha='left', va='top', color='red')
    
    ax.plot([pixval_3sigma]*2, ax.get_ylim(), ls='dotted', color='C0', alpha=0.9)
    ax.text(pixval_3sigma + 0.02 * (ax.get_xlim()[1]-ax.get_xlim()[0]), 
            ax.get_ylim()[0] + 0.97 * (ax.get_ylim()[1]-ax.get_ylim()[0]), 
            '3$\\sigma$={:.8g}'.format(pixval_3sigma), 
            ha='left', va='top', color='C0')
    
    #header['PIXDYN'] = (dynamical_range, 'dynamical range')
    header['PIXMIN'] = bin_min
    header['PIXMAX'] = bin_max
    header['FITMIN'] = fit_min
    header['FITMAX'] = fit_max
    header['PIXMODE'] = pixval_mode
    header['PIXMED'] = pixval_median
    header['PIXMEAN'] = pixval_mean_fitted
    header['PIXSTD'] = pixval_stddev_fitted
    header['PIX1SIG'] = pixval_1sigma
    header['PIX3SIG'] = pixval_3sigma
    if conv_kern > 0.0:
        header['CONVKERN'] = (conv_kern, 'Gaussian2DKernel stddev in pixel')
    
    # add model fitting line into the plot
    ax.plot(bin_centers, hist_fitted, ls='solid', color='C3')
    ax.grid(True, color='#aaaaaa', alpha=0.5, lw=0.5, ls='dotted')
    
    # save figure
    fig.savefig(output_histogram_figure)
    print('Output to "{}"'.format(output_histogram_figure))
    
    # get reduced chisq
    chisq_fitted = (hist_fitted - hist) / np.max(hist) / nbin
    
    # mask sources by 3-sigma and invalid data
    mask = np.logical_or(image > pixval_3sigma, dqmask == 0)
    arr = mask.astype(int).astype(float)
    if conv_kern > 0.0:
        kernel = Gaussian2DKernel(x_stddev=conv_kern) # 1 pixel stddev kernel
        arr = apy_convolve(arr, kernel)
        arr[arr<1e-6] = 0.0
    
    # save fits file
    hdu = fits.PrimaryHDU(data=arr, header=header)
    hdu.writeto(output_fits_image)
    print('Output to "{}"'.format(output_fits_image))




@click.command()
@click.argument('fits_image', type=click.Path(exists=True))
@click.option('--sigma', type=float, default=3.0)
@click.option('--smooth', type=float, default=1.0)
@click.option('--nbin', type=int, default=91)
@click.option('--min', 'bin_min', type=float, default=None, help='range to get the histogram')
@click.option('--max', 'bin_max', type=float, default=None, help='range to get the histogram')
@click.option('--fit-min', type=float, default=None, help='range to fit the histogram')
@click.option('--fit-max', type=float, default=None, help='range to fit the histogram')
@click.option('--fit-half/--no-fit-half', type=bool, default=True, is_flag=True, help='fit only the left half of the histogram?')
@click.option('--fit-core/--no-fit-core', type=bool, default=False, is_flag=True, help='fit only the core part (min to +1sigma) of the histogram?')
@click.option('--dynamical-range', type=float, default=2.0)
@click.option('--check-rate-image/--no-check-rate-image', is_flag=True, default=False)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=True)
def main(fits_image, sigma, smooth, nbin, 
         bin_min, 
         bin_max, 
         fit_min, 
         fit_max, 
         fit_half, 
         fit_core, 
         dynamical_range, check_rate_image, overwrite):
    
    make_seed_image_for_rate_image(
        fits_image, 
        sigma = sigma, 
        smooth = smooth, 
        nbin = nbin, 
        bin_min = bin_min, 
        bin_max = bin_max, 
        fit_min = fit_min, 
        fit_max = fit_max, 
        fit_half = fit_half, 
        fit_core = fit_core, 
        dynamical_range = dynamical_range, 
        check_rate_image = check_rate_image, 
        overwrite = overwrite, 
    )



if __name__ == '__main__':
    main()

    
