#!/usr/bin/env python
# 
import os, sys, re, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.modeling import models as apy_models
from astropy.modeling import fitting as apy_fitting
from astropy.convolution import convolve as apy_convolve
from astropy.convolution import Gaussian2DKernel
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

from jwst.datamodels import ImageModel


def make_seed_image_for_rate_image(
        fits_image, 
        sigma = 3.0, 
        smooth = 1.0, 
        nbin = 50, 
        dynamical_range = 2.0, 
        check_rate_image = True, 
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
    
    # deal with masked data
    image_model = ImageModel(
                fits_image
            )
    dqmask = bitfield_to_boolean_mask(
                image_model.dq,
                interpret_bit_flags('~DO_NOT_USE+NON_SCIENCE', # SkyMatchStep().dqbits is '~DO_NOT_USE+NON_SCIENCE' in default
                    flag_name_map=dqflags_pixel), 
                good_mask_value=1,
                dtype=np.uint8
            ) * np.isfinite(image_model.data) # following "jwst/skymatch/skymatch_step.py", this is valid data
    
    #print('dqmask', dqmask, np.count_nonzero(dqmask))
    
    # get valid pixels
    valid_data = image[dqmask].ravel()
    
    # analyze pixel histogram
    nbin = nbin # 50
    dynamical_range = dynamical_range # 2.0 #<TODO>#
    clip_sigma = sigma # 1.0 #<TODO>
    conv_kern = smooth # 0 #<TODO># do a convolution to the not-used-for-background mask
    pixval_median = np.nanpercentile(valid_data, 50)
    pixval_min = max(1.0/dynamical_range*pixval_median, np.nanmin(valid_data))
    pixval_max = min(dynamical_range*pixval_median, np.nanmax(valid_data))
    hist, bin_edges = np.histogram(valid_data, bins=nbin, range=(pixval_min, pixval_max))
    bin_centers = (bin_edges[1:]+bin_edges[0:-1])/2.0
    
    # plot pixel histogram
    fig, ax = plt.subplots(ncols=1, nrows=1)
    ax.bar(bin_edges[0:-1], hist, width=(bin_edges[1:]-bin_edges[0:-1]), alpha=0.8, color='C0', align='edge')
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())
    ax.plot([pixval_median]*2, ax.get_ylim(), ls='dashed', color='C1')
    
    # model fitting 1D Gaussian histogram
    model_init = apy_models.Gaussian1D(amplitude=np.max(hist), mean=pixval_median, stddev=pixval_median*0.2)
    fitter = apy_fitting.LevMarLSQFitter()
    model_fitted = fitter(model_init, bin_centers, hist)
    hist_fitted = model_fitted(bin_centers)
    pixval_mean_fitted = model_fitted.mean.value
    pixval_stddev_fitted = model_fitted.stddev.value
    pixval_3sigma = pixval_mean_fitted + 3.0 * pixval_stddev_fitted
    ax.plot([pixval_mean_fitted]*2, ax.get_ylim(), ls='dotted', color='C2')
    header['PIXMED'] = pixval_median
    header['PIXMEAN'] = pixval_mean_fitted
    header['PIXSTD'] = pixval_stddev_fitted
    header['PIX3SIG'] = pixval_3sigma
    if conv_kern > 0.0:
        header['CONVKERN'] = (conv_kern, 'Gaussian2DKernel stddev in pixel')
    
    # add model fitting line into the plot
    ax.plot(bin_centers, hist_fitted, ls='solid', color='C3')
    
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
@click.option('--nbin', type=int, default=50)
@click.option('--dynamical-range', type=float, default=2.0)
@click.option('--check-rate-image/--no-check-rate-image', default=True)
@click.option('--overwrite/--no-overwrite', default=True)
def main(fits_image, sigma, smooth, nbin, dynamical_range, check_rate_image, overwrite):
    
    make_seed_image_for_rate_image(
        fits_image, 
        sigma = sigma, 
        smooth = smooth, 
        nbin = nbin, 
        dynamical_range = dynamical_range, 
        check_rate_image = check_rate_image, 
        overwrite = overwrite, 
    )



if __name__ == '__main__':
    main()

    
