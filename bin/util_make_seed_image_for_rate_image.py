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



def make_seed_image_for_rate_image(fits_image, check_rate_image = True, overwrite = True):
    
    # check fits_image
    if check_rate_image:
        if not fits_image.endswith('_rate.fits'):
            raise Exception('Error! The input FITS image file name should ends with \"_rate.fits\"!')
        fits_name = re.sub(r'^(.*)_rate\.fits$', r'\1', fits_image)
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
    
    # analyze pixel histogram
    nbin = 40
    dynamic_range = 2.0
    pixval_median = np.nanpercentile(image, 50)
    pixval_min = max(1.0/dynamic_range*pixval_median, np.nanmin(image))
    pixval_max = min(dynamic_range*pixval_median, np.nanmax(image))
    hist, bin_edges = np.histogram(image.ravel(), bins=nbin, range=(pixval_min, pixval_max))
    bin_centers = (bin_edges[1:]+bin_edges[0:-1])/2.0
    
    # plot pixel histogram
    fig, ax = plt.subplots(ncols=1, nrows=1)
    ax.bar(bin_edges[0:-1], hist, width=(bin_edges[1:]-bin_edges[0:-1]), alpha=0.8, color='C0', align='edge')
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())
    ax.plot([pixval_median]*2, ax.get_ylim(), ls='dashed', color='C1')
    
    # model fitting
    model_init = apy_models.Gaussian1D(amplitude=np.max(hist), mean=pixval_median, stddev=pixval_median*0.2)
    fitter = apy_fitting.LevMarLSQFitter()
    model_fitted = fitter(model_init, bin_centers, hist)
    hist_fitted = model_fitted(bin_centers)
    pixval_mean_fitted = model_fitted.mean.value
    pixval_stddev_fitted = model_fitted.stddev.value
    pixval_3sigma = pixval_mean_fitted + 3.0 * pixval_stddev_fitted
    header['PIXMED'] = pixval_median
    header['PIXMEAN'] = pixval_mean_fitted
    header['PIXSTD'] = pixval_stddev_fitted
    header['PIX3SIG'] = pixval_3sigma
    conv_kern = 1
    header['CONVKERN'] = (conv_kern, 'Gaussian2DKernel stddev in pixel')
    
    # add model fitting line into the plot
    ax.plot(bin_centers, hist_fitted, ls='solid', color='C3')
    
    # save figure
    fig.savefig(output_histogram_figure)
    print('Output to "{}"'.format(output_histogram_figure))
    
    # get reduced chisq
    chisq_fitted = (hist_fitted - hist) / np.max(hist) / nbin
    
    # mask sources by 3-sigma
    mask = (image > pixval_3sigma)
    arr = mask.astype(int).astype(float)
    kernel = Gaussian2DKernel(x_stddev=conv_kern) # 1 pixel stddev kernel
    arr = apy_convolve(arr, kernel)
    arr[arr<1e-6] = 0.0
    
    # save fits file
    hdu = fits.PrimaryHDU(data=arr, header=header)
    hdu.writeto(output_fits_image)
    print('Output to "{}"'.format(output_fits_image))




@click.command()
@click.argument('fits_image', type=click.Path(exists=True))
@click.option('--check-rate-image/--no-check-rate-image', default=True)
@click.option('--overwrite/--no-overwrite', default=True)
def main(fits_image, check_rate_image, overwrite):
    
    make_seed_image_for_rate_image(fits_image, check_rate_image, overwrite)



if __name__ == '__main__':
    main()

    
