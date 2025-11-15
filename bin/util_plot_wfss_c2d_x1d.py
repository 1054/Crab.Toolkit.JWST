#!/usr/bin/env python
#
# Find SLTNAME 
# 
import os, sys, re, json, copy, glob, shutil
import click
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from pprint import pprint
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt



@click.command()
@click.argument('c2d_file', type=click.Path(exists=True))
@click.argument('x1d_file', type=click.Path(exists=True), required=False, default=None)
@click.option('--lsigma', type=float, default=-5.0, help='lower N sigma for plotting, default is -5.0')
@click.option('--usigma', type=float, default=10.0, help='upper N sigma for plotting, default is +10.0')
@click.option('--xflip', type=bool, default=None, help='flip image x axis? default is True')
@click.option('--makefig', type=bool, is_flag=True, default=False, help='make output figure?')
@click.option('--figfile', type=click.Path(), default=None, help='output name for saving the figure')
def main(c2d_file, x1d_file, lsigma, usigma, xflip, makefig, figfile):
    if x1d_file is None:
        x1d_file = re.sub(r'(_cal|_s2d)\.fits$', r'', c2d_file) + '_x1d.fits'
        if not os.path.exists(x1d_file):
            raise Exception("Error! Please also input x1d file!")
    print('Input c2d/cal file: {}'.format(c2d_file))
    print('Input c1d/x1d file: {}'.format(x1d_file))
    image_data, image_header = fits.getdata(c2d_file, header=True)
    image_header0 = fits.getheader(c2d_file, 0)
    spec_data, spec_header = fits.getdata(x1d_file, header=True)
    wcs = WCS(image_header, naxis=2)
    if image_header['CTYPE1'] == 'WAVELENGTH':
        pixsc = image_header['CDELT2'] * 3600.0
        if xflip is None:
            xflip = False
    else:
        pixsc = np.sqrt(proj_plane_pixel_area(wcs)) * 3600.0
        if xflip is None:
            xflip = True
    extent = np.repeat(np.array(image_data.data.shape[::-1]) * pixsc, 2) * np.array([-1, 1, -1, 1])
    print('pixsc', pixsc, 'extent', extent)
    if xflip:
        image_data = image_data[:, ::-1]
        extent[0] *= -1
        extent[1] *= -1
    wave = spec_data['WAVELENGTH']
    flux = spec_data['FLUX']
    print('spec_data', type(spec_data), spec_data._coldefs.names)
    try:
        fluxerr = spec_data['FLUX_ERROR']
    except:
        fluxerr = spec_data['FLUX_ERR']
    fluxunit = spec_header['TUNIT2']

    wmin = np.nanmin(wave)
    wmax = np.nanmax(wave)
    wave1 = np.nanpercentile(wave, 15.0)
    wave2 = np.nanpercentile(wave, 85.0)
    wave_mask = np.logical_and(wave>=wave1, wave<=wave2)

    image_mask = np.repeat(wave_mask[np.newaxis, :], image_data.shape[0], axis=0)
    mean, med, sig = sigma_clipped_stats(image_data[image_mask])
    vmin = max(med+lsigma*sig, np.nanmin(image_data[image_mask]))
    vmax = min(mean+usigma*sig, np.nanmax(image_data[image_mask]))

    mean, med, sig = sigma_clipped_stats(flux[wave_mask])
    fmin = max(med+lsigma*sig, np.nanmin(flux[wave_mask]))
    fmax = min(mean+usigma*sig, np.nanmax(flux[wave_mask]))

    fluxbase = np.zeros(flux.shape)
    signal_mask = np.logical_and(wave_mask, np.abs(flux-med) > 3.0*sig)
    flux_masked = copy.copy(flux)
    signal_mask = (convolve(signal_mask.astype(int).astype(float), Gaussian1DKernel(4.0)) > 0.1)
    flux_smoothed = convolve(flux_masked, Box1DKernel(10), nan_treatment='fill', )
    fluxbase = flux_smoothed

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 5))

    if 'SLTNAME' in image_header:
        sid = image_header['SLTNAME']
        ra = image_header['SRCRA']
        dec = image_header['SRCDEC']
        coord = SkyCoord(ra*u.deg, dec*u.deg)
        radecstr = coord.to_string('hmsdms', sep=':', precision=3)
        pa = image_header['PA_APER']
        filt = image_header0['FILTER']
        pupil = image_header0['PUPIL']
        fig.suptitle(f'SLIT ID {sid} - {radecstr} - PA: {pa:.1f} - {filt} - {pupil}')
    elif 'TARGET' in image_header0 and 'RA' in image_header0 and 'DEC' in image_header0 and 'GRATING' in image_header0:
        sid = image_header0['TARGET']
        ra = image_header0['RA']
        dec = image_header0['DEC']
        coord = SkyCoord(ra*u.deg, dec*u.deg)
        radecstr = coord.to_string('hmsdms', sep=':', precision=3)
        filt = image_header0['FILTER']
        grat = image_header0['GRATING']
        fig.suptitle(f'SLIT ID {sid} - {radecstr} - {filt} - {grat}')

    ax1.imshow(image_data, origin='lower', interpolation='nearest', 
               extent=extent, vmin=vmin, vmax=vmax)
    ax1.set_aspect((image_data.shape[1]/image_data.shape[0])/7.0)
    ax1.set_ylabel('Offset [arcsec]')

    ax2.plot(wave, flux, alpha=0.8, zorder=30)
    ax2.fill_between(wave, fluxbase-fluxerr, fluxbase+fluxerr, 
                     color='k', alpha=0.2, lw=0, zorder=40)
    ax2.set_ylim(fmin, fmax)
    ax2.set_xlim(wmin, wmax)
    ax2.xaxis.set_major_locator(mpl.ticker.MaxNLocator(12))
    ax2.set_xlabel('Wavelength [micron]')
    ax2.set_ylabel(f'Flux [{fluxunit}]')
    
    fig.tight_layout()

    if makefig and figfile is None:
        figfile = c2d_file.replace('_cal.fits', '') + '_plot.png'

    if figfile is not None:
        figdir = os.path.dirname(figfile)
        if figdir != '' and figdir != '.' and not os.path.exists(figdir):
            os.makedirs(figdir)
        if os.path.exists(figfile):
            shutil.move(figfile, figfile+'.backup')
        fig.savefig(figfile, transparent=False, dpi=300)
        print('Output to {!r}'.format(figfile))
    
    plt.show(block=True)


if __name__ == '__main__':
    main()


