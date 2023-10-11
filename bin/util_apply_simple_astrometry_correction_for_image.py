#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os, sys, re, shutil
import glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, proj_plane_pixel_area
#from astropy.convolution import convolve, Gaussian2DKernel, Box1DKernel
import astropy.units as u
import click


# code name and version
CODE_NAME = 'util_apply_simple_astrometry_correction_for_catalog.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20231010'
CODE_HOMEPAGE = ''


# main
@click.command()
@click.argument('input_image_file', type=click.Path(exists=True))
@click.argument('output_image_file', type=click.Path(exists=False))
@click.option('--old-ra-dec', nargs=2, type=str, required=True, help='Input the old RA Dec for an anchor point, in degrees.')
@click.option('--new-ra-dec', nargs=2, type=str, required=True, help='Input the new RA Dec for an anchor point, in degrees.')
def main(
        input_image_file, 
        output_image_file, 
        old_ra_dec, 
        new_ra_dec,
    ):
    """
    This script will update the CRPIX1 CRPIX2 to achieve the simple astrometry correction. 
    
    The image needs to have the north in the up direction. 
    
    """
    
    print('Reading image %r'%(input_image_file))
    image, header = fits.getdata(input_image_file, header=True)
    wcs = WCS(header, naxis=2)
    wcs.sip = None
    
    old_ra, old_dec = old_ra_dec
    new_ra, new_dec = new_ra_dec
    try:
        old_ra_dec_scoord = SkyCoord(float(old_ra)*u.deg, float(old_dec)*u.deg)
    except:
        old_ra_dec_scoord = SkyCoord(old_ra, old_dec, unit=(u.hour, u.deg))
    try:
        new_ra_dec_scoord = SkyCoord(float(new_ra)*u.deg, float(new_dec)*u.deg)
    except:
        new_ra_dec_scoord = SkyCoord(new_ra, new_dec, unit=(u.hour, u.deg))
    
    (old_x,), (old_y,) = wcs.wcs_world2pix([old_ra_dec_scoord.ra.deg,], [old_ra_dec_scoord.dec.deg,], 1)
    (new_x,), (new_y,) = wcs.wcs_world2pix([new_ra_dec_scoord.ra.deg,], [new_ra_dec_scoord.dec.deg,], 1)
    
    d_x = old_x - new_x
    d_y = old_y - new_y
    header['CRPIX1'] = header['CRPIX1'] + d_x
    header['CRPIX2'] = header['CRPIX2'] + d_y
    header['HISTORY'] = 'Applied astrometry correction (CRPIX {:+.3f} {:+.3f}) with the anchor point at old RA Dec {} {} and new RA Dec {} {} for the input image "{}".'.format(d_x, d_y, old_ra, old_dec, new_ra, new_dec, input_image_file)

    hdu = fits.PrimaryHDU(data=image, header=header)

    if os.path.isfile(output_image_file):
        print('Found existing output file, backing up as %r'%(output_image_file+'.backup'))
        shutil.move(output_image_file, output_image_file+'.backup')
    if output_image_file.find(os.sep)>=0:
        if not os.path.isdir(os.path.dirname(output_image_file)):
            os.makedirs(os.path.dirname(output_image_file))
    print('Writing image %r'%(output_image_file))
    hdu.writeto(output_image_file)



########
# MAIN #
########

if __name__ == '__main__':
    
    main()



