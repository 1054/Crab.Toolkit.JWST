#!/usr/bin/env python
# 
"""
Comparing images from different data reduction by drawing random circles and analyzing the statistics.

Usage: 
    ./util_compare_image_reduction.py image_1.fits image_2.fits image_3.fits output_name

Input: 
    FITS images

Output:
    Analysis files: 
        {output_name}_flux_vs_flux.pdf
        {output_name}_pos_vs_pos.pdf

By Daizhong Liu @MPE. 

Last updates: 
    2022-12-14 

"""
import os, sys, re, shutil
import astropy.units as u
import astropy.constants as const
import click
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, proj_plane_pixel_area
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300

# code name and version
CODE_NAME = 'util_compare_image_reduction.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221214'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)


# defaults
default_aperture_size = 1.5 # arcsec
default_aperture_number = 1000 


# define class
class ImageComparer(object):
    """docstring for ImageComparer"""
    def __init__(self, input_images):
        self.image_files = input_images
        self.hdulists = []
        self.imhdus = []
        self.extnames = []
        self.wcses = []
        self.pixscs = []
        self.rects = []
        for image_file in self.image_files:
            hdul = fits.open(image_file)
            if 'SCI' in hdul:
                extname = 'SCI'
            else:
                extname = ''
                for ihdu in range(len(hdul)):
                    if hdul[ihdu].header['NAXIS'] == 2:
                        extname = hdul[ihdu].name
                        break
                if extname == '':
                    raise Exception('Error! Could not find 2D image in fits files: {!r}'.format(image_file))
            logger.info('Reading {!r} extension {!r}'.format(image_file, extname))
            imhdu = hdul[extname]
            header = hdul[extname].header
            wcs = WCS(header, naxis=2)
            pixsc = np.sqrt(proj_plane_pixel_area(wcs))*3600.
            ny, nx = header['NAXIS2'], header['NAXIS1']
            rect = wcs.wcs_pix2world([1, nx, nx, 1], [1, 1, ny, ny], 1)
            self.hdulists.append(hdul)
            self.imhdus.append(imhdu)
            self.extnames.append(extname)
            self.wcses.append(wcs)
            self.pixscs.append(pixsc)
            self.rects.extend(np.array(rect).T.tolist()) # make rect = [ (x1, y1), (x2, y2), ... ]
        self.rects = np.array(self.rects)
    # 
    def __enter__(self):
        return self
    # 
    def __exit__(self, type, value, traceback):
        for hdul in self.hdulists:
            hdul.close()
    # 
    def draw_apertures(self, 
            output_name, 
            aperture_size = default_aperture_size, 
            aperture_number = default_aperture_number, 
        ): 
        apertures_dict = OrderedDict()
        apertures_dict['id'] = []
        apertures_dict['ra'] = []
        apertures_dict['dec'] = []
        for i in range(len(self.imhdus)):
            apertures_dict[f'x_{i}'] = []
            apertures_dict[f'y_{i}'] = []
            apertures_dict[f'rad_{i}'] = []
            apertures_dict[f'ix_{i}'] = []
            apertures_dict[f'iy_{i}'] = []
            apertures_dict[f'irad_{i}'] = []
            apertures_dict[f'fluxconv_{i}'] = []
            apertures_dict[f'mean_{i}'] = []
            apertures_dict[f'median_{i}'] = []
            apertures_dict[f'stddev_{i}'] = []
            apertures_dict[f'weighted_x_{i}'] = []
            apertures_dict[f'weighted_y_{i}'] = []
            apertures_dict[f'weighted_ra_{i}'] = []
            apertures_dict[f'weighted_dec_{i}'] = []
        chunk = 100
        count = 0
        # calc boundary
        max_ra = np.min([self.rects[0, 0], self.rects[3, 0]])
        min_ra = np.max([self.rects[1, 0], self.rects[2, 0]])
        min_dec = np.max([self.rects[0, 1], self.rects[1, 1]])
        max_dec = np.min([self.rects[2, 1], self.rects[3, 1]])
        while count < aperture_number:
            logger.info('Processing aperture {}'.format(count))
            rand_ra = np.random.random_sample()
            rand_dec = np.random.random_sample()
            ra = (max_ra - min_ra) * rand_ra + min_ra
            dec = (max_dec - min_dec) * rand_dec + min_dec
            ixyrad_list = []
            fxyrad_list = []
            for i in range(len(self.imhdus)):
                imhdu = self.imhdus[i]
                image = imhdu.data
                header = imhdu.header
                wcs = self.wcses[i]
                pixsc = self.pixscs[i]
                (fx,), (fy,) = wcs.wcs_world2pix([ra], [dec], 0)
                ix, iy = int(np.round(fx)), int(np.round(fy))
                ny, nx = header['NAXIS2'], header['NAXIS1']
                frad = aperture_size/2.0/pixsc
                irad = int(np.ceil(frad))
                if ix > irad and ix < nx-irad-1 and iy > irad and iy < ny-irad-1:
                    if np.all(np.isfinite(image[iy-irad:iy+irad+1, ix-irad:ix+irad+1])): 
                        if np.all(~np.isclose(image[iy-irad:iy+irad+1, ix-irad:ix+irad+1], 0)):
                            ixyrad_list.append([ix, iy, irad])
                            fxyrad_list.append([fx, fy, frad])
            # need all image to have valid pixels
            if len(ixyrad_list) != len(self.imhdus):
                continue
            # draw aperture
            apertures_dict['id'].append(count)
            apertures_dict['ra'].append(ra)
            apertures_dict['dec'].append(dec)
            for i in range(len(self.imhdus)):
                imhdu = self.imhdus[i]
                #image = imhdu.data
                header = imhdu.header
                wcs = self.wcses[i]
                pixsc = self.pixscs[i]
                ny, nx = header['NAXIS2'], header['NAXIS1']
                ix, iy, irad = ixyrad_list[i]
                fx, fy, frad = fxyrad_list[i]
                image = imhdu.data[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                #ggy, ggx = np.mgrid[0:ny, 0:nx]
                #gx = ggx[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                #gy = ggy[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                gy, gx = np.mgrid[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                mask = (np.sqrt((gx-fx)**2 + (gy-fy)**2) < frad)
                if header['BUNIT'] == 'MJy/sr':
                    fluxconv = 1
                elif header['BUNIT'] == '10.0*nanoJansky':
                    fluxconv = (10.0*u.nJy / (pixsc*u.arcsec)**2).to(u.MJy/u.sr).value
                else:
                    raise Exception('Error! BUNIT {!r} is not implemented!'.format(header['BUNIT']))
                aperdata = image[mask] * fluxconv
                fmask = mask.astype(int).astype(float)
                #print(gx.shape, image.shape, fmask.shape)
                weighted_x = np.nansum(gx*image*fmask) / np.nansum(image*fmask)
                weighted_y = np.nansum(gy*image*fmask) / np.nansum(image*fmask)
                (weighted_ra,), (weighted_dec,) = wcs.wcs_pix2world([weighted_x], [weighted_y], 0)
                apertures_dict[f'x_{i}'].append(fx)
                apertures_dict[f'y_{i}'].append(fy)
                apertures_dict[f'rad_{i}'].append(frad)
                apertures_dict[f'ix_{i}'].append(ix)
                apertures_dict[f'iy_{i}'].append(iy)
                apertures_dict[f'irad_{i}'].append(irad)
                apertures_dict[f'fluxconv_{i}'].append(fluxconv)
                apertures_dict[f'mean_{i}'].append(np.nanmean(aperdata))
                apertures_dict[f'median_{i}'].append(np.nanmedian(aperdata))
                apertures_dict[f'stddev_{i}'].append(np.nanstd(aperdata))
                apertures_dict[f'weighted_x_{i}'].append(weighted_x)
                apertures_dict[f'weighted_y_{i}'].append(weighted_y)
                apertures_dict[f'weighted_ra_{i}'].append(weighted_ra)
                apertures_dict[f'weighted_dec_{i}'].append(weighted_dec)
            # 
            # dump table
            if count > 0 and aperture_number > 10 and count % int(aperture_number/10) == 0 :
                apertures_table = Table(apertures_dict)
                apertures_table.write(f'{output_name}.dump.{count}.csv', overwrite=True)
                logger.info('Output to {!r}'.format(f'{output_name}.dump.{count}.csv'))
            # 
            # next count
            count += 1
        # 
        apertures_table = Table(apertures_dict)
        apertures_table.write(output_name+'.csv', overwrite=True)
        logger.info('Output to {!r}'.format(output_name+'.csv'))
        



# 
@click.command()
@click.argument('input_images', type=click.Path(exists=True), nargs=-1)
@click.argument('output_name', type=str, nargs=1)
@click.option('--aperture-size', type=float, default=default_aperture_size, help='aperture diameter in arcsec.')
@click.option('--aperture-number', type=int, default=default_aperture_number, help='aperture number.')
@click.option('--output-dir', type=click.Path(exists=False), default=None)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
@click.option('--verbose/--no-verbose', is_flag=True, default=True)
def main(
        input_images, 
        output_name, 
        aperture_size, 
        aperture_number,
        output_dir,
        overwrite,
        verbose, 
    ):
    
    with ImageComparer(input_images) as imcmp:
        imcmp.draw_apertures(output_name, aperture_size=aperture_size, aperture_number=aperture_number)
    



# 
if __name__ == '__main__':
    
    main()



