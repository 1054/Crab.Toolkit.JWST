#!/usr/bin/env python
# 
# A utility to make cutout image. 
# 
import os, sys, re, copy, datetime, shutil, time
import numpy as np
import click
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from reproject import reproject_interp


def make_cutout_using_reproject(
        input_image_file, 
        output_image_file, 
        FoV_RA, 
        FoV_Dec, 
        pixel_size = None, 
        center_RA = None, 
        center_Dec = None, 
        recenter_to_RA = None, 
        recenter_to_Dec = None, 
        set_output_bitpix_32 = False, 
    ):
    # 
    # read input fits image
    hdul = fits.open(input_image_file)
    ihdu = 0
    header = None
    while ihdu < len(hdul):
        try:
            ndim = len(hdul[ihdu].data.shape)
        except:
            ndim = 0
        if ndim >= 2:
            break
        ihdu += 1
    
    if ndim > 2:
        while ndim > 2:
            if hdul[ihdu].data.shape[0] == 1:
                hdul[ihdu].data = hdul[ihdu].data[0]
                if f'NAXIS{ndim}' in hdul[ihdu].header:
                    del hdul[ihdu].header[f'NAXIS{ndim}']
                hdul[ihdu].header['NAXIS'] -= 1
                ndim -= 1
            else:
                raise Exception('Error! The input data is not a 2D image nor could be reduced to a 2D image!')
    
    image = hdul[ihdu].data
    if ihdu == 0:
        header = hdul[ihdu].header
    else:
        header = copy.copy(hdul[0].header)
        header['COMMENT'] = ''
        header['COMMENT'] = 'Header for extension {} (1-{})'.format(ihdu, len(hdul))
        header['COMMENT'] = ''
        for key in hdul[ihdu].header:
            if key in header:
                try:
                    header[key+str(ihdu)] = hdul[ihdu].header[key]
                except:
                    pass
            else:
                header[key] = hdul[ihdu].header[key]
    hdul.close()
    # 
    # determine wcs, pixscale, x_size y_size
    wcs = WCS(header, naxis=2)
    wcs.sip = None
    if pixel_size is None:
        pixscales = proj_plane_pixel_scales(wcs) * 3600.0 # arcsec
        print('pixscales', pixscales, 'arcsec')
    else:
        pixscales = [-pixel_size, pixel_size]
    x_pixsc = np.abs(pixscales[0])
    y_pixsc = np.abs(pixscales[1])
    x_sizeF = FoV_RA / x_pixsc # arcsec
    y_sizeF = FoV_Dec / y_pixsc # arcsec
    x_size = int(np.ceil(FoV_RA / x_pixsc))
    y_size = int(np.ceil(FoV_Dec / y_pixsc))
    print('x_size', x_size, 'y_size', y_size)

    ny, nx = image.shape
    cra, cdec = wcs.wcs_pix2world((nx-1)/2.0, (ny-1)/2.0, 0)


    # convert RA Dec sexagesimal format
    if center_RA is None:
        center_RA = float(cra)
    if center_Dec is None:
        center_Dec = float(cdec)
    if isinstance(center_RA, str) or isinstance(center_Dec, str):
        center_skycoord = SkyCoord(str(center_RA), str(center_Dec), unit=(u.hour, u.deg))
        center_RA = center_skycoord.ra.deg
        center_Dec = center_skycoord.dec.deg


    # prepare cutout header
    #cutout_header = copy.deepcopy(header)
    cutout_header = fits.Header()
    cutout_header['BITPIX'] = -32
    cutout_header['NAXIS'] = 2
    cutout_header['NAXIS1'] = x_size
    cutout_header['NAXIS2'] = y_size
    cutout_header['CTYPE1'] = 'RA---TAN'
    cutout_header['CTYPE2'] = 'DEC--TAN'
    cutout_header['CUNIT1'] = 'deg'
    cutout_header['CUNIT2'] = 'deg'
    cutout_header['CDELT1'] = -x_pixsc / 3600.0
    cutout_header['CDELT2'] = y_pixsc / 3600.0
    cutout_header['CRPIX1'] = (x_size+1)/2. # 1-based number
    cutout_header['CRPIX2'] = (y_size+1)/2. # 1-based number
    cutout_header['CRVAL1'] = center_RA
    cutout_header['CRVAL2'] = center_Dec
    cutout_header['RADESYS'] = 'ICRS'
    cutout_header['EQUINOX'] = 2000
    #'NAXIS','NAXIS1','NAXIS2','CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2',
    for key in ['BUNIT','BMAJ','BMIN','BPA','TELESCOP','INSTRUME','FILTER','EXPTIME','PA_V3',
                'DATE-OBS','TIME-OBS','PHOTMODE','PHOTFLAM','PHTFLAM1','PHTFLAM2','PHTRATIO','PHOTFNU','PHOTZPT','PHOTPLAM','PHOTBW',
                'S_REGION',]:
        if key in header:
            cutout_header[key] = header[key]
    cutout_header['HISTORY'] = ''
    cutout_header['HISTORY'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ' + time.strftime('%Z')
    cutout_header['HISTORY'] = 'Created cutout image at center RA Dec %s %s with FoV size RA Dec %s %s [arcsec] from the input image file "%s"'%(\
                               center_RA, center_Dec, FoV_RA, FoV_Dec, input_image_file)
    cutout_header['HISTORY'] = ''
    #print(cutout_header)


    cutout_image, cutout_footprint = reproject_interp((image, wcs), cutout_header)


    # optionally recenter the image to some RA Dec for astrometry correction
    if recenter_to_RA is not None and recenter_to_Dec is not None:
        cutout_header['CRVAL1'] = recenter_to_RA
        cutout_header['CRVAL2'] = recenter_to_Dec


    # generate fits HDU
    output_hdu = fits.PrimaryHDU(data=cutout_image, header=cutout_header)

    if set_output_bitpix_32:
        output_hdu.header['BITPIX'] = -32
        output_hdu.data[np.isnan(output_hdu.data)] = 0.0
        output_hdu.data = output_hdu.data.astype(np.float32)


    # prepare to write output file
    output_file = output_image_file
    if output_file.find(os.sep) >= 0:
        if not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')

    output_hdu.writeto(output_file, overwrite=True)
    print('Written to "%s"'%(output_file))




@click.command()
@click.argument('input_fits_file', type=click.Path(exists=True))
@click.argument('output_fits_file', type=click.Path(exists=False), required=False, default=None)
@click.option('--fov', type=float, nargs=2, default=(None, None), help='Field of view in arcsec, two float values.')
@click.option('--pixsc', type=float, default=None, help='Pixel size in arcsec, one float value.')
@click.option('--center', type=float, nargs=2, default=(None, None), help='Center RA Dec, two float values.')
@click.option('--set-output-bitpix-32', type=bool, is_flag=True, default=False, help='Set output BITPIX to -32 for compatibility.')
def main(
        input_fits_file, 
        output_fits_file, 
        fov, 
        pixsc, 
        center, 
        set_output_bitpix_32, 
    ):
    
    make_cutout_using_reproject(
        input_fits_file, 
        output_fits_file, 
        fov[0], 
        fov[1], 
        pixel_size = pixsc, 
        center_RA = center[0], 
        center_Dec = center[1], 
        set_output_bitpix_32 = set_output_bitpix_32, 
    )



if __name__ == '__main__':
    
    main()



