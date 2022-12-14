#!/usr/bin/env python
# 
"""
Cut an input catalog to image field of view.

Usage: 
    ./util_cut_catalog_to_image_field_of_view.py some_catalog.fits some_image.fits output_subcatalog.fits

By Daizhong Liu @MPE. 

Last updates: 
    2022-11-30 

"""
import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area


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





def cut_catalog_to_image_field_of_view(
        input_catalog_file, 
        fits_image_file, 
        output_catalog_file, 
        buffer = 1.0, # arcsec at each side
        save_region_file = True, 
        meta_name = None, 
        overwrite = False, 
        verbose = True, 
    ):
    
    if verbose:
        logger.info('Cutting catalog {!r} to image {!r} field of view'.format(input_catalog_file, fits_image_file))
    
    # check output file
    if os.path.isfile(output_catalog_file) and not overwrite:
        logger.info('Found existing output catalog file {!r} and overwrite is False. Skipping.'.format(output_catalog_file))
        return
    
    # read input file
    catalog = Table.read(input_catalog_file)
    
    sci_data = None
    rms_data = None
    wht_data = None
    sci_header = None
    rms_header = None
    wht_header = None
    main_header = None
    with fits.open(fits_image_file) as hdul:
        main_header = copy.deepcopy(hdul[0].header)
        for hdu in hdul:
            # read first valid hdu as sci data
            if sci_data is None and hdu.data is not None and len(hdu.data.shape) >= 2:
                sci_data = copy.copy(hdu.data)
                sci_header = copy.deepcopy(hdu.header)
            # if the hdu has an extension name SCI, use it as the sci data
            if 'EXTNAME' in hdu.header:
                if hdu.header['EXTNAME'] == 'SCI':
                    sci_data = copy.copy(hdu.data)
                    sci_header = copy.deepcopy(hdu.header)
                elif hdu.header['EXTNAME'] == 'ERR':
                    rms_data = copy.copy(hdu.data)
                    rms_header = copy.deepcopy(hdu.header)
                elif hdu.header['EXTNAME'] == 'WHT':
                    wht_data = copy.copy(hdu.data)
                    wht_header = copy.deepcopy(hdu.header)
    
    # 
    if sci_data is None:
        raise Exception('Could not read SCI data from the input fits file {!r}'.format(fits_image_file))
    
    # mask out zero-weight pixels
    if wht_data is not None:
        mask_zero_weight = (wht_data==0)
        if np.count_nonzero(mask_zero_weight) > 0:
            sci_data[mask_zero_weight] = np.nan
            rms_data[mask_zero_weight] = np.nan
    
    # get wcs
    wcs = WCS(sci_header, naxis=2)
    
    # get RA Dec
    colname_RA = None
    for possible_colname in ['RA', 'ra', 'ALPHA_J2000']:
        try:
            idx = catalog.colnames.index(possible_colname)
        except:
            idx = -1
        if idx >= 0:
            colname_RA = catalog.colnames[idx]
            break
    if colname_RA is None:
        raise Exception('Could not find RA column in the input catalog {!r} (colnames: {})'.format(input_catalog_file, repr(catalog.colnames)))
    
    colname_DEC = None
    for possible_colname in ['DEC', 'Dec', 'dec', 'DELTA_J2000']:
        try:
            idx = catalog.colnames.index(possible_colname)
        except:
            idx = -1
        if idx >= 0:
            colname_DEC = catalog.colnames[idx]
            break
    if colname_DEC is None:
        raise Exception('Could not find DEC column in the input catalog {!r} (colnames: {})'.format(input_catalog_file, repr(catalog.colnames)))
    
    # convert input catalog RA DEC to pixel coordinates
    ra, dec = catalog[colname_RA].data, catalog[colname_DEC].data # must be in degrees
    px, py = wcs.wcs_world2pix(ra, dec, 0)
    pixsc = np.sqrt(proj_plane_pixel_area(wcs))*3600.0
    nx, ny = sci_header['NAXIS1'], sci_header['NAXIS2']
    buf = np.abs(buffer / pixsc) # 1.0 arcsec buffer
    select = np.logical_and.reduce((px > -buf, px < nx-1+buf,
                                    py > -buf, py < ny-1+buf))
    
    # select sources
    if verbose:
        logger.info('Selecting {} sources'.format(np.count_nonzero(select)))
    output_catalog = catalog[select]
    if colname_RA != 'RA': 
        output_catalog.rename_column(colname_RA, 'RA')
    if colname_DEC != 'DEC': 
        output_catalog.rename_column(colname_DEC, 'DEC')
    
    # add meta
    if meta_name is None:
        meta_name = os.path.splitext(os.path.basename(input_catalog_file))[0]
    if meta_name != '':
        output_catalog.meta['name'] = meta_name
    
    # save output catalog
    if os.path.isfile(output_catalog_file):
        shutil.move(output_catalog_file, output_catalog_file+'.backup')
    else:
        output_dir = os.path.dirname(output_catalog_file)
        if output_dir != '' and not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    output_catalog.write(output_catalog_file)
    if verbose:
        logger.info('Output to "{}"'.format(output_catalog_file))
    
    # save region file
    if save_region_file:
        output_region_file = os.path.splitext(output_catalog_file)[0] + '.ds9.reg'
        if os.path.isfile(output_region_file):
            shutil.move(output_region_file, output_region_file+'.backup')
        with open(output_region_file, 'w') as fp:
            fp.write('# DS9 Region file\n')
            fp.write('fk5\n')
            for i in range(len(output_catalog)):
                fp.write('circle({:.8f},{:.8f},{:.3f}")\n'.format(
                    output_catalog['RA'][i], 
                    output_catalog['DEC'][i], 
                    0.5, # default radius 0.5 arcsec
                ))
        if verbose:
            logger.info('Output to "{}"'.format(output_region_file))




@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('fits_image_file', type=click.Path(exists=True))
@click.argument('output_catalog_file', type=click.Path(exists=False))
@click.option('--buffer', type=float, default=1.0, help='Buffer in arcsec.')
@click.option('--meta-name', type=str, default=None, help='Meta name to be written into the output catalog fits header.')
@click.option('--save-region-file/--no-save-region-file', is_flag=True, default=True)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
@click.option('--verbose/--no-verbose', is_flag=True, default=True)
def main(
        input_catalog_file, 
        fits_image_file, 
        output_catalog_file, 
        buffer, 
        meta_name, 
        save_region_file, 
        overwrite, 
        verbose, 
    ):
    
    cut_catalog_to_image_field_of_view(
        input_catalog_file, 
        fits_image_file, 
        output_catalog_file, 
        buffer = buffer, 
        meta_name = meta_name, 
        save_region_file = save_region_file, 
        overwrite = overwrite, 
        verbose = verbose, 
    )




if __name__ == '__main__':
    
    main()



