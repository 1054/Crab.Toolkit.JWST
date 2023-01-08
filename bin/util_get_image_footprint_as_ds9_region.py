#!/usr/bin/env python
# 
"""
Get footprints of images as DS9 region file.

By Daizhong Liu @MPE. 

Last updates: 
    2023-01-08

"""
import os, sys, re, copy, shutil
import click
import itertools
import numpy as np
import multiprocessing as mp
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict
from matplotlib import cm
from matplotlib import colors as mpl_colors
from shapely.geometry import Point, Polygon

from jwst.associations import load_asn

# code name and version
CODE_NAME = 'util_get_image_footprint_as_ds9_region.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230108'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)



COLOR_CYCLER = itertools.cycle(cm.rainbow(np.linspace(0, 1, 7)))


def write_ds9_region_file(
        filenames:list, 
        footprints:list, 
        region_file:str = 'ds9.reg',
    ):
    with open(region_file, 'w') as fp:
        fp.write('# Region file format: DS9 version 4.1\n')
        fp.write('global color=green\n')
        fp.write('fk5\n')
        for i in range(len(footprints)):
            color = mpl_colors.to_hex(next(COLOR_CYCLER))
            vertices = footprints[i]
            #vertices.extend([vertices[0], vertices[1]])
            fp.write('polygon({}) # text={{[{}]{}}} color={}\n'.format(
                    ', '.join(map(str, vertices)),
                    i, filenames[i], 
                    color
                )
            )

def get_one_image_footprint(
        image_file,
    ):
    with fits.open(image_file) as hdul:
        header = None
        if 'SCI' in hdul:
            header = hdul['SCI'].header
        else:
            for hdu in hdul:
                if hdu.header['NAXIS'] >= 2:
                    header = hdu.header
                    break
        if header is None:
            raise Exception('Error! Could not read an image data from {!r}'.format(image_file))
    wcs = WCS(header, naxis=2)
    wcs.sip = None
    nx, ny = header['NAXIS1'], header['NAXIS2']
    ra, dec = wcs.wcs_pix2world([1, nx, nx, 1, 1], [1, 1, ny, ny, 1], 1)
    footprint = np.column_stack((ra, dec)).ravel().tolist()
    return footprint


def get_image_footprint_as_ds9_region(
        image_files, 
        output_region_file, 
        asn_file = None, 
    ):
    
    # check if an asn_file has been given
    if asn_file is not None:
        with open(asn_file, 'r') as fp:
            asn_dict = load_asn(fp)
            image_files = [t['expname'] for t in asn_dict['products'][0]['members'] if t['exptype']=='science']
        print('Processing {} images in the asn file {!r}'.format(len(image_files), asn_file))
    else:
        print('Processing {} images from the command line input'.format(len(image_files)))
    
    # run in parallel
    n_parallel = mp.cpu_count()
    pool = mp.Pool(n_parallel)
    k = 0
    rets = []
    for k in range(len(image_files)):
        ret = pool.apply_async(get_one_image_footprint, 
                               args = (image_files[k],
                                      ),
        )
        rets.append(ret)
    pool.close()
    pool.join()
    footprints = []
    for k in range(len(image_files)):
        footprint = rets[k].get()
        footprints.append(footprint)
    
    # write ds9 region file
    write_ds9_region_file(
        image_files,
        footprints,
        output_region_file,
    )
    print('Output to {!r}'.format(region_file))
    
    # return
    return




@click.command()
@click.argument('image_files', nargs=-1, type=click.Path(exists=True))
@click.argument('output_region_file', nargs=1, type=click.Path(exists=False))
@click.option('--asn-file', type=click.Path(exists=True), default=None, help='We can input an asn file instead of image files.')
def main(
        image_files, 
        output_region_file, 
        asn_file, 
    ):

    get_image_footprint_as_ds9_region(
        image_files, 
        output_region_file, 
        asn_file = asn_file, 
    )




if __name__ == '__main__':
    
    main()



