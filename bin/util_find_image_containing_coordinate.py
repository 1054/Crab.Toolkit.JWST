#!/usr/bin/env python
# 
"""
Find the image(s) containing an input coordinate (RA Dec).

By Daizhong Liu @MPE. 

Last updates: 
    2023-06-03

"""
import os, sys, re, copy, json, shutil
import click
import itertools
import numpy as np
import multiprocessing as mp
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict
from matplotlib import cm
from matplotlib import colors as mpl_colors
from shapely.geometry import Point, Polygon

#from jwst.associations.asn_from_list import asn_from_list
#from jwst.associations import load_asn

import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

# code name and version
CODE_NAME = 'util_find_image_containing_coordinate.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230603'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)






def check_one_image_containing_coordinate(
        image_file, 
        scoord, 
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
    polygon = Polygon(np.column_stack((ra, dec)))
    return polygon.contains(Point(scoord.ra.deg, scoord.dec.deg))


def find_image_containing_coordinate(
        image_files, 
        scoord, 
        asn_file = None, 
    ):
    
    # check if an asn_file has been given
    if asn_file is not None:
        from jwst.associations import load_asn
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
        ret = pool.apply_async(check_one_image_containing_coordinate, 
                               args = (image_files[k],
                                       scoord,
                                      ),
        )
        rets.append(ret)
    pool.close()
    pool.join()
    matched_images = []
    for k, image_file in enumerate(image_files):
        result = rets[k].get()
        if result == True:
            matched_images.append(image_file)
    
    if len(matched_images) > 0:
        print('Found {} images containing the coordinate RA Dec {} {}:'.format(
            len(matched_images), scoord.ra.deg, scoord.dec.deg))
        for image_file in matched_images:
            print('  '+image_file)
    else:
        print('Found no image containing the coordinate RA Dec {} {}:'.format(
            scoord.ra.deg, scoord.dec.deg))
    
    # return
    return matched_images




@click.command()
@click.argument('image_files', nargs=-1, type=click.Path(exists=True))
#@click.argument('coordinate', nargs=2, type=str)
@click.option('--radec', 'coordinate', nargs=2, type=str)
@click.option('--asn-file', type=click.Path(exists=True), default=None, help='We can input an asn file instead of image files.')
def main(
        image_files, 
        coordinate, 
        asn_file, 
    ):

    if len(image_files) == 0 and asn_file is None:
        raise Exception('Please input image_files!')

    ra, dec = coordinate
    try: 
        float(ra)
        float(dec)
        scoord = SkyCoord(float(ra)*u.deg, float(dec)*u.deg)
    except:
        scoord = SkyCoord(ra + ' ' + dec, unit=(u.hour, u.deg))

    find_image_containing_coordinate(
        image_files, 
        scoord, 
        asn_file = asn_file, 
    )




if __name__ == '__main__':
    
    main()



