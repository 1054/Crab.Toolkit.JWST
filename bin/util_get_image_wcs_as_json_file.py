#!/usr/bin/env python
# 
"""
Get WCS keywords of images as json file.

By Daizhong Liu @MPE. 

Last updates: 
    2023-06-03

"""
import os, sys, re, copy, json, shutil
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

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations import load_asn

import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

# code name and version
CODE_NAME = 'util_get_image_wcs_as_json_file.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230603'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)






def get_one_image_wcs(
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
    footprint = np.column_stack((ra, dec)).ravel()
    wcs_header = wcs.to_header()
    wcs_dict = {}
    wcs_dict['FILENAME'] = image_file
    for key in wcs_header:
        wcs_dict[key] = wcs_header[key]
    wcs_dict['S_REGION'] = 'poylgon ' + ' '.join(['{:.8f}'.format(t) for t in footprint])
    return wcs_dict


def get_image_wcs_as_json(
        image_files, 
        output_json_file, 
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
        ret = pool.apply_async(get_one_image_wcs, 
                               args = (image_files[k],
                                      ),
        )
        rets.append(ret)
    pool.close()
    pool.join()
    wcs_dicts = []
    for k in range(len(image_files)):
        wcs_dict = rets[k].get()
        wcs_dicts.append(wcs_dict)
    
    # write json file
    with open(output_json_file, 'w') as fp:
        if len(wcs_dicts) == 1:
            json.dump(wcs_dicts[0], fp, indent=4)
        else:
            json.dump(wcs_dicts, fp, indent=4)
    
    print('Output to {!r}'.format(output_json_file))
    
    # return
    return




@click.command()
@click.argument('image_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.argument('output_json_file', nargs=1, type=click.Path(exists=False))
@click.option('--asn-file', type=click.Path(exists=True), default=None, help='We can input an asn file instead of image files.')
def main(
        image_files, 
        output_json_file, 
        asn_file, 
    ):

    if len(image_files) == 0:
        raise Exception('Please input image_files!')

    get_image_wcs_as_json(
        image_files, 
        output_json_file, 
        asn_file = asn_file, 
    )




if __name__ == '__main__':
    
    main()



