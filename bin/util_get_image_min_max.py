#!/usr/bin/env python
# 
"""
Get image min max values. 

By Daizhong Liu @MPE. 

Last updates: 
    2023-08-21

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
CODE_NAME = 'util_get_image_min_max.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230814'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)






@click.command()
@click.argument('image_file', type=click.Path(exists=True))
def main(
        image_file, 
    ):

    with fits.open(image_file) as hdul:
        iext = 0
        for hdu in hdul:
            if hdu.header['NAXIS'] >= 2:
                print('{} {}'.format(np.nanmin(hdu.data), np.nanmax(hdu.data)))
                break
            iext += 1




if __name__ == '__main__':
    
    main()



