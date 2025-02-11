#!/usr/bin/env python
# 
"""
Get footprints of images as DS9 region file.

By Daizhong Liu @MPE. 

Last updates: 
    2023-01-08
    2025-01-19 edge_detection

"""
import os, sys, re, copy, shutil
import click
import itertools
import numpy as np
import multiprocessing as mp
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict
from matplotlib import cm
from matplotlib import colors as mpl_colors
from shapely.geometry import Point, Polygon
from reproject import reproject_interp

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations import load_asn

import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

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
        n_footprints:list, 
        region_file:str = 'ds9.reg',
    ):
    with open(region_file, 'w') as fp:
        fp.write('# Region file format: DS9 version 4.1\n')
        fp.write('global color=green\n')
        fp.write('fk5\n')
        for k in range(len(footprints)):
            color = mpl_colors.to_hex(next(COLOR_CYCLER))
            n_footprint = n_footprints[k]
            footprint = footprints[k]
            if n_footprint == 1:
                vertices = footprint
                #vertices.extend([vertices[0], vertices[1]])
                fp.write('polygon({}) # text={{[{}]{}}} color={}\n'.format(
                        ', '.join(map(str, vertices)),
                        k, filenames[k], 
                        color
                    )
                )
            else:
                for i in range(n_footprint):
                    vertices = footprint[i]
                    #vertices.extend([vertices[0], vertices[1]])
                    fp.write('polygon({}) # text={{[{}-{}]{}}} color={}\n'.format(
                            ', '.join(map(str, vertices)),
                            k, i, filenames[k], 
                            color
                        )
                    )

def get_one_image_footprint(
        image_file,
        edge_detection,
        invalid_pixels,
        resample_pixel_size,
        edge_diff_min,
        edge_diff_max,
        output_masked_image,
        output_edge_image,
    ):
    header = None
    data = None
    with fits.open(image_file) as hdul:
        if 'SCI' in hdul:
            header = hdul['SCI'].header
            data = hdul['SCI'].data
        else:
            for hdu in hdul:
                if hdu.header['NAXIS'] >= 2:
                    header = hdu.header
                    data = hdu.data
                    break
        if header is None:
            raise Exception('Error! Could not read an image data from {!r}'.format(image_file))
    # 
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=Warning)
        wcs = WCS(header, naxis=2)
        wcs.sip = None
        pixscale = np.sqrt(proj_plane_pixel_area(wcs))
    nx, ny = header['NAXIS1'], header['NAXIS2']
    # 
    # the corners
    ra, dec = wcs.wcs_pix2world([1, nx, nx, 1, 1], [1, 1, ny, ny, 1], 1)
    footprint = np.column_stack((ra, dec)).ravel().tolist()
    n_footprint = 1
    # 
    if edge_detection:
        # mask out invalid pixels
        if invalid_pixels is None:
            invalid_pixels = []
        if len(invalid_pixels) > 0:
            data = data.astype(float)
            for invalid_pixel in invalid_pixels:
                print('Masking out invalid pixels with value {}'.format(invalid_pixel))
                data[data==invalid_pixel] = np.nan
        print('Detecting image edge for precise footprint ... (very time consuming)')
        masked = np.zeros(data.shape)
        masked[np.isfinite(data)] = 1.0
        # 
        # resample grid to about 200x200
        header2 = wcs.to_header()
        if resample_pixel_size is None:
            if nx > ny:
                rescale = float(ny)/200
            else:
                rescale = float(nx)/200
        else:
            rescale = (resample_pixel_size/3600.0) / pixscale
        nx2 = int(np.ceil(nx/rescale))
        ny2 = int(np.ceil(ny/rescale))
        pixscale2 = pixscale * rescale
        print('Resampling to pixel size {:.5g} -> {:.5g}, image size {}x{} -> {}x{}'.format(pixscale*3600.0, pixscale2*3600.0, nx, ny, nx2, ny2))
        header2['NAXIS'] = 2
        header2['NAXIS1'] = nx2
        header2['NAXIS2'] = ny2
        header2['CRPIX1'] = (nx2+1.0)/2.0
        header2['CRPIX2'] = (ny2+1.0)/2.0
        header2['CRVAL1'] = (ra[0]+ra[1])/2.0
        header2['CRVAL2'] = (dec[0]+dec[3])/2.0
        header2['CDELT1'] = -pixscale2
        header2['CDELT2'] = pixscale2
        header2['PC1_1'] = 1.0
        header2['PC1_2'] = 0.0
        header2['PC2_1'] = 0.0
        header2['PC2_2'] = 1.0
        masked2 = reproject_interp((masked, wcs), header2, return_footprint=False)
        masked2[masked2<0.5] = 0.0
        #fits.PrimaryHDU(data=masked2, header=header2).writeto('tmp_masked2.fits', overwrite=True)
        #print('Output to {!r}'.format('tmp_masked2.fits'))
        if output_masked_image is not None:
            check_dir_path = os.path.dirname(output_masked_image)
            if check_dir_path != '':
                if not os.path.isdir(check_dir_path):
                    os.makedirs(check_dir_path)
            fits.writeto(output_masked_image, masked2, header2, overwrite=True)
            print('Output to {!r}'.format(output_masked_image))
        # 
        # smooth and detect edge
        smoothed = convolve(masked2, Tophat2DKernel(3))
        #smoothed[smoothed < 0.95] = 0.0
        #smoothed[smoothed > 1.05] = 1.0
        diff2 = (smoothed-masked2)
        diff2[diff2<0] = 0.0
        #fits.PrimaryHDU(data=diff2, header=header2).writeto('tmp_diff2.fits', overwrite=True)
        #print('Output to {!r}'.format('tmp_diff2.fits'))
        if output_edge_image is not None:
            check_dir_path = os.path.dirname(output_edge_image)
            if check_dir_path != '':
                if not os.path.isdir(check_dir_path):
                    os.makedirs(check_dir_path)
            fits.writeto(output_edge_image, diff2, header2, overwrite=True)
            print('Output to {!r}'.format(output_edge_image))
        # 
        gy, gx = np.mgrid[0:ny2, 0:nx2]
        edge = np.logical_and(diff2 > edge_diff_min, diff2 < edge_diff_max)
        x = gx[edge].ravel().tolist()
        y = gy[edge].ravel().tolist()
        xsorted = [[]]
        ysorted = [[]]
        i = 0
        k = 0
        n = len(x)
        print('n {}'.format(n))
        while len(x) > 0:
            if i == 0:
                imin = 0
                dmin = 99.
            else:
                distances = ((xsorted[k][-1]-np.array(x))**2 + (ysorted[k][-1]-np.array(y))**2)
                imin = np.argmin(distances)
                dmin = distances[imin]
                # if i > 0.5*n:
                #     distances2 = ((np.array(xsorted)-x[imin])**2 + (np.array(ysorted)-y[imin])**2)
                #     if dmin > np.min(distances2):
                #         break
            x1 = x.pop(imin)
            y1 = y.pop(imin)
            if dmin < 2.0:
                continue
            if dmin > 100.0: # at most 10 pixels distance
                k += 1
                xsorted.append([])
                ysorted.append([])
            xsorted[k].append(x1)
            ysorted[k].append(y1)
            i += 1
        wcs2 = WCS(header2, naxis=2)
        footprint = []
        for k in range(len(xsorted)):
            if len(xsorted[k]) >= 10: # at least 10 points for a polyon
                ra, dec = wcs2.wcs_pix2world(xsorted[k], ysorted[k], 0)
                footprint.append(np.column_stack((ra, dec)).ravel().tolist())
        for k, v in enumerate(footprint):
            print('footprint {} has {} vertices'.format(k, len(v)))
        n_footprint = len(footprint)
        if n_footprint == 1:
            footprint = footprint[0]
    # 
    return footprint, n_footprint


def get_image_footprint_as_ds9_region(
        image_files, 
        output_region_file, 
        asn_file = None, 
        edge_detection = False, 
        invalid_pixels = None, 
        resample_pixel_size = None, 
        edge_diff_min = 0.1, 
        edge_diff_max = 0.9, 
        output_masked_image = None, 
        output_edge_image = None, 
        n_parallel = None, 
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
    if n_parallel is None:
        n_parallel = 4 # mp.cpu_count() // 4
    pool = mp.Pool(n_parallel)
    k = 0
    rets = []
    for k in range(len(image_files)):
        ret = pool.apply_async(get_one_image_footprint, 
                               args = (image_files[k],
                                       edge_detection,
                                       invalid_pixels,
                                       resample_pixel_size, 
                                       edge_diff_min, 
                                       edge_diff_max, 
                                       output_masked_image, 
                                       output_edge_image, 
                                      ),
        )
        rets.append(ret)
    
    pool.close()
    #pool.join()

    footprints = []
    n_footprints = []
    for k in range(len(image_files)):
        footprint, n_footprint = rets[k].get()
        footprints.append(footprint)
        n_footprints.append(n_footprint)
    
    # write ds9 region file
    write_ds9_region_file(
        image_files,
        footprints,
        n_footprints, 
        output_region_file,
    )
    print('Output to {!r}'.format(output_region_file))
    
    # return
    return




@click.command()
@click.argument('image_files', nargs=-1, required=True, type=click.Path(exists=True))
@click.argument('output_region_file', nargs=1, type=click.Path(exists=False))
@click.option('--asn-file', type=click.Path(exists=True), default=None, help='We can input an asn file instead of image files.')
@click.option('--edge-detection', type=bool, is_flag=True, default=False, help='Use slow edge detection to get the precise footprint. Default is False.')
@click.option('--invalid-pixels', type=float, multiple=True, default=[], help='Invalid pixel values. Default is just NaN pixels.')
@click.option('--resample-pixel-size', type=float, default=None, help='Resampling at this pixel size for faster processing. In arcsec unit. In default the resampling is done at an image size of 300x300.')
@click.option('--edge-diff-min', type=float, default=0.1, help='Smooth - original, value 0.0 - 1.0.')
@click.option('--edge-diff-max', type=float, default=0.9, help='Smooth - original, value 0.0 - 1.0.')
@click.option('--output-masked-image', type=click.Path(exists=False), default=None, help='Save masked image as FITS file.')
@click.option('--output-edge-image', type=click.Path(exists=False), default=None, help='Save edge image as FITS file.')
@click.option('--ncpu', 'n_parallel', type=int, default=None, help='Parallel processing.')
def main(
        image_files, 
        output_region_file, 
        asn_file, 
        edge_detection, 
        invalid_pixels, 
        resample_pixel_size, 
        edge_diff_min, 
        edge_diff_max, 
        output_masked_image, 
        output_edge_image, 
        n_parallel, 
    ):

    if len(image_files) == 0:
        raise Exception('Please input image_files!')

    get_image_footprint_as_ds9_region(
        image_files, 
        output_region_file, 
        asn_file = asn_file, 
        edge_detection = edge_detection, 
        invalid_pixels = invalid_pixels, 
        resample_pixel_size = resample_pixel_size, 
        edge_diff_min = edge_diff_min, 
        edge_diff_max = edge_diff_max, 
        output_masked_image = output_masked_image, 
        output_edge_image = output_edge_image, 
        n_parallel = n_parallel, 
    )




if __name__ == '__main__':
    
    main()



