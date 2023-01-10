#!/usr/bin/env python
# 
"""
Make groups of overlaping images based on their spatial coverages. 

Usage: 
    (Used in Python as a module)
    from util_run_mosaic_image_subgrouping import run_mosaic_image_subgrouping
    run_mosaic_image_subgrouping(image_files)
    
    (Run from the command line)
    ./util_run_mosaic_image_subgrouping.py input_images*.fits
    (or)
    ./util_run_mosaic_image_subgrouping.py --asn-file asn_file.json
    
Outputs:
    Output some directories:
        run_mosaic_image_subgroup_1
        run_mosaic_image_subgroup_2
        run_mosaic_image_subgroup_3
        ...

By Daizhong Liu @MPE. 

Last updates: 
    2023-01-09

"""
import os, sys, re, copy, json, shutil
import click
import itertools
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
import warnings
from astropy.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

from collections import OrderedDict
from matplotlib import cm
from matplotlib import colors as mpl_colors

from stpipe import Pipeline
#from jwst import datamodels
#from jwst.datamodels import ModelContainer 
#from jwst.outlier_detection.outlier_detection_step import OutlierDetectionStep
#from jwst.pipeline.calwebb_image3 import Image3Pipeline
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations import load_asn

import multiprocessing as mp

from shapely.geometry import Point, Polygon

# code name and version
CODE_NAME = 'util_run_mosaic_image_subgrouping.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230109'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)



COLOR_CYCLER = itertools.cycle(cm.rainbow(np.linspace(0, 1, 7)))


def write_ds9_region_file(
        group_table, 
        regionfile = 'ds9.reg',
    ):
    with open(regionfile, 'w') as fp:
        fp.write('# Region file format: DS9 version 4.1\n')
        fp.write('global color=green\n')
        fp.write('fk5\n')
        for i in range(len(group_table)):
            color = mpl_colors.to_hex(next(COLOR_CYCLER))
            vertices = [group_table['ra1'][i], group_table['dec1'][i], 
                        group_table['ra2'][i], group_table['dec1'][i], 
                        group_table['ra2'][i], group_table['dec2'][i], 
                        group_table['ra1'][i], group_table['dec2'][i]]
            fp.write('polygon({}) # text={{[{}]{}}} color={}\n'.format(
                    ','.join(map(str, vertices)),
                    i, group_table['group_dir'][i], 
                    color
                )
            )


class ImageFootPrint():
    """docstring for ImageFootPrint"""
    def __init__(self, ra_array, dec_array):
        self.ra_array = np.array(ra_array)
        self.dec_array = np.array(dec_array)
        self.polygon = Polygon(self.footprint_vertices())
    def footprint_vertices(self):
        ra_array = np.concatenate([self.ra_array, [self.ra_array[0]]]).ravel()
        dec_array = np.concatenate([self.dec_array, [self.dec_array[0]]]).ravel()
        footprint_vertices = np.column_stack((ra_array, dec_array))
        return footprint_vertices
    def min_ra(self):
        return np.min(self.ra_array)
    def min_dec(self):
        return np.min(self.dec_array)
    def max_ra(self):
        return np.max(self.ra_array)
    def max_dec(self):
        return np.max(self.dec_array)
    def intersects(self, shape_object):
        return self.polygon.intersects(shape_object)


def get_one_image_footprint(
        image_file,
    ):
    header = None
    with fits.open(image_file) as hdul:
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
    ra, dec = wcs.wcs_pix2world([1, nx, nx, 1], [1, 1, ny, ny], 1) # lower left, lower right, upper right, upper left
    footprint = ImageFootPrint(ra, dec)
    del wcs, header, ra, dec
    return footprint


def get_min_max_ra_dec_for_footprints(footprints):
    min_ra = +np.inf
    max_ra = -np.inf
    min_dec = +np.inf
    max_dec = -np.inf
    for k in range(len(footprints)):
        min_ra = min(min_ra, footprints[k].min_ra())
        max_ra = max(max_ra, footprints[k].max_ra())
        min_dec = min(min_dec, footprints[k].min_dec())
        max_dec = max(max_dec, footprints[k].max_dec())
    return min_ra, max_ra, min_dec, max_dec



def run_mosaic_image_subgrouping(
        image_files, 
        output_name = 'run_mosaic_image_subgroup', 
        asn_file = None, 
        cat_file = None, 
        grid_step = 10.0, # arcmin
        anchor_ra = None, 
        anchor_dec = None, 
        pixel_size = 1.0, # arcsec
        cpu_cores = 'all', 
        save_group_meta_file = None, # if None then '{output_name}_meta.txt'
        save_group_table_file = None, # if None then '{output_name}_table.txt'
        save_group_region_file = None, # if None then '{output_name}_regions.reg'
        overwrite = False, 
        verbose = True, 
    ):
    """
    Make groups of overlaping images based on their spatial coverages. 
    
    We will group the images in this way. 
    - First, we find the full RA Dec range of all images.
    - Then, 
    - Then, 
    
    Args:
        
        cpu_cores: 'all', 'none', or an integer
    
    """
    
    
    # check if an asn_file has been given
    if asn_file is not None:
        with open(asn_file, 'r') as fp:
            asn_dict = load_asn(fp)
            image_files = [t['expname'] for t in asn_dict['products'][0]['members'] if t['exptype']=='science']
        logger.info('Processing {} images in the asn file {!r}'.format(len(image_files), asn_file))
    else:
        logger.info('Processing {} images from the command line input'.format(len(image_files)))
    
    # check if image file names have duplicates (TODO)
    #check_duplicate_names = (len(list(set(list(map(os.path.basename, image_files))))) == len(image_files))
    #if check_duplicate_names:
    #    raise Exception('Error! The input image_files contain duplicated file names! We cannot process that for now.')
    
    # run in parallel
    if cpu_cores == 'all':
        n_parallel = mp.cpu_count()
    elif cpu_cores == 'none':
        n_parallel = 0
    elif isinstance(cpu_cores, int) or re.match(r'^[0-9]+$', cpu_cores):
        n_parallel = int(cpu_cores)
    else:
        n_parallel = mp.cpu_count()
    
    n_parallel = 0 #<TODO># multiprocessing/resource_tracker.py:216: UserWarning: resource_tracker: There appear to be 24 leaked semaphore objects to clean up at shutdown
    if n_parallel > 0:
        logger.info('Extracting footprints in {} parallel threads'.format(n_parallel))
        # 
        footprints = []
        with mp.Pool(n_parallel) as pool:
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
            # 
            for k in range(len(image_files)):
                footprint = rets[k].get()
                footprints.append(footprint)
    else:
        logger.info('Extracting footprints ...')
        # 
        footprints = []
        for k in range(len(image_files)):
            footprint = get_one_image_footprint(image_files[k])
            footprints.append(footprint)
    
    min_ra, max_ra, min_dec, max_dec = get_min_max_ra_dec_for_footprints(footprints)
    
    # define lower left (1,1) pixel as origin (1-base)
    origin_ra = max_ra
    origin_dec = min_dec
    
    # define step in integer pixels
    step_nx = int(np.ceil(grid_step * 60.0 / pixel_size))
    step_ny = int(np.ceil(grid_step * 60.0 / pixel_size))
    
    # define dec grid and anchor dec
    step_dec = float(step_ny) * pixel_size / 3600.0 # in degrees
    grid_dec = np.arange(min_dec, max_dec, step_dec) # increasing order, +Y
    if anchor_dec is None:
        #anchor_dec = (min_dec + max_dec) / 2.0
        anchor_dec = min_dec + 0.5 * float(len(grid_dec)) * step_dec
    
    # define ra grid and anchor ra
    step_ra = -float(step_nx) * pixel_size / 3600.0 / np.cos(np.deg2rad(anchor_dec)) # in degrees, negative
    grid_ra = np.arange(max_ra, min_ra, step_ra) # decreasing order, +X
    if anchor_ra is None:
        #anchor_ra = (min_ra + max_ra) / 2.0
        anchor_ra = max_ra + 0.5 * float(len(grid_ra)) * step_ra
    
    # get reference point pixel coordinate (1-base)
    mosaic_crpix1 = -(anchor_ra - origin_ra) * np.cos(np.deg2rad(anchor_dec)) * 3600.0 / pixel_size + 1
    mosaic_crpix2 = (anchor_dec - origin_dec) * 3600.0 / pixel_size + 1
    
    # 
    group_meta_dict = OrderedDict()
    group_meta_dict['n_groups'] = len(grid_ra)*len(grid_dec)
    group_meta_dict['n_col'] = len(grid_ra)
    group_meta_dict['n_row'] = len(grid_dec)
    group_meta_dict['naxis1'] = step_nx * len(grid_ra)
    group_meta_dict['naxis2'] = step_ny * len(grid_dec)
    group_meta_dict['crval1'] = anchor_ra
    group_meta_dict['crval2'] = anchor_dec
    group_meta_dict['crpix1'] = mosaic_crpix1
    group_meta_dict['crpix2'] = mosaic_crpix2
    group_meta_dict['cdelt1'] = -pixel_size / 3600.0 # negative
    group_meta_dict['cdelt2'] = pixel_size / 3600.0 # 
    group_meta_dict['crota2'] = 0.0
    group_meta_dict['ra1'] = origin_ra
    group_meta_dict['dec1'] = origin_dec
    group_meta_dict['ra2'] = origin_ra + step_ra
    group_meta_dict['dec2'] = origin_dec + step_dec
    group_meta_dict['step_nx'] = step_nx # pixels
    group_meta_dict['step_ny'] = step_ny # pixels
    group_meta_dict['image_width'] = step_nx * len(grid_ra) * pixel_size # arcsec
    group_meta_dict['image_height'] = step_ny * len(grid_dec) * pixel_size # arcsec
    group_meta_dict['pixel_size'] = pixel_size # arcsec
    group_meta_dict['n_images'] = len(image_files)
    group_meta_dict['image_indices_files'] = [(t, image_files[t]) for t in range(len(image_files))]
    
    group_table_dict = OrderedDict()
    group_table_dict['id'] = []
    group_table_dict['col'] = []
    group_table_dict['row'] = []
    group_table_dict['naxis1'] = []
    group_table_dict['naxis2'] = []
    group_table_dict['crval1'] = []
    group_table_dict['crval2'] = []
    group_table_dict['crpix1'] = []
    group_table_dict['crpix2'] = []
    group_table_dict['cdelt1'] = []
    group_table_dict['cdelt2'] = []
    group_table_dict['crota2'] = []
    group_table_dict['ra1'] = []
    group_table_dict['dec1'] = []
    group_table_dict['ra2'] = []
    group_table_dict['dec2'] = []
    group_table_dict['pad_ra1'] = []
    group_table_dict['pad_ra2'] = []
    group_table_dict['pad_dec1'] = []
    group_table_dict['pad_dec2'] = []
    group_table_dict['pad_x1'] = []
    group_table_dict['pad_x2'] = []
    group_table_dict['pad_y1'] = []
    group_table_dict['pad_y2'] = []
    group_table_dict['in_x1'] = []
    group_table_dict['in_x2'] = []
    group_table_dict['in_y1'] = []
    group_table_dict['in_y2'] = []
    group_table_dict['out_x1'] = []
    group_table_dict['out_x2'] = []
    group_table_dict['out_y1'] = []
    group_table_dict['out_y2'] = []
    group_table_dict['group_dir'] = []
    group_table_dict['n_images'] = []
    group_table_dict['image_indices'] = []
    group_table_dict['image_files'] = []
    #kdigit = int(np.log10(len(grid_ra)*len(grid_dec)))+1
    jdigit = int(np.log10(len(grid_dec)))+1+1
    idigit = int(np.log10(len(grid_ra)))+1+1
    #kfmt = '{:0%dd}'%(kdigit)
    jfmt = '{:0%dd}'%(jdigit)
    ifmt = '{:0%dd}'%(idigit)
    group_id = 0
    for j in range(len(grid_dec)):
        for i in range(len(grid_ra)):
            group_id += 1
            origin_ra = grid_ra[i]
            origin_dec = grid_dec[j]
            cell_vertices = [
                [origin_ra,         origin_dec],
                [origin_ra+step_ra, origin_dec],
                [origin_ra+step_ra, origin_dec+step_dec],
                [origin_ra,         origin_dec+step_dec],
                [origin_ra,         origin_dec]
            ]
            cell_polygon = Polygon(
                cell_vertices
            )
            col = jfmt.format(i+1) # str
            row = ifmt.format(j+1) # str
            group_dir = f'{output_name}_{group_id}_col{col}_row{row}'
            where_intersects = np.array([t.intersects(cell_polygon) for t in footprints])
            image_indices = np.argwhere(where_intersects).ravel().tolist()
            image_files_in_group = (np.array(image_files)[where_intersects]).tolist()
            footprints_in_group = np.take(footprints, image_indices)
            max_ra_in_group = np.max([t.max_ra() for t in footprints_in_group])
            min_ra_in_group = np.min([t.min_ra() for t in footprints_in_group])
            min_dec_in_group = np.min([t.min_dec() for t in footprints_in_group])
            max_dec_in_group = np.max([t.max_dec() for t in footprints_in_group])
            # compute the padding of image at left, bottom and right, top, 
            # to make sure we enclosed all images in this group so that 
            # the jwst pipeline resample step can drizzle all of them.
            # I added 5.0 arcsec (`5.0 / pixel_size`) to make sure all images are enclosed, 
            # because the jwst pipeline resample step seems will discard images 
            # which are partially enclosed!
            pad_ra1 = max(0.0, (max_ra_in_group - origin_ra))
            pad_ra2 = max(0.0, ((origin_ra + step_ra) - min_ra_in_group))
            pad_dec1 = max(0.0, (origin_dec - min_dec_in_group))
            pad_dec2 = max(0.0, (max_dec_in_group - (origin_dec + step_dec)))
            pad_x1 = int(np.ceil(pad_ra1 * np.cos(np.deg2rad(anchor_dec)) / pixel_size * 3600.0) + 5.0 / pixel_size)
            pad_x2 = int(np.ceil(pad_ra2 * np.cos(np.deg2rad(anchor_dec)) / pixel_size * 3600.0) + 5.0 / pixel_size)
            pad_y1 = int(np.ceil(pad_dec1 / pixel_size * 3600.0) + 5.0 / pixel_size)
            pad_y2 = int(np.ceil(pad_dec2 / pixel_size * 3600.0) + 5.0 / pixel_size) 
            group_table_dict['id'].append(group_id) # kfmt.format(group_id)
            group_table_dict['col'].append(col) # str
            group_table_dict['row'].append(row) # str
            group_table_dict['naxis1'].append(step_nx + pad_x1 + pad_x2)
            group_table_dict['naxis2'].append(step_ny + pad_y1 + pad_y2)
            group_table_dict['crval1'].append(anchor_ra)
            group_table_dict['crval2'].append(anchor_dec)
            group_table_dict['crpix1'].append(mosaic_crpix1 - i*step_nx + pad_x1)
            group_table_dict['crpix2'].append(mosaic_crpix2 - j*step_ny + pad_y1)
            group_table_dict['cdelt1'].append(-pixel_size / 3600.0) # negative
            group_table_dict['cdelt2'].append(pixel_size / 3600.0) # 
            group_table_dict['crota2'].append(0.0)
            group_table_dict['ra1'].append(origin_ra)
            group_table_dict['dec1'].append(origin_dec)
            group_table_dict['ra2'].append(origin_ra+step_ra)
            group_table_dict['dec2'].append(origin_dec+step_dec)
            group_table_dict['pad_ra1'].append(pad_ra1) # degrees
            group_table_dict['pad_ra2'].append(pad_ra2) # degrees
            group_table_dict['pad_dec1'].append(pad_dec1) # degrees
            group_table_dict['pad_dec2'].append(pad_dec2) # degrees
            group_table_dict['pad_x1'].append(pad_x1)
            group_table_dict['pad_x2'].append(pad_x2)
            group_table_dict['pad_y1'].append(pad_y1)
            group_table_dict['pad_y2'].append(pad_y2)
            group_table_dict['in_x1'].append(pad_x1)            # array indexing
            group_table_dict['in_x2'].append(pad_x1 + step_nx)  # array indexing
            group_table_dict['in_y1'].append(pad_y1)            # array indexing
            group_table_dict['in_y2'].append(pad_y1 + step_ny)  # array indexing, in[in_y1:in_y2,in_x1:in_x2]
            group_table_dict['out_x1'].append(i * step_nx)      # array indexing
            group_table_dict['out_x2'].append((i+1) * step_nx)  # array indexing
            group_table_dict['out_y1'].append(j * step_ny)      # array indexing
            group_table_dict['out_y2'].append((j+1) * step_ny)  # array indexing, out[out_y1:out_y2,out_x1:out_x2] = in_array[in_y1:in_y2,in_x1:in_x2]
            group_table_dict['group_dir'].append(group_dir)
            group_table_dict['n_images'].append(len(image_indices))
            group_table_dict['image_indices'].append(','.join(map(str,image_indices)))
            group_table_dict['image_files'].append(','.join(image_files_in_group))
    
    if save_group_meta_file is None:
        save_group_meta_file = f'{output_name}_meta.json'
    if save_group_meta_file != '':
        if os.path.isfile(save_group_meta_file):
            shutil.move(save_group_meta_file, save_group_meta_file+'.backup')
        else:
            save_group_table_dir = os.path.dirname(save_group_meta_file)
            if save_group_table_dir != '' and not os.path.isdir(save_group_table_dir):
                os.makedirs(save_group_table_dir)
        logger.info('Preparing mosaic image subgrouping meta and saving it to {!r}.'.format(save_group_meta_file))
        with open(save_group_meta_file, 'w') as fp:
            json.dump(group_meta_dict, fp, indent=4)
        logger.info('Saved group meta file as {!r}'.format(save_group_meta_file))
    
    group_table = Table(group_table_dict)
    
    # save gropu table file
    if save_group_table_file is None:
        save_group_table_file = f'{output_name}_table.txt'
    if save_group_table_file != '':
        if os.path.isfile(save_group_table_file):
            shutil.move(save_group_table_file, save_group_table_file+'.backup')
        else:
            save_group_table_dir = os.path.dirname(save_group_table_file)
            if save_group_table_dir != '' and not os.path.isdir(save_group_table_dir):
                os.makedirs(save_group_table_dir)
        logger.info('Building mosaic image subgrouping table and saving it to {!r}.'.format(save_group_table_file))
        if re.match(r'^(.*)\.(txt|dat)', save_group_table_file):
            group_table.write(save_group_table_file, format='ascii.commented_header')
            group_table.write(re.sub(r'^(.*)\.(txt|dat)', r'\1.csv', save_group_table_file), format='csv', overwrite=True) # also save csv file
        else:
            group_table.write(save_group_table_file)
        logger.info('Saved group table file as {!r}'.format(save_group_table_file))
    
    
    # save group region file
    if save_group_region_file is None:
        save_group_region_file = f'{output_name}_regions.reg'
    if save_group_region_file != '':
        write_ds9_region_file(
            group_table, 
            save_group_region_file,
        )
    
    
    # read catfile and prepare to copy/link to subgroup folders
    cat_dict = OrderedDict()
    if cat_file is not None and cat_file != '':
        with open(cat_file, 'r') as fp:
            for line in fp:
                if line.strip().startswith('#') or line.strip() == '':
                    continue
                key, val = line.split()
                cat_dict[key] = val
    
    
    # prepare subfolders and asn file
    group_asn_files = []
    group_cat_files = []
    for i in range(len(group_table)):
        group_id = group_table['id'][i]
        col = group_table['col'][i]
        row = group_table['row'][i]
        group_dir = group_table['group_dir'][i] # f'{output_name}_{group_id}_col{col}_row{row}'
        logger.info('Copying/linking image files into the subgroup directory {!r}'.format(group_dir))
        
        if not os.path.isdir(group_dir):
            os.makedirs(group_dir)
        
        image_files_in_group = group_table['image_files'][i].split(',')
        image_names_in_group = []
        cat_dict_in_group = OrderedDict()
        
        for image_file in image_files_in_group:
            image_name = os.path.basename(image_file)
            image_path = os.path.join(group_dir, image_name)
            if os.path.islink(image_path):
                os.remove(image_path)
            os.symlink(os.path.relpath(image_file, group_dir), image_path)
            image_names_in_group.append(image_name)
        
            if image_name in cat_dict:
                cat_name = os.path.basename(cat_dict[image_name])
                cat_path = os.path.join(group_dir, cat_name)
                if os.path.islink(cat_path):
                    os.remove(cat_path)
                os.symlink(os.path.relpath(cat_dict[image_name], group_dir), cat_path)
                cat_dict_in_group[image_name] = cat_name
        
        # write catfile.txt
        if len(cat_dict_in_group) > 0:
            cat_file_in_group = os.path.join(group_dir, 'catfile.txt')
            if os.path.isfile(cat_file_in_group):
                shutil.move(cat_file_in_group, cat_file_in_group+'.backup')
            with open(cat_file_in_group, 'w') as fp:
                for key in cat_dict_in_group:
                    fp.write('{} {}\n'.format(key, cat_dict_in_group[key]))
        else:
            cat_file_in_group = ''
        
        # write asn.json
        asn_file_in_group = os.path.join(group_dir, 'asn.json')
        #asn_obj = asn_from_list(
        #    list(zip(image_names_in_group, ['science']*len(image_files_in_group))), 
        #    product_name=group_dir,
        #    with_exptype=True,
        #)
        asn_obj = asn_from_list(
            image_names_in_group,
            product_name=group_dir,
            with_exptype=False,
        )
        _file_name, serialized = asn_obj.dump()
        if os.path.isfile(asn_file_in_group):
            shutil.move(asn_file_in_group, asn_file_in_group+'.backup')
        with open(asn_file_in_group, 'w') as fp:
            fp.write(serialized)
        group_asn_files.append(asn_file_in_group)
        group_cat_files.append(cat_file_in_group)
    
    
    # return
    return group_meta_dict, group_table, group_asn_files, group_cat_files




@click.command()
@click.argument('image_files', nargs=-1, type=click.Path(exists=True))
@click.argument('output_name', nargs=1, type=click.Path(exists=False))
@click.option('--asn-file', type=click.Path(exists=True), default=None, help='We can input an asn file instead of image files.')
@click.option('--cat-file', '--cat-file', 'cat_file', type=click.Path(exists=True), default=None, help='We can input a tweakreg catfile.')
@click.option('--grid-step', type=float, default=10.0, help='Gridding step in arcminutes.')
@click.option('--anchor-ra', type=float, default=None, help='Anchor point RA (degrees) in tangential projection. None for full mosaic image center.')
@click.option('--anchor-dec', type=float, default=None, help='Anchor point Dec (degrees) in tangential projection. None for full mosaic image center.')
@click.option('--pixel-size', type=float, default=1.0, help='Imaging pixel size in arcsec. Default is 1.0.')
@click.option('--cpu-cores', type=str, default='all', help='Multiprocessing cpu cores.')
@click.option('--save-group-meta-file', type=click.Path(exists=False), default=None, help='Output group meta. If None the output to \'{output_name}_meta.json\'')
@click.option('--save-group-table-file', type=click.Path(exists=False), default=None, help='Output group table. If None the output to \'{output_name}_table.txt\'')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False, help='Overwrite existing files.')
@click.option('--verbose/--no-verbose', is_flag=True, default=True, help='Verbose screen output.')
def main(
        image_files, 
        output_name, 
        asn_file, 
        cat_file, 
        grid_step, 
        anchor_ra, 
        anchor_dec, 
        pixel_size, 
        cpu_cores, 
        save_group_meta_file, 
        save_group_table_file, 
        overwrite, 
        verbose, 
    ):


    run_mosaic_image_subgrouping(
        image_files, 
        output_name, 
        asn_file = asn_file, 
        cat_file = cat_file, 
        grid_step = grid_step, 
        anchor_ra = anchor_ra, 
        anchor_dec = anchor_dec, 
        pixel_size = pixel_size, 
        cpu_cores = cpu_cores, 
        save_group_meta_file = save_group_meta_file, 
        save_group_table_file = save_group_table_file, 
        overwrite = overwrite, 
        verbose = verbose, 
    )




if __name__ == '__main__':
    
    main()



