#!/usr/bin/env python
# 
"""
Run outlier detection in parallel for groups of overlaping images. 

Usage: 
    (must be used in Python as a module)
    from util_run_outlier_detection_in_parallel import run_outlier_detection_in_parallel
    run_outlier_detection_in_parallel(pipeline_object, image_models)

By Daizhong Liu @MPE. 

Last updates: 
    2023-01-05

"""
import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict

from stpipe import Pipeline
from jwst import datamodels
from jwst.datamodels import ModelContainer 
from jwst.outlier_detection.outlier_detection_step import OutlierDetectionStep
from jwst.pipeline.calwebb_image3 import Image3Pipeline

import multiprocessing as mp

from shapely.geometry import Point, Polygon

# code name and version
CODE_NAME = 'util_run_outlier_detection_in_parallel.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230105'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)



def write_ds9_region_file(
        image_models:ModelContainer, 
        core_indices:list, 
        main_indices:list, 
        full_indices:list, 
        region_file:str = 'ds9.reg',
    ):
    with open(region_file, 'w') as fp:
        fp.write('# Region file format: DS9 version 4.1\n')
        fp.write('global color=green\n')
        fp.write('fk5\n')
        for i in core_indices:
            vertices = image_models[i].meta.wcs.footprint().ravel().tolist()
            vertices.extend([vertices[0], vertices[1]])
            fp.write('polygon({}) # text={{[{}]({}){}}} color={}\n'.format(
                    ', '.join(map(str, vertices)),
                    i, 'core', image_models[i].meta.filename, 
                    'cyan'
                )
            )
        for i in list(set(main_indices)-set(core_indices)):
            vertices = image_models[i].meta.wcs.footprint().ravel().tolist()
            vertices.extend([vertices[0], vertices[1]])
            fp.write('polygon({}) # text={{[{}]({}){}}} color={}\n'.format(
                    ', '.join(map(str, vertices)),
                    i, 'main', image_models[i].meta.filename, 'green'
                )
            )
        for i in list(set(full_indices)-set(main_indices)):
            vertices = image_models[i].meta.wcs.footprint().ravel().tolist()
            vertices.extend([vertices[0], vertices[1]])
            fp.write('polygon({}) # text={{[{}]({}){}}} color={} dash=1\n'.format(
                    ', '.join(map(str, vertices)),
                    i, 'edge', image_models[i].meta.filename, 'orange'
                )
            )
    


def run_one_outlier_detection_group(
        pipeline_object:Pipeline, 
        image_models:ModelContainer, 
        core_indices:list, 
        main_indices:list, 
        full_indices:list, 
    ):
    """
    Function to run one outlier_detection group. 
    """
    main_indices_str = '_'.join(map(str, main_indices[0:10]))
    working_dir = 'run_outlier_detection_' + main_indices_str
    pipeline_object.log.info('Running run_one_outlier_detection_group in ' + working_dir)
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)
    current_dir = os.getcwd()
    os.chdir(working_dir)
    outlier_detection = OutlierDetectionStep('outlier_detection_' + main_indices_str) # parent = pipeline_object
    outlier_detection.pixfrac = pipeline_object.outlier_detection.pixfrac
    outlier_detection.in_memory = False
    outlier_detection.output_file = os.path.basename(pipeline_object.output_file)
    outlier_detection.save_results = True
    outlier_detection.suffix = 'crf'
    outlier_detection.search_output_file = False
    write_ds9_region_file(image_models, core_indices, main_indices, full_indices, 'footprint.ds9.reg')
    running_indices = full_indices
    running_image_models = ModelContainer()
    updating_which = []
    for i in running_indices:
        running_image_models.append(copy.copy(image_models[i])) # TODO: copy.copy or not?
        updating_which.append(i in main_indices)
    updating_which = np.array(updating_which)
    message_str = ('Running outlier_detection_{} with {} image models: {} '
                   '(indices: {})'
                   .format(
                     main_indices_str, len(running_image_models), running_image_models, 
                     running_indices
                   )
                  )
    pipeline_object.log.info(message_str)
    running_image_models = outlier_detection.run(running_image_models)
    np.argwhere(full_indices)
    updating_image_models = ModelContainer()
    for kk in range(len(running_image_models)):
        if updating_which[kk]:
            updating_image_models.append(running_image_models[kk])
    os.chdir(current_dir)
    return updating_image_models




def run_outlier_detection_in_parallel(
        pipeline_object:Pipeline, 
        image_models:ModelContainer, 
        core_image_per_group = 1, # minimum core image(s) of each group to start with
        main_image_min_overlap = 0.5, # expand the main image list of each group with images overlapping by more than this fraction with 
        use_pipeline_logger = True, 
        max_cores = 'all', 
        save_group_table_file = '', 
        return_group_table = False, 
        overwrite = False, 
        verbose = True, 
    ):
    """
    Run outlier detection in groups of images by their overlapping fractions.
    
    We will group the images in this way. 
    - First, we start from the first index of the input `image_models` list, then 
      we continue to include more 'core' images depending on the `core_per_group` number. 
      In this way we build a `core_footprint` and a `core_indicies` list.
    - Next, we go over the remain image list for each group and find out images 
      that are overlapped with the `core_footprint`. Those overlapped by at least this 
      `main_image_min_overlap` fraction are considered as 'main' images, and those with 
      less overlapping are sorted as 'edge' images. We build a `main_indicies` list 
      by extending the `core_indicies` with the main overlapped images, and a `edge_indicies` 
      list with all other overlapped images. 
    - Then, we send the main+core+edge list to the outlier_detection
    
    Args:
        
        core_image_per_group :
            minimum number of core image(s) in each group to start with
        
        main_image_min_overlap :
            expand the main image list of each group with images overlapping 
            by more than this fraction with the core images' footprint
    
    """
    
    # set logger
    if use_pipeline_logger:
        this_logger = pipeline_object.log
    else:
        this_logger = logger
    
    
    # load existing 'save_group_table_file'
    group_table = None
    if save_group_table_file != '':
        if os.path.isfile(save_group_table_file) and not overwrite: 
            this_logger.info('Loading existing group table file {!r}'.format(save_group_table_file))
            if re.match(r'^.*\.(txt|dat)', save_group_table_file):
                group_table = Table.read(save_group_table_file, format='ascii.commented_header')
            else:
                group_table = Table.read(save_group_table_file)
            if 'core_indices' not in group_table.colnames or \
               'main_indices' not in group_table.colnames or \
               'edge_indices' not in group_table.colnames or \
               'full_indices' not in group_table.colnames:
                group_table = None
            group_core_indices = [list(map(int, t.split(','))) for t in group_table['core_indices']]
            group_main_indices = [list(map(int, t.split(','))) for t in group_table['main_indices']]
            group_edge_indices = [list(map(int, t.split(','))) for t in group_table['edge_indices']]
            group_full_indices = [list(map(int, t.split(','))) for t in group_table['full_indices']]
    
    if group_table is None:
        
        # make footprints
        image_footprints = []
        for i in range(len(image_models)):
            message_str = ('Extracting image_models[{}] @{} {!r} footprint'
                           .format(
                             i, hex(id(image_models[i])), image_models[i].meta.filename, 
                           )
                          )
            this_logger.info(message_str)
            image_footprint = Polygon(image_models[i].meta.wcs.footprint())
            image_footprints.append(image_footprint)
        
        # determine groups of images
        # We start from the first index of the image list, go through the list to take 
        # 'core_image_per_group'-number of overlapping images, and build a 'core_footprint'.
        # Then we go over the remain image list to find all overlapping images, and
        # put images that are >50% overlaped also into a main-image list. 
        # Then we base on the main image footprint and go over again the remain image list
        # to find all overlapped images, thus building a full-image list to process in this group. 
        # Then we will call the outlier_detection, but we will only update the results
        # for images in the main-image list. 
        group_core_indices = []
        group_main_indices = []
        group_edge_indices = []
        group_full_indices = []
        all_indices = np.arange(len(image_models)).tolist()
        g = 0
        while len(all_indices) > 0:
            first_index = all_indices.pop(0)
            core_indices = []
            core_indices.append(first_index)
            core_footprint = copy.copy(image_footprints[first_index])
            k = 0
            while k < len(all_indices) and len(all_indices) > 0 and len(core_indices) < core_image_per_group:
                next_index = all_indices[k]
                next_footprint = image_footprints[next_index]
                if core_footprint.intersects(next_footprint):
                    all_indices.pop(k)
                    core_indices.append(next_index)
                    core_footprint = core_footprint.union(next_footprint)
                else:
                    k+=1
            # 
            # now we have core_indices and core_footprint for this group
            # we need to find all other overlapped images
            if len(core_indices) > 0:
                # 
                # go through the remain `all_indices` list to sort out main overlapped images
                # whose overlapping fraction with the `core_footprint` is larger than the `main_image_min_overlap` fraction.
                main_indices = copy.copy(core_indices)
                main_footprint = copy.copy(core_footprint)
                k = 0
                while k < len(all_indices) and len(all_indices) > 0:
                    next_index = all_indices[k]
                    next_footprint = image_footprints[next_index]
                    overlap_fraction = core_footprint.intersection(next_footprint).area / next_footprint.area
                    this_logger.info('checking overlap fraction between core_footprint ({}) and next_footprint ({}): {}'.format(
                        repr(core_indices), next_index, overlap_fraction))
                    if overlap_fraction >= main_image_min_overlap:
                        all_indices.pop(k)
                        main_indices.append(next_index)
                        main_footprint = main_footprint.union(next_footprint)
                        this_logger.info('growing main indices ({}): {}'.format(len(main_indices), repr(main_indices)))
                    else:
                        k+=1
                # 
                # go through the remain `all_indices` list to sort out edge overlapped images
                # which overlap with the `main_footprint`.
                edge_indices = []
                full_indices = copy.copy(main_indices)
                k = 0
                while k < len(all_indices) and len(all_indices) > 0:
                    next_index = all_indices[k]
                    next_footprint = image_footprints[next_index]
                    if main_footprint.intersects(next_footprint):
                        all_indices.pop(k)
                        edge_indices.append(next_index)
                        full_indices.append(next_index)
                    else:
                        k+=1
                # 
                group_core_indices.append(core_indices)
                group_main_indices.append(main_indices)
                group_edge_indices.append(edge_indices)
                group_full_indices.append(full_indices)
                # 
                if verbose:
                    this_logger.info('group_core_indices[{}] = {}'.format(g, repr(core_indices)))
                    this_logger.info('group_main_indices[{}] = {}'.format(g, repr(main_indices)))
                    this_logger.info('group_full_indices[{}] = {}'.format(g, repr(full_indices)))
                # 
                g+=1
        
        
        # if return_group_table then do not run any processing but just return the group table
        group_table_dict = OrderedDict()
        group_table_dict['group_id'] = []
        group_table_dict['core_indices'] = []
        group_table_dict['main_indices'] = []
        group_table_dict['main_filenames'] = []
        group_table_dict['edge_indices'] = []
        group_table_dict['edge_filenames'] = []
        group_table_dict['full_indices'] = []
        for k in range(len(group_full_indices)):
            core_indices = group_core_indices[k]
            main_indices = group_main_indices[k]
            edge_indices = group_edge_indices[k]
            full_indices = group_full_indices[k]
            main_filenames = [image_models[t].meta.filename for t in main_indices]
            edge_filenames = [image_models[t].meta.filename for t in edge_indices]
            group_table_dict['group_id'].append(k+1)
            group_table_dict['core_indices'].append(','.join(map(str,core_indices)))
            group_table_dict['main_indices'].append(','.join(map(str,main_indices)))
            group_table_dict['main_filenames'].append(','.join(main_filenames))
            group_table_dict['edge_indices'].append(','.join(map(str,edge_indices)))
            group_table_dict['edge_filenames'].append(','.join(edge_filenames))
            group_table_dict['full_indices'].append(','.join(map(str,full_indices)))
        group_table = Table(group_table_dict)
        if save_group_table_file != '':
            if os.path.isfile(save_group_table_file):
                if overwrite:
                    shutil.move(save_group_table_file, save_group_table_file+'.backup')
            else:
                save_group_table_dir = os.path.dirname(save_group_table_file)
                if save_group_table_dir != '' and not os.path.isdir(save_group_table_dir):
                    os.makedirs(save_group_table_dir)
            this_logger.info('Building outlier_detection group table and saving it to {!r}.'.format(save_group_table_file))
            if re.match(r'^.*\.(txt|dat)', save_group_table_file):
                group_table.write(save_group_table_file, format='ascii.commented_header')
            else:
                group_table.write(save_group_table_file)
    
    if return_group_table:
        this_logger.info('Building outlier_detection group table and returning it without actually running outlier_detection.')
        return group_table
    
    
    # determine number of parallel processes
    n_parallel = mp.cpu_count()
    if isinstance(max_cores, int):
        n_parallel = max_cores
    if re.match(r'^[0-9]+$', max_cores):
        n_parallel = int(max_cores)
    elif max_cores == 'all':
        n_parallel = mp.cpu_count()
    elif max_cores == 'half':
        n_parallel = mp.cpu_count()//2
    elif max_cores == 'quarter':
        n_parallel = mp.cpu_count()//4
    elif max_cores == 'none':
        n_parallel = 1
    if n_parallel == 0:
        n_parallel = 1
    this_logger.info('Running outlier_detection in parallel with n_parallel = {}'.format(n_parallel))
    
    
    # iterate run_one_outlier_detection_group
    using_multiprocessing = False # TODO can not be True because something can not be pickled -- CANNOT DO PARALLEL!
    if using_multiprocessing:
        pool = mp.Pool(n_parallel)
        k = 0
        rets = []
        for k in range(len(group_table)):
            core_indices = group_core_indices[k]
            main_indices = group_main_indices[k]
            edge_indices = group_edge_indices[k]
            full_indices = group_full_indices[k]
            message_str = ('Running outlier_detection group {}/{} '
                           '(main indices: {}, edge indices: {})'
                           .format(
                             k+1, len(group_table), 
                             main_indices, edge_indices,
                           )
                          )
            this_logger.info(message_str)
            ret = pool.apply_async(run_one_outlier_detection_group, 
                                   args = (pipeline_object,
                                           image_models,
                                           core_indices,
                                           main_indices,
                                           full_indices),
            )
            rets.append(ret)
        pool.close()
        pool.join()
        for k in range(len(group_table)):
            updated_image_models = rets[k].get()
            for kk in range(len(main_indices)):
                message_str = ('Updating image_models[{}] @{} {!r} from outlier_detection group {}/{}'
                               .format(
                                 main_indices[kk], hex(id(image_models[main_indices[kk]])), 
                                 image_models[main_indices[kk]].meta.filename, 
                                 k+1, len(group_table), 
                               )
                              )
                this_logger.info(message_str)
                image_models[main_indices[kk]] = updated_image_models[kk]
    else:
        for k in range(len(group_table)):
            core_indices = group_core_indices[k]
            main_indices = group_main_indices[k]
            edge_indices = group_edge_indices[k]
            full_indices = group_full_indices[k]
            message_str = ('Running outlier_detection group {}/{} '
                           '(main indices: {}, edge indices: {})'
                           .format(
                             k+1, len(group_table), 
                             main_indices, edge_indices,
                           )
                          )
            this_logger.info(message_str)
            updated_image_models = \
                run_one_outlier_detection_group(pipeline_object,
                                                image_models,
                                                core_indices,
                                                main_indices,
                                                full_indices)
            
            for kk in range(len(main_indices)):
                message_str = ('Updating image_models[{}] @{} {!r} from outlier_detection group {}/{}'
                               .format(
                                 main_indices[kk], hex(id(image_models[main_indices[kk]])), 
                                 image_models[main_indices[kk]].meta.filename, 
                                 k+1, len(group_table), 
                               )
                              )
                this_logger.info(message_str)
                image_models[main_indices[kk]] = updated_image_models[kk]
        
    
    # return
    return image_models




@click.command()
@click.argument('input_asn_file', type=click.Path(exists=True))
@click.argument('output_name', type=click.Path(exists=False))
@click.option('--core-image-per-group', type=int, default=1, help='The number of core images to start with. We will grow the coverage by overlapped fraction.')
@click.option('--main-image-min-overlap', type=float, default=0.5, help='The minimum overlapping fraction of each main image to the union of core images.')
@click.option('--max-cores', type=str, default='none', help='This parameter is not used because multiprocessing does not work. Something cannot be pickled.')
@click.option('--save-group-table-file', type=click.Path(exists=False), default='run_outlier_detection_group_table.txt', help='An output table listing the determined groups of images.')
@click.option('--return-group-table', is_flag=True, default=False, help='Only return the group table then exit. Will not actually run the outlier_detection step.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False, help='Overwrite existing files.')
@click.option('--verbose/--no-verbose', is_flag=True, default=True, help='Verbose screen output.')
def main(
        input_asn_file, 
        output_name, 
        core_image_per_group, 
        main_image_min_overlap, 
        max_cores, 
        save_group_table_file, 
        return_group_table, 
        overwrite, 
        verbose, 
    ):

    pipeline_object = Image3Pipeline()
    
    pipeline_object.output_file = output_name # os.path.splitext(output_file)[0]
    pipeline_object.output_ext = ".fits" # default
    pipeline_object.save_results = True
    
    pipeline_object.outlier_detection.pixfrac = 0.48
    pipeline_object.outlier_detection.in_memory = False
    
    asn_filename = input_asn_file
    asn_exptypes = ['science']
    with datamodels.open(asn_filename, asn_exptypes=asn_exptypes) as input_models:
        run_outlier_detection_in_parallel(
            pipeline_object, 
            input_models, 
            core_image_per_group = core_image_per_group, 
            main_image_min_overlap = main_image_min_overlap, 
            max_cores = max_cores, 
            save_group_table_file = save_group_table_file, 
            return_group_table = return_group_table, 
            overwrite = overwrite, 
            verbose = verbose, 
        )




if __name__ == '__main__':
    
    main()



