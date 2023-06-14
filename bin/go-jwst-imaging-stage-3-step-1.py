#!/usr/bin/env python
#
"""
Process a number of JWST cal images, output a mosaic image. 

This script runs the JWST calwebb_image3 pipeline.

Inputs
    
    jw*_cal.fits

Outputs

    jw*_cal.fits
    jw*_i2d.fits
    
Last update:
    
    2022-09-10 DZLIU.
    2022-12-03 DZLIU. tweakreg with "*_cat_for_tweakreg.csv".


More notes from CEERS in "ceers_nircam_reduction.ipynb":
calwebb_image3 - Ensemble Calibrations

    The Stage 3 calwebb_image3 pipeline takes one or more calibrated slope images 
    (*_cal.fits files) and combines them into a final mosaic image. It then creates 
    a source catalog from this mosaic. Several steps are performed in order to prepare 
    the data for the mosaic creation. These steps largely mirror what is done by 
    DrizzlePac software when working with HST data.

    First, using common sources found across the input images, the WCS of each image is 
    refined. Background levels are then matched across the inputs. Spurious sources 
    (e.g. cosmic rays that were not flagged in the Jump step during Stage 1 processing) 
    are removed by comparing each individual input image to a median image. The indivudal 
    images are combined into a single mosaic image. A source catalog is created based on 
    the mosaic image. And finally, the individual exposures are updated using the 
    information from the preceding steps. New versions of the individual calibrated slope 
    images are produced that contain matched backgrounds, flagged spurious sources, and 
    improved WCS objects.

    All JWST imaging mode data, regardless of instrument, are processed through the 
    calwebb_image3 pipeline. The steps and the order in which they are performed is 
    the same for all data. The pipeline is a wrapper which will string together all 
    of the appropriate steps in the proper order. See Figure 1 on the calwebb_image3 
    algorithm page for a map of the steps that are performed on the input data.

Inputs

    2D calibrated images (*_cal.fits), organized in an ASN file

Outputs

    2D cosmic-ray flagged images (*_crf.fits), created during the OutlierDetection step
    2D resampled, combined mosaic image (*_i2d.fits) including all exposures in the association, created during the Resample step
    2D segmentation map (*_segm.fits) based on the *_i2d.fits image, created by the SourceCatalog step
    Catalog of photometry (*_cat.escv) saved as an ASCII file in ecsv format, created by the SourceCatalog step
    
From ceers_nircam_reduction.ipynb

"""

# Packages that allow us to get information about objects:
import os, sys, re, gc, json, copy, datetime, time, glob, shutil, gc
import click
from collections import OrderedDict
try:
    from packaging.version import parse as LooseVersion
except:
    from distutils.version import LooseVersion

# Numpy library:
import numpy as np

# Astropy tools:
from astropy.io import fits
from astropy.table import Table

# Import JWST pipeline-related modules

# List of possible data quality flags
from jwst.datamodels import dqflags

# The entire calwebb_image3 pipeline
from jwst.pipeline import calwebb_image3

# Individual steps that make up calwebb_image3
from jwst.tweakreg import TweakRegStep
from jwst.skymatch import SkyMatchStep
from jwst.outlier_detection import OutlierDetectionStep
from jwst.resample import ResampleStep
from jwst.source_catalog import SourceCatalogStep
from jwst import datamodels
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.associations import load_asn

# Import jwst package itself
import jwst

# Import stpipe.Pipeline
from stpipe import Pipeline

# Import AsdfFile
from asdf import AsdfFile

# Setup logging
import logging

# Define utility functions
def get_script_dir():
    """Get current script file's directory path."""
    return os.path.abspath(os.path.dirname(__file__))

def get_script_name():
    """Get current script file name without the suffix and replaced some characters to underscores."""
    return re.sub(r'[^a-zA-Z0-9_]', r'_', os.path.splitext(os.path.basename(__file__))[0])

def setup_logger():
    logger_streamhandler = logging.StreamHandler()
    logger_streamhandler_formatter = logging.Formatter("[%(asctime)-8s] %(message)s", "%H:%M:%S")
    logger_streamhandler.setFormatter(logger_streamhandler_formatter)
    logger_streamhandler.setLevel(logging.DEBUG)

    log_file = get_script_name()
    log_time = datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    logger_filehandler = logging.FileHandler(f"log_{log_file}_{log_time}.txt", mode='a')
    logger_filehandler_formatter = logging.Formatter("[%(asctime)-15s] %(message)s", "%Y-%m-%d %H:%M:%S")
    logger_filehandler.setFormatter(logger_filehandler_formatter)

    logger = logging.getLogger()
    while len(logger.handlers) > 0:
        del logger.handlers[0]
    logger.addHandler(logger_streamhandler)
    logger.addHandler(logger_filehandler)
    logger.setLevel(logging.DEBUG)
    
    return logger

def expand_list_with_comma_items(input_list):
    if isinstance(input_list, str):
        input_list = [input_list]
    output_list = []
    for item in input_list:
        if not isinstance(item, str):
            item = str(item)
        if item.find(',') >= 0:
            output_list.extend(item.split(','))
        else:
            output_list.append(item)
    return output_list


def link_files_to_current_directory(
        file_list, 
        current_directory = None, 
    ):
    """Input abspath, make symlinks in current directory, output the file names without paths.
    """
    if current_directory is None:
        current_directory = os.getcwd()
    name_list = []
    for file_path in file_list:
        if os.path.dirname(file_path) != '': # contains a path
            file_name = os.path.basename(file_path)
            if os.path.islink(file_name):
                os.remove(file_name)
            os.symlink(os.path.relpath(file_path, current_directory),
                       os.path.join(current_directory, file_name))
        else:
            file_name = file_path
        name_list.append(file_name)
    return name_list


def asn_from_list_to_file(
        image_list, 
        asn_file,
        product_name = None,
        asn_id = None, 
        clear_asn_id = False, 
        rule = DMS_Level3_Base, 
        #TODO: exptype
    ):
    name_list = link_files_to_current_directory(image_list)
    if product_name is None: # if product_name is not given, use current directory name.
        product_name = os.path.basename(os.getcwd())
    #asn_object = asn_from_list(
    #    list(zip(name_list, ['science']*len(name_list))), 
    #    product_name = product_name,
    #    with_exptype = True
    #)
    asn_object = asn_from_list(
        name_list, 
        product_name = product_name,
        rule = rule, 
    )
    if clear_asn_id and 'asn_id' in asn_object:
        del asn_object['asn_id']
    elif asn_id is not None:
        asn_object['asn_id'] = asn_id
    _file_name, serialized = asn_object.dump()
    if os.path.isfile(asn_file):
        shutil.move(asn_file, asn_file+'.backup')
    with open(asn_file, 'w') as fp:
        fp.write(serialized)
    return asn_file


def asdf_from_step_to_file(
        step_object, 
        asdf_file,
    ):
    asdf_object = AsdfFile(step_object.get_pars())
    if os.path.isfile(asdf_file):
        shutil.move(asdf_file, asdf_file+'.backup')
    asdf_object.write_to(asdf_file)
    
    # 20230318
    # see -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/stpipe/config_asdf.html#config-asdf-files
    asdf_file2 = os.path.splitext(asdf_file)[0] + '_export_config.asdf'
    if os.path.isfile(asdf_file2):
        shutil.move(asdf_file2, asdf_file2+'.backup')
    step_object.export_config(asdf_file2) # include_meta=True
    
    return asdf_file


def imcat_dict_to_catfile(
        imcat_dict, 
        image_list,
        cat_file,
    ):
    catfile_lines = []
    for image_file in image_list:
        if image_file in imcat_dict:
            imcat_file = imcat_dict[image_file]
            (imcat_name,) = link_files_to_current_directory((imcat_file,))
            (image_name,) = link_files_to_current_directory((image_file,))
            catfile_lines.append('{} {}'.format(image_name, imcat_name))
    if os.path.isfile(cat_file):
        shutil.move(cat_file, cat_file+'.backup')
    if len(catfile_lines) > 0:
        with open(cat_file, 'w') as fp:
            for line in catfile_lines:
                fp.write(line+'\n')
        return cat_file
    else:
        return ''


def clean_up_intermediate_outlier_i2d_files(
        image_list, 
        output_name, 
        asn_id = None, 
        reload_image_models = False, 
    ):
    # move some output files to subdirectories
    # for example the "outlier_i2d.fits", so that we do not find them when doing an `ls *_i2d.fits`
    reloading_image_model_files = []
    for i in range(len(image_list)):
        image_file = image_list[i]
        # add all possible files to move
        files_to_move = []
        files_to_move.append(f'{output_name}_{i}_outlier_i2d.fits')
        files_to_move.append(f'asn_{output_name}_{i}_blot.fits')
        if asn_id is not None:
            files_to_move.append(f'{output_name}_{i}_{asn_id}_outlier_i2d.fits')
            files_to_move.append(f'{output_name}_{asn_id}_{i}_outlier_i2d.fits')
            files_to_move.append(f'asn_{output_name}_{i}_{asn_id}_blot.fits')
            files_to_move.append(f'asn_{asn_id}_{output_name}_{i}_blot.fits')
        dir_move_to = f'{output_name}_{i}_outlier'
        for file_to_move in files_to_move:
            if os.path.isfile(file_to_move):
                if not os.path.isdir(dir_move_to):
                    os.makedirs(dir_move_to)
                shutil.move(file_to_move, 
                    os.path.join(dir_move_to, os.path.basename(file_to_move)))
        # 
        # also rename files which have f'{output_name}_{asn_id}_{i}_outlierdetection.fits'
        # as f'{output_name}_{i}_{asn_id}_outlierdetection.fits'
        if asn_id is not None:
            if os.path.isfile(f'{output_name}_{asn_id}_{i}_outlierdetection.fits'):
                shutil.move(f'{output_name}_{asn_id}_{i}_outlierdetection.fits', 
                            f'{output_name}_{i}_{asn_id}_outlierdetection.fits')
        if reload_image_models:
            reloading_image_model_files.append(f'{output_name}_{i}_{asn_id}_outlierdetection.fits')
    if reload_image_models: # not tested 20230412
        asn_from_list_to_file(reloading_image_model_files, 'asn_outlier_detection_reloading.json')
        reloaded_image_models = datamodels.ModelContainer()
        asn_data = reloaded_image_models.read_asn('asn_outlier_detection_reloading.json')
        reloaded_image_models.from_asn(asn_data)
        return reloaded_image_models


def run_individual_steps_for_one_asn_file(
        pipeline_object, 
        asn_filename, 
        output_name, 
    ):
    
    pipeline_object.log.info(
        f'*-*-*\n*-*-* run_individual_steps_for_one_asn_file *-*-*\n*-*-*'
    )
    pipeline_object.log.info(
        f'Step {pipeline_object.name} running with asn file {asn_filename}.'
    )
    pipeline_object.log.info(
        f'Step {pipeline_object.name} parameters are: {pipeline_object.get_pars()}'
    )
    
    pipeline_object.output_name = output_name
    
    with datamodels.open(asn_filename, asn_exptypes=['science']) as input_models:
        # 
        # 1. tweakreg
        pipeline_object.tweakreg.output_file = output_name
        pipeline_object.tweakreg.save_results = True
        pipeline_object.tweakreg.search_output_file = False
        image_models = pipeline_object.tweakreg(input_models)
        # TODO: sometimes tweakreg update_fits_wcsinfo fails because `max_pix_error` is too small, 
        # we can run this manually with a larger `max_pix_error`, then save_model.
        from jwst.assign_wcs.util import update_fits_wcsinfo
        for image_model in image_models:
            update_fits_wcsinfo(image_model,
                max_pix_error=1.,
                npoints=16)
        pipeline_object.tweakreg.save_model(image_models, 
            idx=None, suffix='tweakreg',
            format=pipeline_object.tweakreg.name_format, force=True)
        # 
        # 2. skymatch bkgmatch
        pipeline_object.skymatch.output_file = output_name
        pipeline_object.skymatch.save_results = True
        pipeline_object.skymatch.search_output_file = False
        image_models = pipeline_object.skymatch(image_models)
        # 
        # 3. outlier detection
        pipeline_object.outlier_detection.output_file = output_name
        pipeline_object.outlier_detection.maskpt = 0.2 # default is 0.7
        pipeline_object.outlier_detection.save_intermediate_results = True # False
        pipeline_object.outlier_detection.save_results = True # will save to '{asn_name}_{i}_{asn_id}_crf.fits'
        pipeline_object.outlier_detection.suffix = 'outlierdetection' # default is crf
        pipeline_object.outlier_detection.search_output_file = False
        asn_id = 'a3001' # the default, see ...
        image_models = pipeline_object.outlier_detection(image_models)
        image_models = clean_up_intermediate_outlier_i2d_files(image_models, output_name, asn_id = asn_id, 
            reload_image_models = True)
        # 20230412 it seems pipeline does not return the outlier-detected-pixel-masked image models
        #          have to reload image_models
        # 
        # 4. resample
        pipeline_object.resample.output_file = output_name
        pipeline_object.resample.save_results = True
        pipeline_object.resample.suffix = 'i2d'
        result = pipeline_object.resample(image_models)
        # 
        # 5. source_catalog
        pipeline_object.source_catalog.output_file = output_name
        pipeline_object.source_catalog.search_output_file = False
        pipeline_object.source_catalog.save_results = True
        if isinstance(result, datamodels.ImageModel) and result.meta.cal_step.resample == 'COMPLETE':
            pipeline_object.source_catalog(result)
        # 
    # 
    # save pars as asdf file
    asdf_filepath = output_name + '_calwebb_image3.asdf'
    if os.path.isfile(asdf_filepath):
        shutil.move(asdf_filepath, asdf_filepath+'.backup')
    asdf_object = AsdfFile(pipeline_object.get_pars())
    asdf_object.write_to(asdf_filepath)
    pipeline_object.log.info('Parameters are saved into {}'.format(asdf_filepath))
    
    # 20230318
    #asdf_filepath = output_name + '_calwebb_image3_export_config.asdf'
    #if os.path.isfile(asdf_filepath):
    #    shutil.move(asdf_filepath, asdf_filepath+'.backup')
    #pipeline_object.export_config('cube_build.asdf')
    #pipeline_object.log.info('Parameters are saved into {}'.format(asdf_filepath))
    


def run_individual_steps_for_image_files(
        pipeline_object, 
        image_files, 
        output_name, 
        overwrite = False, 
    ):
    
    pipeline_object.log.info(
        f'*-*-*\n*-*-* run_individual_steps_for_image_files *-*-*\n*-*-*'
    )
    pipeline_object.log.info(
        f'Step {pipeline_object.name} parameters are: {pipeline_object.get_pars()}'
    )
    
    pipeline_object.output_name = output_name
    
    # 
    # 1. tweakreg
    pipeline_object.tweakreg.output_file = output_name
    pipeline_object.tweakreg.save_results = True # will save as '{dataset_name}_tweakreg.fits' (without '_cal')
    pipeline_object.tweakreg.search_output_file = False
    processing_image_files = image_files
    processed_image_files = [re.sub(r'_cal(\.fits|)$', r'', t)+'_tweakreg.fits' for t in processing_image_files] # (without '_cal')
    pipeline_object.log.info('Checking tweakreg output file existence: {}'.format(repr(processed_image_files)))
    if not np.all(list(map(os.path.exists, processed_image_files))):
        asn_from_list_to_file(processing_image_files, 'asn_tweakreg.json')
        asdf_from_step_to_file(pipeline_object.tweakreg, 'asdf_tweakreg.txt')
        image_models = pipeline_object.tweakreg('asn_tweakreg.json')
        # TODO: sometimes tweakreg update_fits_wcsinfo fails because `max_pix_error` is too small, 
        # we can run this manually with a larger `max_pix_error`, then save_model.
        from jwst.assign_wcs.util import update_fits_wcsinfo
        for image_model in image_models:
            update_fits_wcsinfo(image_model,
                max_pix_error=1.,
                npoints=16)
        pipeline_object.tweakreg.save_model(image_models, 
            idx=None, suffix='tweakreg', # will save as '{dataset_name}_tweakreg.fits' (without '_cal')
            format=pipeline_object.tweakreg.name_format, force=True)
        # 
        # Notes:
        #   dzliu fixing a potential tweak_step bug "images.from_asn(input)" --> "images.from_asn(asn_data)"
        # 
        # 
        # TODO: 20230318: 
        #   tweakreg was skipped for some images, error messages as below:
        #     ERROR - Number of output coordinates exceeded allocation
        #     ERROR - Multiple sources within specified tolerance matched to a single reference source. 
        #             Try to adjust 'tolerance' and/or 'separation' parameters.
        #   which will make tweakreg to skip the action at the step of aligning JWST images relatively
        #     (this is even before doing the absolute alignment to the absrefcat, i.e., 
        #       before `if align_to_abs_refcat` in "tweakreg_step.py"), 
        #     marking `model.meta.cal_step.tweakreg = "SKIPPED"`, 
        #     and leaving a log message `Skipping 'TweakRegStep'`.
        #   After extensive testing, it seems changing `tolerance` from 1.0 to 0.7 will do the trick...
        # 
    else:
        pipeline_object.log.info('Step tweakreg is skipped because all output files exist: {}'.format(repr(processed_image_files)))
    # 
    # 2. skymatch bkgmatch
    pipeline_object.skymatch.output_file = output_name
    pipeline_object.skymatch.save_results = True # will save as '{output_name}_{index}_cal_skymatch.fits'
    pipeline_object.skymatch.search_output_file = False
    processing_image_files = processed_image_files
    if len(processing_image_files) == 1:
        processed_image_files = [f'{output_name}_skymatch.fits']
    else:
        processed_image_files = [f'{output_name}_{i}_skymatch.fits' for i in range(len(processing_image_files))]
    pipeline_object.log.info('Checking skymatch output file existence: {}'.format(repr(processed_image_files)))
    if not np.all(list(map(os.path.exists, processed_image_files))):
        asn_from_list_to_file(processing_image_files, 'asn_skymatch.json')
        asdf_from_step_to_file(pipeline_object.skymatch, 'asdf_skymatch.txt')
        image_models = pipeline_object.skymatch('asn_skymatch.json')
        del image_models
        gc.collect()
    else:
        pipeline_object.log.info('Step skymatch is skipped because all output files exist: {}'.format(repr(processed_image_files)))
    # 
    # 3. outlier detection
    pipeline_object.outlier_detection.output_file = output_name
    pipeline_object.outlier_detection.save_intermediate_results = True # False
    pipeline_object.outlier_detection.save_results = True # will save as '{output_name}_{index}_{asn_id}_outlierdetection.fits'
    pipeline_object.outlier_detection.suffix = 'outlierdetection' # default is crf
    pipeline_object.outlier_detection.search_output_file = False
    asn_id = 'a3001' # the default, see ...
    processing_image_files = processed_image_files
    clean_up_intermediate_outlier_i2d_files(processing_image_files, output_name, asn_id=asn_id) #<TODO># temporary fix 20230110
    if len(processing_image_files) == 1:
        processed_image_files = [f'{output_name}_{asn_id}_outlierdetection.fits']
    else:
        processed_image_files = [f'{output_name}_{i}_{asn_id}_outlierdetection.fits' for i in range(len(processing_image_files))]
    pipeline_object.log.info('Checking outlier_detection output file existence: {}'.format(repr(processed_image_files)))
    if not np.all(list(map(os.path.exists, processed_image_files))):
        asn_from_list_to_file(processing_image_files, 'asn_outlier_detection.json', asn_id=asn_id)
        asdf_from_step_to_file(pipeline_object.outlier_detection, 'asdf_outlier_detection.txt')
        image_models = pipeline_object.outlier_detection('asn_outlier_detection.json')
        clean_up_intermediate_outlier_i2d_files(image_models, output_name, asn_id=asn_id)
        del image_models
        gc.collect()
        # if there is only a single file, make a copy
        if len(processing_image_files) == 1:
            pipeline_object.log.info('Step outlier_detection is skipped because there is only one file to process. '+\
                                     'Copying {!r} -> {!r}'.format(processing_image_files[0], processed_image_files[0]))
            shutil.copy2(processing_image_files[0], processed_image_files[0])
    else:
        pipeline_object.log.info('Step outlier_detection is skipped because all output files exist: {}'.format(repr(processed_image_files)))
    # 
    # 4. resample
    pipeline_object.resample.output_file = output_name
    pipeline_object.resample.save_results = True # will save as '{output_name}_i2d.fits'
    pipeline_object.resample.suffix = 'i2d'
    pipeline_object.resample.search_output_file = False # do not search the `output_file` variable for basepath
    processing_image_files = processed_image_files
    processed_image_files = [output_name+'_i2d.fits']
    pipeline_object.log.info('Checking resample output file existence: {}'.format(repr(processed_image_files)))
    if not np.all(list(map(os.path.exists, processed_image_files))):
        asn_from_list_to_file(processing_image_files, 'asn_resample.json', asn_id=asn_id)
        asdf_from_step_to_file(pipeline_object.resample, 'asdf_resample.txt')
        #result_model = pipeline_object.resample('asn_resample.json') #-- NOT WORKING? WHY?
        image_models = datamodels.ModelContainer()
        asn_data = image_models.read_asn('asn_resample.json')
        image_models.from_asn(asn_data) # see "tweakreg_step.py"
        result_model = pipeline_object.resample(image_models)
        del result_model, image_models, asn_data
        gc.collect()
    else:
        pipeline_object.log.info('Step outlier_detection is skipped because all output files exist: {}'.format(repr(processed_image_files)))
    # 
    # 5. source_catalog
    pipeline_object.source_catalog.output_file = output_name
    pipeline_object.source_catalog.search_output_file = False
    pipeline_object.source_catalog.save_results = True
    processing_image_files = processed_image_files
    processed_image_files = [output_name+'_cat.ecsv']
    pipeline_object.log.info('Checking source_catalog output file existence: {}'.format(repr(processed_image_files)))
    if not np.all(list(map(os.path.exists, processed_image_files))):
        asn_from_list_to_file(processing_image_files, 'asn_source_catalog.json')
        asdf_from_step_to_file(pipeline_object.source_catalog, 'asdf_source_catalog.txt')
        image_models = datamodels.ModelContainer()
        asn_data = image_models.read_asn('asn_source_catalog.json')
        image_models.from_asn(asn_data) # see "tweakreg_step.py"
        result_catalog = pipeline_object.source_catalog(image_models[0])
        del result_catalog, image_models, asn_data
        gc.collect()
    else:
        pipeline_object.log.info('Step source_catalog is skipped because all output files exist: {}'.format(repr(processed_image_files)))
    # 
    # save pars as asdf file
    asdf_filepath = output_name + '_calwebb_image3.asdf'
    if os.path.isfile(asdf_filepath):
        shutil.move(asdf_filepath, asdf_filepath+'.backup')
    asdf_object = AsdfFile(pipeline_object.get_pars())
    asdf_object.write_to(asdf_filepath)
    pipeline_object.log.info('Parameters are saved into {}'.format(asdf_filepath))



def join_a_string_list_with_omitted_parts(str_list, sep='+', omit='many', maxN=4):
    if len(str_list) <= maxN:
        return sep.join(str_list)
    else:
        return sep.join(str_list[0:maxN-1])+sep+omit+sep+str_list[-1]






# Defaults
# see -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
DEFAULT_KERNEL = 'square'
DEFAULT_PIXFRAC = 1.0 # 0.5
DEFAULT_PIXEL_SCALE_RATIO = 0.48
DEFAULT_PIXEL_SCALE = None



# Main
@click.command()
@click.argument('input_cal_files', nargs=-1, type=str)
@click.argument('output_dir', type=click.Path(exists=False))
@click.option('--program', 'select_program', type=str, multiple=True, default=None, help='Specify a program.')
@click.option('--obsnum', 'select_obsnum', type=str, multiple=True, default=None, help='Specify a obsnum.')
@click.option('--visitnum', 'select_visitnum', type=str, multiple=True, default=None, help='Specify a visitnum.')
@click.option('--instrument', 'select_instrument', type=str, multiple=True, default=None, help='Specify an instrument.')
@click.option('--filter', 'select_filter', type=str, multiple=True, default=None, help='Specify a filter.')
@click.option('--kernel', type=click.Choice(['point', 'square', 'gaussian', 'tophat', 'turbo', 'lanczos2', 'lanczos3']), 
                          default=DEFAULT_KERNEL, 
                          help='Drizzle kernel.'+
                               'The form of the kernel function used to distribute flux onto the output image. '+
                               'Available kernels are square, gaussian, point, tophat, turbo, lanczos2, and lanczos3.')
@click.option('--pixfrac', type=float, 
                           default=DEFAULT_PIXFRAC, 
                           help='Drizzle pixfrac.')
@click.option('--pixel-scale-ratio', type=float, 
                                     default=DEFAULT_PIXEL_SCALE_RATIO, 
                                     help='Drizzle pixel scale ratio. Ratio of input to output pixel scale. '+
                                          'A value of 0.5 means the output image would have 4 pixels sampling '+
                                          'each input pixel. Ignored when --pixel-scale are provided.')
@click.option('--pixel-scale', type=float, 
                               default=DEFAULT_PIXEL_SCALE, 
                               help='Drizzle pixel scale. '+
                                    'Absolute pixel scale in arcsec. '+
                                    'When provided, overrides pixel_scale_ratio.')
@click.option('--abs-refcat', type=click.Path(exists=True), 
                              default=None, 
                              help='Absolute reference catalog, must contain `RA` and `DEC` columns, optionally `weight`. See `tweakwcs/imalign.py` `align_wcs`.')
@click.option('--save-info-table-dir', type=click.Path(exists=False), 
                                       default=None, 
                                       help='Save the dataset-grouped info table to disk. Default directory is the `output_dir`.')
@click.option('--save-info-table-name', type=str, 
                                        default='mosaic_info_table', 
                                        help='Save the dataset-grouped info table to disk. Default file name is "info_table" and two formats are saved, "csv" and "txt".')
@click.option('--use-custom-catalogs', type=bool, 
                                       default=True, 
                                       help='Set use_custom_catalogs to True. Only valid if each "*_cal.fits" has a catalog named "*_cal_cat_for_tweakreg.csv" in the same directory.')
@click.option('--enforce-user-order', type=bool, 
                                      default=False, 
                                      help='Set enforce_user_order to True for tweakreg_step. Align images in user specified order.')
@click.option('--grid-step', type=float, 
                             default=3.0, 
                             help='If `--very-big-mosaic` is set, this will be the box size to divide the full mosaic into. In arcminutes.')
@click.option('--combine-program/--no-combine-program', is_flag=True, default=False, help='Combine all programs into one.')
@click.option('--combine-obsnum/--no-combine-obsnum', is_flag=True, default=False, help='Combine all obsnum into one.')
@click.option('--combine-visitnum/--no-combine-visitnum', is_flag=True, default=False, help='Combine all visitnum into one. If `--combine-obsnum` is set then this will also be set to True.')
@click.option('--combine-filter/--no-combine-filter', is_flag=True, default=False, help='Combine all filters into one.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False, help='Overwrite?')
@click.option('--run-individual-steps/--no-run-individual-steps', is_flag=True, default=False, help='Run individual step of JWST stage3 pipeline? This is turned on if abs_refcat is provided!')
@click.option('--very-big-mosaic/--no-very-big-mosaic', is_flag=True, default=False, help='Very big mosaic mode! If Ture, we will divide images into boxes with size set by `grid_step` in arcmin.')
def main(
        input_cal_files, 
        output_dir, 
        select_program, 
        select_visitnum, 
        select_obsnum, 
        select_instrument, 
        select_filter, 
        kernel, 
        pixfrac, 
        pixel_scale_ratio, 
        pixel_scale, 
        abs_refcat, 
        save_info_table_dir, 
        save_info_table_name, 
        use_custom_catalogs, 
        enforce_user_order, 
        grid_step, 
        combine_program, 
        combine_obsnum, 
        combine_visitnum, 
        combine_filter, 
        overwrite, 
        run_individual_steps,
        very_big_mosaic, 
    ):
    
    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Print JWST pipeline version
    logger.info('JWST pipeline version: {}'.format(jwst.__version__))
    
    # Check CRDS 
    try:
        logger.info("CRDS_PATH: {}".format(os.environ['CRDS_PATH']))
    except KeyError:
        logger.error("Error! CRDS_PATH environment variable not set!")
        sys.exit(-1)
        
    try:
        logger.info("CRDS_SERVER_URL: {}".format(os.environ['CRDS_SERVER_URL']))
    except KeyError:
        logger.error("Error! CRDS_SERVER_URL environment variable not set!")
        sys.exit(-1)
        
    try:
        logger.info("CRDS_CONTEXT: {}".format(os.environ['CRDS_CONTEXT']))
    except KeyError:
        logger.error("Error! CRDS_CONTEXT environment variable not set!")
        sys.exit(-1)
    
    
    # Check output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    
    # expand input file list if there is a wildcard
    input_files = []
    for input_file in input_cal_files:
        if input_file.find('*')>=0:
            for found_file in glob.glob(input_file):
                input_files.append(found_file)
        else:
            if not os.path.isfile(input_file):
                raise Exception('Error! File does not exist: {!r}'.format(input_file))
            input_files.append(input_file)
    
    
    # sort input files
    input_files.sort() # key=LooseVersion does not work
    
    
    # get info_list filtered by program, instrument, filter
    info_list = []
    for input_filepath in input_files:
        # read fits header
        header = fits.getheader(input_filepath, 0)
        # store into info_dict
        info_dict = OrderedDict()
        info_dict['program'] = header['PROGRAM'].strip()
        info_dict['obs_num'] = header['OBSERVTN'].strip()
        if 'TARGPROP' in header:
            info_dict['target_group'] = header['TARGPROP'].strip()
        if 'TARGNAME' in header:
            info_dict['target_name'] = '"'+header['TARGNAME'].strip()+'"'
        if 'TARG_RA' in header:
            info_dict['target_RA'] = header['TARG_RA']
        if 'TARG_DEC' in header:
            info_dict['target_Dec'] = header['TARG_DEC']
        info_dict['instrument'] = header['INSTRUME'].strip()
        info_dict['visit_num'] = header['VISIT'].strip()
        info_dict['visit_group'] = header['VISITGRP'].strip()
        info_dict['seq_id'] = header['SEQ_ID'].strip()
        info_dict['act_id'] = header['ACT_ID'].strip()
        info_dict['exposure'] = header['EXPOSURE'].strip()
        info_dict['filter'] = header['FILTER'].strip()
        if 'PUPIL' in header:
            info_dict['pupil'] = header['PUPIL'].strip()
        else:
            info_dict['pupil'] = '""'
        info_dict['file_path'] = input_filepath
        # check selecting instrument and filter
        if len(select_program) > 0:
            select_program = expand_list_with_comma_items(select_program)
            if info_dict['program'] not in select_program:
                continue
        if len(select_obsnum) > 0:
            select_obsnum = expand_list_with_comma_items(select_obsnum)
            if info_dict['obs_num'] not in select_obsnum:
                continue
        if len(select_visitnum) > 0:
            select_visitnum = expand_list_with_comma_items(select_visitnum)
            if info_dict['visit_num'] not in select_visitnum:
                continue
        if len(select_instrument) > 0:
            select_instrument = expand_list_with_comma_items(select_instrument)
            if info_dict['instrument'] not in select_instrument:
                continue
        if len(select_filter) > 0:
            select_filter = expand_list_with_comma_items(select_filter)
            if info_dict['filter'] not in select_filter:
                continue
        info_list.append(info_dict)
    
    
    # check len(info_list)
    if len(info_list) == 0:
        if len(select_program) > 0 or \
           len(select_obsnum) > 0 or \
           len(select_visitnum) > 0 or \
           len(select_instrument) > 0 or \
           len(select_filter) > 0:
            raise Exception('Error! No input file matching the specified ' + \
                'program {!r} obsnum {!r} visitnum {!r} instrument {!r} and filter {!r}'.format(
                select_program, select_obsnum, select_visitnum, select_instrument, select_filter))
        else:
            logger.error('Error! No input file to process!')
            sys.exit(255)
    
    
    # compile info table
    info_table_dict = OrderedDict()
    for i, info_dict in enumerate(info_list):
        if i == 0:
            for key in info_dict:
                info_table_dict[key] = []
        for key in info_dict:
            info_table_dict[key].append(info_dict[key])
    info_table = Table(info_table_dict)
    
    
    # save info table to disk
    if save_info_table_dir is None:
        save_info_table_dir = output_dir
    if not os.path.isdir(save_info_table_dir):
        os.makedirs(save_info_table_dir)
    info_table_file = os.path.join(save_info_table_dir, save_info_table_name)
    if os.path.isfile(f'{info_table_file}.txt'):
        shutil.move(f'{info_table_file}.txt', f'{info_table_file}.txt.backup')
    if os.path.isfile(f'{info_table_file}.csv'):
        shutil.move(f'{info_table_file}.csv', f'{info_table_file}.csv.backup')
    info_table.write(f'{info_table_file}.txt', format='ascii.fixed_width', delimiter=' ', bookend=True)
    with open(f'{info_table_file}.txt', 'r+') as fp:
        fp.seek(0)
        fp.write('#')
    info_table.write(f'{info_table_file}.csv', format='csv')
    
    
    # group by '{instrument}_{filter}' # '{program}_{obs_num}_{instrument}_{filter}'
    #group_by_columns = ['instrument', 'filter', 'program', 'obs_num']
    group_by_columns = ['instrument', 'filter', 'program', 'obs_num', 'visit_num']
    if combine_filter:
        if len(np.unique(info_table['filter'])) > 1:
            logger.warning('We are combining different filters!')
        group_by_columns.remove('filter')
    if combine_program:
        if len(np.unique(info_table['program'])) > 1:
            logger.warning('We are combining different programs!')
        group_by_columns.remove('program')
    if combine_obsnum:
        if len(np.unique(info_table['obs_num'])) > 1:
            logger.warning('We are combining different obs_num!')
        group_by_columns.remove('obs_num')
    if combine_obsnum:
        combine_visitnum = True # if `--combine-obsnum` is set, this will also be set to True.
    if combine_visitnum:
        if len(np.unique(info_table['visit_num'])) > 1:
            logger.warning('We are combining different visit_num!')
        group_by_columns.remove('visit_num')
    logger.info('group_by_columns: {}'.format(group_by_columns))
    groupped_table = info_table.group_by(group_by_columns)
    
    
    # prepare association file to process all rate files into one single output file
    output_lookup_dict = OrderedDict()
    for subgroup_key, subgroup_table in zip(groupped_table.groups.keys, groupped_table.groups):
        
        unique_programs = np.unique(subgroup_table['program'])
        groupped_program_str = join_a_string_list_with_omitted_parts(unique_programs) # = '+'.join(unique_programs)
        
        unique_obsnums = np.unique(subgroup_table['obs_num'])
        groupped_obsnum_str = join_a_string_list_with_omitted_parts(unique_obsnums) # '+'.join(unique_obsnums)
        
        unique_visitnums = np.unique(subgroup_table['visit_num'])
        groupped_visitnum_str = join_a_string_list_with_omitted_parts(unique_visitnums) # = '+'.join(unique_visitnums)
        
        unique_instruments = np.unique(subgroup_table['instrument'])
        groupped_instrument_str = join_a_string_list_with_omitted_parts(unique_instruments) # = '+'.join(unique_instruments)
        
        unique_filters = np.unique(subgroup_table['filter'])
        groupped_filter_str = join_a_string_list_with_omitted_parts(unique_filters) # = '+'.join(unique_filters)
        
        #unique_obsnums = np.unique(subgroup_table['obs_num'])
        #groupped_obsnum_str = join_a_string_list_with_omitted_parts(unique_obsnums) # = '+'.join(unique_obsnums)
        
        #subgroup_files = subgroup_table['file_path'].data.tolist()
        subgroup_files = [os.path.abspath(t) for t in subgroup_table['file_path']]
        
        # set output directory name
        output_name = 'jw{}_obs{}_{}_{}'.format(
            groupped_program_str, 
            groupped_obsnum_str, 
            groupped_instrument_str,
            groupped_filter_str,
        )
        
        groupped_asn_str = groupped_obsnum_str
        
        # check if this obs has multiple visits or not, we include visit num in the output directory name
        # if there are multiple visits for this single obs. 
        if len(unique_obsnums) == 1:
            if (not combine_obsnum) and (len(np.unique(info_table['visit_num'])) > len(unique_visitnums)):
                
                output_name = 'jw{}_obs{}_visit{}_{}_{}'.format(
                    groupped_program_str,
                    groupped_obsnum_str,
                    groupped_visitnum_str,
                    groupped_instrument_str,
                    groupped_filter_str,
                )
                
                groupped_asn_str = groupped_obsnum_str + groupped_visitnum_str
        
        
        # set output full path
        output_subdir = os.path.join(output_dir, output_name)
        output_file = output_name + '_i2d.fits' # dose not contain the path
        #if very_big_mosaic:
        #    output_file = output_name + '_very_big_mosaic.fits'
        output_filepath = os.path.join(output_subdir, output_file)
        
        # add to lookup dict
        output_lookup_dict[output_filepath] = subgroup_files
        
        # check existing output file, skip the processing or backup and overwrite it. 
        if not os.path.isdir(output_subdir):
            os.makedirs(output_subdir)
        else:
            if os.path.isfile(output_filepath):
                if not overwrite:
                    logger.warning(f"Found processed data file {output_filepath!r} and overwrite is set to False. " + 
                                    "We will skip the processing.")
                    continue
        
        # check if another script is processing this data
        if os.path.exists(output_subdir + '.touch'):
            logger.warning(f"WARNING WARNING WARNING! Found another script processing data dir {output_subdir!r} (with a *.touch file). " + 
                            "We will skip processing it for now!")
            continue
        
        # create touch file
        with open(output_subdir + '.touch', 'a'):
            os.utime(output_subdir + '.touch', None)
        
        # backup existing file
        if os.path.isfile(output_filepath):
            shutil.move(output_filepath, output_filepath+'.backup')
        
        # check existing data set catalog files
        # if there is a '{dataset_cal}_cat_for_tweakreg.csv' file, 
        # use it as catfile for tweakreg
        imcat_dict = OrderedDict() # see "jwst/tweakreg/tweakreg_step.py"
        if use_custom_catalogs:
            for subgroup_file in subgroup_files:
                imcat_file = os.path.splitext(subgroup_file)[0]+'_cat_for_tweakreg.csv'
                if os.path.isfile(imcat_file):
                    imcat_dict[subgroup_file] = imcat_file # all are abspath
        
        # # prepare the json-format "association" file
        # asn_dict = OrderedDict()
        # asn_dict['asn_type'] = 'None'
        # asn_dict['asn_rule'] = 'DMS_Level3_Base'
        # asn_dict['version_id'] = None
        # asn_dict['code_version'] = jwst.__version__
        # asn_dict['degraded_status'] = 'No known degraded exposures in association.'
        # asn_dict['program'] = groupped_program_str
        # asn_dict['constraints'] = 'No constraints'
        # asn_dict['asn_id'] = groupped_asn_str
        # asn_dict['asn_pool'] = 'none'
        # asn_dict['products'] = []
        # product_dict = OrderedDict()
        # product_dict['name'] = output_name
        # product_dict['members'] = []
        # asn_dict['products'].append(product_dict)
        # for subgroup_file in subgroup_files:
        #     # get 'expname'
        #     # relative to asn file dir path, 
        #     # see "jwst/datamodels/container.py"
        #     # see also "jwst/tweakreg/tweakreg_step.py", 
        #     # ```
        #     # filename = member['expname']
        #     # member['expname'] = path.join(asn_dir, filename)
        #     # ```
        #     #expname = os.path.relpath(subgroup_file, os.path.dirname(output_subdir))
        #     # 
        #     # 20221203: now using symlink and run the processing from the output_subdir!
        #     expname = os.path.basename(subgroup_file)
        #     explink = os.path.join(output_subdir, expname)
        #     product_dict['members'].append(
        #         {'expname': expname, 
        #          'exptype': 'science'
        #         }
        #     )
        #     if not os.path.exists(explink):
        #         if os.path.islink(explink):
        #             os.remove(explink)
        #         os.symlink(os.path.relpath(subgroup_file, output_subdir), 
        #                    explink)
        #     # 
        #     # if there is a '{dataset_cal}_cat_for_tweakreg.csv' file, 
        #     # use it as catfile for tweakreg
        #     calcat_file = os.path.splitext(subgroup_file)[0]+'_cat_for_tweakreg.csv'
        #     if os.path.isfile(calcat_file):
        #         calcat_name = os.path.basename(calcat_file)
        #         calcat_link = os.path.join(output_subdir, calcat_name)
        #         if not os.path.exists(calcat_link):
        #             if os.path.islink(calcat_link):
        #                 os.remove(calcat_link)
        #             os.symlink(os.path.relpath(calcat_file, output_subdir), 
        #                        calcat_link)
        #         catkey = expname
        #         catdict[catkey] = calcat_name
        
        # asn_filename = 'asn.json'
        # asn_file = os.path.join(output_subdir, asn_filename)
        
        # if os.path.isfile(asn_file):
        #     shutil.move(asn_file, asn_file+'.backup')
        
        # with open(asn_file, 'w') as fp:
        #     json.dump(asn_dict, fp, indent=4)
        
        # catfile = None
        # if use_custom_catalogs and len(catdict) > 0:
        #     catfilename = 'catfile.txt'
        #     catfilepath = os.path.join(output_subdir, catfilename)
        #     if os.path.isfile(catfilepath):
        #         shutil.move(catfilepath, catfilepath+'.backup')
        #     with open(catfilepath, 'w') as fp:
        #         for catkey in catdict:
        #             fp.write('{} {}\n'.format(catkey, catdict[catkey]))
        #     catfile = catfilename # we will run the process inside the output_subdir, so set only filename here
        
        
        if abs_refcat is not None and abs_refcat != '':
            abs_refcat = os.path.abspath(abs_refcat) # get abspath
            run_individual_steps = True
        
        
        # Print progress
        logger.info("Processing {} -> {}".format(subgroup_files, output_filepath))
        
        
        # chdir
        current_dir = os.getcwd()
        logger.info("chdir {}".format(output_subdir))
        os.chdir(output_subdir)
        
        
        # prepare "asn.json" under 'output_subdir'
        asn_filename = asn_from_list_to_file(subgroup_files, 'asn.json')
        
        
        # prepare "catfile.txt" under 'output_subdir'
        catfile = imcat_dict_to_catfile(imcat_dict, subgroup_files, 'catfile.txt') # this will produce nothing if imcat_dict is empty.
        
        
        # prepare to run
        pipeline_object = calwebb_image3.Image3Pipeline()
        
        
        # Set some parameters that pertain to the entire pipeline
        pipeline_object.output_file = output_name
        pipeline_object.output_ext = ".fits"
        pipeline_object.save_results = True
        
        
        # Set some parameters that pertain to some of the individual steps
        pipeline_object.tweakreg.enforce_user_order = enforce_user_order
        pipeline_object.tweakreg.minobj = 7 # default is 15
        pipeline_object.tweakreg.searchrad = 1.0 # default is 2.0
        pipeline_object.tweakreg.separation = 1.0 # default is 0.1 arcsec
        pipeline_object.tweakreg.tolerance = 0.7 #<20230318># 0.7 # 1.0 # default is 0.7 arcsec
        pipeline_object.tweakreg.save_catalogs = True # will save as "./*_cal_cat.ecsv", but always in current directory, see "tweakreg_step.py"
        pipeline_object.tweakreg.save_results = True # will save as "{output_subdir}/*_cal_tweakreg.fits"
        pipeline_object.tweakreg.search_output_file = False # do not use output_file define in parent step
        pipeline_object.tweakreg.output_use_model = True # use DataModel.meta.filename as output_file
        # set catfile if found
        if catfile is not None and catfile != '':
            pipeline_object.tweakreg.use_custom_catalogs = use_custom_catalogs
            pipeline_object.tweakreg.catfile = catfile
        # set abs_refcat if user has input that
        if abs_refcat is not None and abs_refcat != '':
            pipeline_object.tweakreg.abs_refcat = abs_refcat
            #pipeline_object.abs_fitgeometry = 'rshift' # 'shift', 'rshift', 'rscale', 'general' (shift, rotation, and scale)
            pipeline_object.tweakreg.abs_minobj = 1 # default is 15
            pipeline_object.tweakreg.abs_searchrad = 1.0 # default is 6.0
            pipeline_object.tweakreg.abs_separation = 1.0 # default is 0.1 arcsec
            pipeline_object.tweakreg.abs_tolerance = 1.0 # default is 0.7 arcsec
            #pipeline_object.tweakreg.abs_use2dhist = False # default is True, but ...
            pipeline_object.tweakreg.abs_use2dhist = True # default is True, but ...
            #pipeline_object.tweakreg.abs_fitgeometry = 'shift' # 'rshift'
            pipeline_object.tweakreg.save_abs_catalog = True
            #pipeline_object.tweakreg.save_abs_catalog = True # only if abs_refcat `gaia_cat_name in SINGLE_GROUP_REFCAT`
        # do manual 2dhist for a better tweakreg
        if catfile is not None and catfile != '' and \
           abs_refcat is not None and abs_refcat != '':
            from util_match_cat_file_to_abs_refcat_with_2dhist import match_cat_file_to_abs_refcat_with_2dhist
            pipeline_object.tweakreg.catfile, pipeline_object.tweakreg.abs_refcat = \
                match_cat_file_to_abs_refcat_with_2dhist(catfile, abs_refcat,
                    output_dir = os.getcwd())
        
        
        # SkyMatchStep
        pipeline_object.skymatch.subtract = True
        pipeline_object.skymatch.skymethod = "global+match"
        #pipeline_object.skymatch.lsigma = 2.0
        #pipeline_object.skymatch.usigma = 2.0
        #pipeline_object.skymatch.nclip = 10
        #pipeline_object.skymatch.upper = 1.0
        pipeline_object.skymatch.save_results = True # will save as '*_skymatch.fits'
        
        
        # OutlierDetectionStep
        pipeline_object.outlier_detection.pixfrac = pixfrac
        pipeline_object.outlier_detection.in_memory = True # must be True otherwise outlier_detection flag_cr sci_image.dq updates are not saved to the same object! 20230307
        #pipeline_object.outlier_detection.suffix = 'crf' # just use the default
        #pipeline_object.outlier_detection.save_results = True # will save as '*_crf.fits'
        
        # 20230307
        pipeline_object.outlier_detection.snr = '3.0 2.5'
        pipeline_object.outlier_detection.scale = '10. 5'
        
        
        # ResampleStep
        # Set the ratio of input to output pixels to create an output mosaic 
        # on a 0.015"/pixel scale
        # see -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
        #pipeline_object.resample_data = True # True is the default
        pipeline_object.resample.kernel = kernel
        pipeline_object.resample.pixfrac = pixfrac
        pipeline_object.resample.pixel_scale_ratio = pixel_scale_ratio
        pipeline_object.resample.pixel_scale = pixel_scale
        #pipeline_object.resample.suffix = 'i2d' # just use the default
        #pipeline_object.resample.save_results = True
        
        
        # SourceCatalogStep
        #pipeline_object.source_catalog.save_results = True
        
        
        # run
        if not run_individual_steps:
            
            asdf_from_step_to_file(pipeline_object, 'asdf_calwebb_image3.asdf')
            
            pipeline_object.run(asn_filename)
        
        # we can also run individual steps following "jwst/pipeline/calwebb_image3.py" `process()`
        else: # TODO
            
            # 20230109 very big mosaic
            
            if very_big_mosaic:
                pipeline_object.log.info('*-*-*\n*-*-* Very big mosaic mode! *-*-*\n*-*-*')
                
                # get pixel_size if only pixel_scale_ratio is given
                # and pa_v3 # no need 
                with datamodels.open(subgroup_files[0]) as refmodel:
                    ref_pa_v3 = refmodel.meta.pointing.pa_v3 # no need 
                    if pipeline_object.resample.pixel_scale is None:
                        # see "jwst/assign_wcs/util.py"
                        from jwst.assign_wcs.util import compute_scale
                        ref_fiducial = np.array([
                            refmodel.meta.wcsinfo.ra_ref, 
                            refmodel.meta.wcsinfo.dec_ref
                        ])
                        ref_pixel_size = compute_scale(
                            refmodel.meta.wcs, 
                            ref_fiducial,
                            pscale_ratio = pipeline_object.resample.pixel_scale_ratio
                        ) * 3600.0
                    else:
                        ref_pixel_size = pipeline_object.resample.pixel_scale
                
                # make mosaic subgrouping
                from util_run_mosaic_image_subgrouping import run_mosaic_image_subgrouping
                mosaic_group_meta, mosaic_group_table, group_asn_files, group_cat_files = \
                    run_mosaic_image_subgrouping(
                        image_files = None, # we input asn_file, no need to input image_files.
                        output_name = 'run_mosaic_image_subgroup', # hard-coded output prefix
                        asn_file = asn_filename, 
                        grid_step = grid_step, 
                        cat_file = catfile, # it cannot contain a path otherwise things will be complicated.
                        pixel_size = ref_pixel_size, 
                    )
                
                # loop each mosaic group in 'mosaic_group_dir'
                first_valid_idx = None
                for mosaic_group_idx in range(len(mosaic_group_table)):
                    
                    if mosaic_group_table['n_images'][mosaic_group_idx] == 0:
                        continue
                    
                    if first_valid_idx is None:
                        first_valid_idx = mosaic_group_idx
                    
                    mosaic_group_dir = mosaic_group_table['group_dir'][mosaic_group_idx]
                    mosaic_group_images = mosaic_group_table['image_files'][mosaic_group_idx].split(',')
                    mosaic_group_output_file = os.path.join(mosaic_group_dir, 
                                                            mosaic_group_dir+'_i2d.fits')
                    if os.path.isfile(mosaic_group_output_file) and not overwrite:
                        logger.warning(f"Found processed very big mosaic tile {mosaic_group_output_file!r} and overwrite is set to False. " + 
                                        "We will skip the processing.")
                        pass
                    else:
                        pipeline_object_copy = copy.deepcopy(pipeline_object)
                        pipeline_object_copy.tweakreg.catfile = os.path.basename(group_cat_files[mosaic_group_idx])
                        pipeline_object_copy.resample.crpix = (
                            mosaic_group_table['crpix1'][mosaic_group_idx],
                            mosaic_group_table['crpix2'][mosaic_group_idx]
                        ) # IMPORTANT: resample crpix is 0-based according to ...
                        pipeline_object_copy.resample.crval = (
                            mosaic_group_table['crval1'][mosaic_group_idx],
                            mosaic_group_table['crval2'][mosaic_group_idx]
                        )
                        pipeline_object_copy.resample.rotation = 0.0 # A value of 0.0 would orient the final output image to be North up. 
                        #pipeline_object_copy.resample.rotation = ref_pa_v3-360. # A value of 0.0 would orient the final output image to be North up. 
                        pipeline_object_copy.resample.output_shape = (
                            int(mosaic_group_table['naxis1'][mosaic_group_idx]),
                            int(mosaic_group_table['naxis2'][mosaic_group_idx])
                        )
                        pipeline_object_copy.resample.pixel_scale = ref_pixel_size
                        # 
                        pipeline_object_copy.log.info(
                            '*-*-*\n*-*-* Processing very big mosaic tile {!r} ({} x {} @ {:.3f} mas) *-*-*\n*-*-*'.format(
                            mosaic_group_dir, 
                            pipeline_object_copy.resample.output_shape[0], 
                            pipeline_object_copy.resample.output_shape[1], 
                            ref_pixel_size*1000.0)
                        )
                        mosaic_parent_dir = os.getcwd()
                        os.chdir(mosaic_group_dir)
                        run_individual_steps_for_image_files(
                            pipeline_object_copy, 
                            mosaic_group_images,
                            mosaic_group_dir, # use this as output file name
                        )
                        os.chdir(mosaic_parent_dir)
                        pipeline_object_copy.log.info(
                            '*-*-*\n*-*-* Processed very big mosaic tile {!r} ({} x {} @ {:.3f} mas) *-*-*\n*-*-*'.format(
                            mosaic_group_dir, 
                            pipeline_object_copy.resample.output_shape[0], 
                            pipeline_object_copy.resample.output_shape[1], 
                            ref_pixel_size*1000.0)
                        )
                        # 
                        del pipeline_object_copy
                
                # 
                gc.collect()
                
                # initialize the output very big mosaic fits file using memmap
                pipeline_object.log.info('Building final very big mosaic {!r}'.format(output_file))
                if not os.path.isfile(output_file) or overwrite:
                    # read first valid tile image
                    imagefile1 = os.path.join(mosaic_group_table['group_dir'][first_valid_idx], 
                                              mosaic_group_table['group_dir'][first_valid_idx]+'_i2d.fits')
                    header1 = fits.getheader(imagefile1, 0)
                    # initialize with a small data array
                    data = np.zeros((100, 100), dtype=np.float32)
                    hdu = fits.PrimaryHDU(data=data, header=header1.copy())
                    header = hdu.header
                    header['BITPIX'] = -32
                    header['NAXIS1'] = mosaic_group_meta['naxis1']
                    header['NAXIS2'] = mosaic_group_meta['naxis2']
                    header.append('', useblanks=False, end=True)
                    header.append(('', 'Mosaic WCS Information'), useblanks=False, end=True)
                    header.append('', useblanks=False, end=True)
                    header.append(('RADESYS', 'FK5'), useblanks=False, end=True)
                    header.append(('EQUINOX', 2000.0), useblanks=False, end=True)
                    header.append(('CTYPE1', 'RA---TAN'), useblanks=False, end=True)
                    header.append(('CTYPE2', 'DEC--TAN'), useblanks=False, end=True)
                    header.append(('CUNIT1', 'deg'), useblanks=False, end=True)
                    header.append(('CUNIT2', 'deg'), useblanks=False, end=True)
                    header.append(('CRVAL1', mosaic_group_meta['crval1'], 'In degrees.'), useblanks=False, end=True)
                    header.append(('CRVAL2', mosaic_group_meta['crval2'], 'In degrees.'), useblanks=False, end=True)
                    header.append(('CRPIX1', mosaic_group_meta['crpix1']), useblanks=False, end=True)
                    header.append(('CRPIX2', mosaic_group_meta['crpix2']), useblanks=False, end=True)
                    header.append(('CDELT1', mosaic_group_meta['cdelt1'], 'In degrees.'), useblanks=False, end=True)
                    header.append(('CDELT2', mosaic_group_meta['cdelt2'], 'In degrees.'), useblanks=False, end=True)
                    header.append(('CROTA2', mosaic_group_meta['crota2'], 'In degrees.'), useblanks=False, end=True)
                    header.append(('PIXSCALE', ref_pixel_size, 'In arcsec.'), useblanks=False, end=True)
                    header.append('', useblanks=False, end=True)
                    header.append(('HISTORY', 'Mosaic made with Crab.Toolkit.JWST/bin/go-jwst-imaging-stage-3-step-1.py'), useblanks=False, end=True)
                    header.append('', useblanks=False, end=True)
                    header['EXTEND'] = True
                    header['EXTNAME'] = 'SCI'
                    #while (len(header))%(2880//80) != 0: # fits header block is padded to N*2880 by standard
                    #    header.append('', useblanks=False, end=True)
                    # 
                    header.tofile(output_file, endcard=True, padding=True, overwrite=overwrite)
                    headerstr = header.tostring(sep='', endcard=True, padding=True)
                    with open(output_file, 'rb+') as fobj:
                        last_byte = len(headerstr) + (header['NAXIS1'] * header['NAXIS2'] * np.abs(header['BITPIX']//8))
                        last_byte_padded = int(np.ceil(float(last_byte)/2880))*2880 # fits data blocks are padded to 2880 by standard
                        fobj.seek(last_byte-1)
                        for iter_byte in range(last_byte-1, last_byte_padded):
                            fobj.write(b'\0')
                        fpos = fobj.tell()
                        # 
                        # write extension 2
                        header2 = fits.Header()
                        header2['XTENSION'] = 'IMAGE'
                        header2['BITPIX'] = -32
                        header2['NAXIS1'] = mosaic_group_meta['naxis1']
                        header2['NAXIS2'] = mosaic_group_meta['naxis2']
                        for key in ['NAXIS', 'NAXIS1', 'NAXIS2', 'RADESYS', 'EQUINOX', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2', 
                                    'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CROTA2']:
                            header2[key] = header[key]
                        header2['PCOUNT'] = 0
                        header2['GCOUNT'] = 1
                        header2['EXTNAME'] = 'ERR'
                        #while (len(header2))%(2880//80) != 0: # fits header block is padded to N*2880 by standard
                        #    header2.append('', useblanks=False, end=True)
                        # 
                        header2str = header2.tostring(sep='', endcard=True, padding=True) # manually control the padding to Nx2880
                        fobj.write(header2str.encode())
                        last_byte = len(header2str) + (header2['NAXIS1'] * header2['NAXIS2'] * np.abs(header2['BITPIX']//8))
                        last_byte_padded = int(np.ceil(float(last_byte)/2880))*2880 # fits data blocks are padded to 2880 by standard
                        fobj.seek(fpos + last_byte-1)
                        for iter_byte in range(last_byte-1, last_byte_padded):
                            fobj.write(b'\0')
                    # 
                    del header2
                    del header1
                    del header
                    del data
                    del hdu
                    
                with fits.open(output_file, mode='update', memmap=True) as out_hdul:
                    for mosaic_group_idx in range(len(mosaic_group_table)):
                        if mosaic_group_table['n_images'][mosaic_group_idx] == 0:
                            continue
                        imagefile1 = os.path.join(mosaic_group_table['group_dir'][mosaic_group_idx], 
                                                  mosaic_group_table['group_dir'][mosaic_group_idx]+'_i2d.fits')
                        with fits.open(imagefile1, mode='readonly', memmap=True) as in_hdul:
                            in_x1 = mosaic_group_table['in_x1'][mosaic_group_idx]
                            in_x2 = mosaic_group_table['in_x2'][mosaic_group_idx]
                            in_y1 = mosaic_group_table['in_y1'][mosaic_group_idx]
                            in_y2 = mosaic_group_table['in_y2'][mosaic_group_idx]
                            out_x1 = mosaic_group_table['out_x1'][mosaic_group_idx]
                            out_x2 = mosaic_group_table['out_x2'][mosaic_group_idx]
                            out_y1 = mosaic_group_table['out_y1'][mosaic_group_idx]
                            out_y2 = mosaic_group_table['out_y2'][mosaic_group_idx]
                            pipeline_object.log.info('Copying image data from {!r} [{}:{},{}:{}] to {!r} [{}:{},{}:{}]'.format(
                                imagefile1, in_y1, in_y2, in_x1, in_x2, 
                                output_file, out_y1, out_y2, out_x1, out_x2, 
                            ))
                            #out_hdul[0].data[out_y1:out_y2,out_x1:out_x2] = in_hdul['SCI'].data[in_y1:in_y2,in_x1:in_x2]
                            #<TODO><20230116># also copy 'ERR' 'CON'
                            
                            out_hdul['SCI'].data[out_y1:out_y2,out_x1:out_x2] = in_hdul['SCI'].data[in_y1:in_y2,in_x1:in_x2]
                            out_hdul['ERR'].data[out_y1:out_y2,out_x1:out_x2] = in_hdul['ERR'].data[in_y1:in_y2,in_x1:in_x2]
                            
                    out_hdul.flush()
                
                # 
                gc.collect()
                
            else:
                
                # 20230412 still has outlier detection pixel issue, already set in_memory = True
                #          it must be that calling `` does not return the updated image models.
                # run_individual_steps_for_one_asn_file(
                #     pipeline_object, 
                #     asn_filename,
                #     output_name, 
                # )
                
                with open(asn_filename, 'r') as fp:
                    asn_dict_tmp = load_asn(fp)
                image_files_tmp = [t['expname'] for t in asn_dict_tmp['products'][0]['members'] if t['exptype']=='science']
                run_individual_steps_for_image_files(
                    pipeline_object, 
                    image_files_tmp,
                    output_name, 
                )
            
            # 
            # 2023-01-09 -- replaced by `run_individual_steps_for_one_asn_file()`
            # 
            # with datamodels.open(asn_filename, asn_exptypes=['science']) as input_models:
            #     # 
            #     # 1. tweakreg
            #     pipeline_object.tweakreg.output_file = output_name
            #     pipeline_object.tweakreg.save_results = True
            #     pipeline_object.tweakreg.search_output_file = False
            #     image_models = pipeline_object.tweakreg(input_models)
            #     # TODO: sometimes tweakreg update_fits_wcsinfo fails because `max_pix_error` is too small, 
            #     # we can run this manually with a larger `max_pix_error`, then save_model.
            #     from jwst.assign_wcs.util import update_fits_wcsinfo
            #     for image_model in image_models:
            #         update_fits_wcsinfo(image_model,
            #             max_pix_error=1.,
            #             npoints=16)
            #     pipeline_object.tweakreg.save_model(image_models, idx=None, suffix='tweakreg',
            #                     format=pipeline_object.tweakreg.name_format, force=True)
            #     # 
            #     # 2. skymatch bkgmatch
            #     pipeline_object.skymatch.output_file = output_name
            #     pipeline_object.skymatch.save_results = True
            #     pipeline_object.skymatch.search_output_file = False
            #     image_models = pipeline_object.skymatch(image_models)
            #     # 
            #     # 3. outlier detection
            #     # here we can split image_models into each obsnum group, 
            #     # then parallelize the outlier_detection process,
            #     # this also saves memory usage
            #     pipeline_object.outlier_detection.output_file = output_name
            #     pipeline_object.outlier_detection.save_results = True
            #     pipeline_object.outlier_detection.suffix = 'crf'
            #     pipeline_object.outlier_detection.search_output_file = False
            #     if len(unique_obsnums) == 1 and len(unique_visitnums) == 1:
            #         image_models = pipeline_object.outlier_detection(image_models)
            #     else:
            #         from util_run_outlier_detection_in_parallel import run_outlier_detection_in_parallel
            #         image_models = run_outlier_detection_in_parallel(pipeline_object, image_models)
            #     pipeline_object.outlier_detection.save_model(image_models, idx=None, suffix='outlierdetection',
            #                     format=pipeline_object.outlier_detection.name_format, force=True)
            #     # 
            #     # 4. resample
            #     pipeline_object.resample.output_file = output_name
            #     pipeline_object.resample.save_results = True
            #     pipeline_object.resample.suffix = 'i2d'
            #     result = pipeline_object.resample(image_models)
            #     # 
            #     # 5. source_catalog
            #     pipeline_object.source_catalog.output_file = output_name
            #     pipeline_object.source_catalog.search_output_file = False
            #     pipeline_object.source_catalog.save_results = True
            #     if isinstance(result, datamodels.ImageModel) and result.meta.cal_step.resample == 'COMPLETE':
            #         pipeline_object.source_catalog(result)
        
        
        # move some output files to subdirectories
        # for example the "outlier_i2d.fits", so that we do not find them when doing an `ls *_i2d.fits`
        for i in range(len(subgroup_files)):
            subgroup_file = subgroup_files[i]
            dataset_name = os.path.basename(subgroup_file).replace('_cal.fits', '')
            if os.path.isfile(f'{output_name}_{i}_outlier_i2d.fits'):
                if not os.path.isdir(f'{output_name}_{i}_outlier'):
                    os.makedirs(f'{output_name}_{i}_outlier')
                shutil.move(f'{output_name}_{i}_outlier_i2d.fits', 
                            f'{output_name}_{i}_outlier/{output_name}_{i}_outlier_i2d.fits')
        
        
        # chdir back
        logger.info("chdir {}".format(current_dir))
        os.chdir(current_dir)
        
        
        # remove touch file
        os.remove(output_subdir + '.touch')
        
        
        # check
        assert os.path.isfile(output_filepath)
        
        
        # clean up
        # for i in range(len(subgroup_files)):
        #     subgroup_file = subgroup_files[i]
        #     dataset_name = subgroup_file.replace('_cal.fits', '')
        #     if os.path.isfile(f'{output_name}_{i}_outlier_i2d.fits'):
        #         if not os.path.isdir(f'{output_subdir}/{output_name}_{i}_outlier'):
        #             os.makedirs(f'{output_subdir}/{output_name}_{i}_outlier')
        #         shutil.move(f'{output_name}_{i}_outlier_i2d.fits', 
        #                     f'{output_subdir}/{output_name}_{i}_outlier/{output_name}_{i}_outlier_i2d.fits')
        #     # 
        #     if os.path.isfile(f'{output_subdir}/{dataset_name}_tweakreg.fits'):
        #         if not os.path.isdir(f'{output_subdir}/{output_name}_{i}_tweakreg'):
        #             os.makedirs(f'{output_subdir}/{output_name}_{i}_tweakreg')
        #         shutil.move(f'{output_subdir}/{dataset_name}_tweakreg.fits', 
        #                     f'{output_subdir}/{output_name}_{i}_tweakreg/{dataset_name}_tweakreg.fits')
        #     # 
        #     if os.path.isfile(f'{dataset_name}_cal_cat.ecsv'):
        #         if not os.path.isdir(f'{output_subdir}/{output_name}_{i}_tweakreg'):
        #             os.makedirs(f'{output_subdir}/{output_name}_{i}_tweakreg')
        #         shutil.move(f'{dataset_name}_cal_cat.ecsv', 
        #                     f'{output_subdir}/{output_name}_{i}_tweakreg/{dataset_name}_cal_cat.ecsv')
        
        # log
        logger.info("Processed {} -> {}".format(subgroup_files, output_filepath))
    
    
    # update info table with final mosaic
    if os.path.isfile(f'{info_table_file}.json'):
        shutil.move(f'{info_table_file}.json', f'{info_table_file}.json.backup')
    with open(f'{info_table_file}.json', 'w') as fp:
        json.dump(output_lookup_dict, fp, indent=4)
    if os.path.isfile(f'{info_table_file}.list'):
        shutil.move(f'{info_table_file}.out', f'{info_table_file}.out.backup')
    with open(f'{info_table_file}.out', 'w') as fp:
        for key in output_lookup_dict.keys():
            fp.write(key+'\n')
    
    
    # log
    logger.info("All done!")




# Entry
if __name__ == '__main__':
    main()



