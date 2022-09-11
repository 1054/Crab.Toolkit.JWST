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
import os, sys, re, json, copy, datetime, time, glob, shutil
from collections import OrderedDict
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
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# Import jwst package itself
import jwst

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




# Main
@click.command()
@click.argument('input_cal_files', nargs=-1, type=str)
@click.argument('output_dir', type=click.Path(exists=False))
@click.option('--program', type=str, default='', help='Specify a program.')
@click.option('--obsnum', type=str, default='', help='Specify a obsnum.')
@click.option('--instrument', type=str, default='', help='Specify an instrument.')
@click.option('--filter', 'infilter', type=str, default='', help='Specify a filter.')
@click.option('--group-filter/--no-group-filter', is_flag=True, default=True, help='Process groups of data of the same filter.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_cal_files, 
        output_dir, 
        program, 
        obsnum, 
        instrument, 
        infilter, 
        check_filter, 
        overwrite, 
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
    
    
    # Check output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    
    # expand input file list if there is a wildcard
    input_files = []
    for input_file in input_cal_files:
        if input_file.find('*')>=0:
            for found_file in glob.glob(input_file):
                input_files.append(os.path.abspath(found_file))
        else:
            if not os.path.isfile(input_file):
                raise Exception('Error! File does not exist: {!r}'.format(input_file))
            input_files.append(os.path.abspath(input_file))
    
    
    # sort input files
    input_files.sort(key=LooseVersion)
    
    
    # get info_list filtered by program, instrument, filter
    info_list = []
    for input_filepath in input_files:
        # read fits header
        header = fits.getheader(input_filepath, 0)
        # store into info_dict
        info_dict = OrderedDict()
        info_dict['program'] = header['PROGRAM'].strip()
        info_dict['obs_num'] = header['OBSERVTN'].strip()
        info_dict['target_group'] = header['TARGPROP'].strip()
        info_dict['target_name'] = '"'+header['TARGNAME'].strip()+'"'
        info_dict['target_RA'] = header['TARG_RA']
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
        # check instrument and infilter
        if program != '':
            if program != info_dict['program']:
                continue
        if obsnum != '':
            if obsnum != info_dict['obs_num']:
                continue
        if instrument != '':
            if instrument != info_dict['instrument']:
                continue
        if infilter != '':
            if infilter != info_dict['filter']:
                continue
        info_dict['file_path'].append(input_filepath)
        info_list.append(info_dict)
    
    
    # check len(info_list)
    if len(info_list) == 0:
        if program != '' or obsnum != '' or instrument != '' or infilter != '':
            raise Exception('Error! No input file matching the specified ' + \
                'program {!r} obsnum {!r} instrument {!r} and filter {!r}'.format(
                program, obsnum, instrument, infilter))
        else:
            logger.error('Error! No input file to process!')
            sys.exit(255)
    
    
    # compile info table
    info_table_dict = OrderedDict()
    for i, info_dict in enumerate(info_list):
        if i == 0:
            for key in info_dict:
                info_table_dict[key] = []
        info_table_dict[key].append(info_dict[key])
    info_table = Table(info_table_dict)
    
    
    # sort by '{instrument}_{filter}_{program}_{obs_num}'
    #info_table.sort(['instrument', 'filter', 'program', 'obs_num', 'visit_num'])
    
    
    # output info table
    info_table_file = os.path.join(output_dir, 'info_table')
    if not os.path.isdir(os.path.dirname(info_table_file)):
        os.makedirs(os.path.dirname(info_table_file))
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
    groupped_table = info_table.group_by(['instrument', 'filter'])
    
    
    # prepare association file to process all rate files into one single output file
    for subgroup_key, subgroup_table in zip(groupped_table.groups.keys, groupped_table.groups):
        
        unique_programs = np.unique(subgroup_table['program'])
        groupped_program_str = '+'.join(unique_programs)
        
        unique_obsnums = np.unique(subgroup_table['obs_num'])
        groupped_obsnum_str = '+'.join(unique_obsnums)
        
        subgroup_files = subgroup_table['file_path'].data.tolist()
        
        output_name = 'jw{}_obs{}_{}_{}'.format(
            groupped_program_str, 
            groupped_obsnum_str, 
            subgroup_key['instrument'],
            subgroup_key['filter'])
        
        output_subdir = output_name
        output_file = output_name + '_i2d.fits'
        output_filepath = os.path.join(output_dir, output_subdir, output_file)
        
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
            logger.warning(f"Found another script processing data dir {output_subdir!r}. " + 
                            "We will skip processing it.")
            continue
        
        # create touch file
        with open(output_subdir + '.touch', 'a'):
            os.utime(output_subdir + '.touch', None)
        
        # backup existing file
        if os.path.isfile(output_filepath):
            shutil.move(output_filepath, output_filepath+'.backup')
        
        
        # prepare the json-format "association" file
        asn_dict = OrderedDict()
        asn_dict['asn_type'] = 'None'
        asn_dict['asn_rule'] = 'DMS_Level3_Base'
        asn_dict['version_id'] = None
        asn_dict['code_version'] = jwst.__version__
        asn_dict['degraded_status'] = 'No known degraded exposures in association.'
        asn_dict['program'] = groupped_program_str
        asn_dict['constraints'] = 'No constraints'
        asn_dict['asn_id'] = groupped_obsnum_str
        asn_dict['asn_pool'] = 'none'
        asn_dict['products'] = []
        product_dict = OrderedDict()
        product_dict['name'] = output_name
        product_dict['members'] = []
        asn_dict['products'].append(product_dict)
        for subgroup_file in subgroup_files:
            product_dict['members'].append(
                {'expname': subgroup_file,
                 'exptype': 'science'
                }
            )
        
        asn_file = output_subdir + '_asn.json'
        
        if os.path.isfile(asn_file):
            shutil.move(asn_file, asn_file+'.backup')
        
        with open(asn_file, 'w') as fp:
            json.dump(asn_dict, fp, indent=4)
        
        
        # Print progress
        logger.info("Processing {} -> {}".format(subgroup_files, output_filepath))
        
        
        # prepare to run
        pipeline_object = calwebb_image3.Image3Pipeline()
        

        # Set some parameters that pertain to the entire pipeline
        pipeline_object.output_dir = output_subdir
        pipeline_object.output_file = output_name # os.path.splitext(output_file)[0]
        pipeline_object.output_ext = ".fits" # default
        pipeline_object.save_results = True

        # Set some parameters that pertain to some of the individual steps
        # Set OutlierDetection
        #pipeline_object.outlier_detection.skip = True
        pipeline_object.outlier_detection.output_dir = output_dir
        # Turn on TweakRegStep
        #pipeline_object.tweakreg.skip = True
        #pipeline_object.tweakreg.save_catalogs = True
        #pipeline_object.tweakreg.save_results = True
        #pipeline_object.tweakreg.output_dir = output_dir
        # Turn on SkyMatchStep
        #pipeline_object.skymatch.skip = True
        pipeline_object.skymatch.subtract = True
        pipeline_object.skymatch.skymethod = "global+match"
        #pipeline_object.skymatch.lsigma = 2.0
        #pipeline_object.skymatch.usigma = 2.0
        #pipeline_object.skymatch.nclip = 10
        #pipeline_object.skymatch.upper = 1.0
        pipeline_object.skymatch.save_results = True
        # Set the ratio of input to output pixels to create an output mosaic 
        # on a 0.015"/pixel scale
        pipeline_object.resample.pixel_scale_ratio = 0.48
        
        
        # run
        pipeline_object.run(asn_file)
        
        
        # remove touch file
        os.remove(output_subdir + '.touch')
        
        
        # check
        assert os.path.isfile(output_filepath)
        
        
        # log
        logger.info("Processed {} -> {}".format(subgroup_files, output_filepath))
    
    
    # log
    logger.info("All done!")




# Entry
if __name__ == '__main__':
    main()



