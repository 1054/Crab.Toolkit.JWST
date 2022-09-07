#!/usr/bin/env python
#
"""
Process JWST rate data under "calibrated/" and named as "jw*_cal.fits"
into calibrated FITS data under "calibrated/" and named as "jw*_i2d.fits". 

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
import asdf
from collections import OrderedDict

# Numpy library:
import numpy as np

# For downloading data
import requests

# Astropy tools:
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch

# matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300

# Set default pipeline context
if 'CRDS_PATH' not in os.environ:
    os.environ['CRDS_PATH'] = os.path.expanduser('~/jwst_crds_cache')
    print('Setting {!r} to {!r}'.format('CRDS_PATH', os.environ['CRDS_PATH']))
if 'CRDS_SERVER_URL' not in os.environ:
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
    print('Setting {!r} to {!r}'.format('CRDS_SERVER_URL', os.environ['CRDS_SERVER_URL']))
if 'CRDS_CONTEXT' not in os.environ:
    os.environ['CRDS_CONTEXT'] = 'jwst_0968.pmap' # TODO
    print('Setting {!r} to {!r}'.format('CRDS_CONTEXT', os.environ['CRDS_CONTEXT']))

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
if __name__ == '__main__':

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
    
    # Read command line input arguments
    overwrite = False
    prefixes = []
    iarg = 1
    while iarg < len(sys.argv):
        if sys.argv[iarg] == '--overwrite' or sys.argv[iarg] == '-overwrite':
            overwrite = True
        else:
            prefixes.append(sys.argv[iarg])
        iarg += 1
    
    # 
    if len(prefixes) <= 1:
        logger.error("Error! Please input at least two prefixes so that we can search for \"{prefix}*/calibrated2_cals/{prefix}*_cal.fits\" files to make mosaic.")
        sys.exit(-1)
    
    suffix_str = '_' + '+'.join(prefixes)
    
    # Find all "jw*/calibrated2_cals/jw*_cal.fits"
    input_files = []
    for prefix in prefixes:
        # find files to make mosaic
        found_files = glob.glob(f"{prefix}*/calibrated2_cals/{prefix}*_cal.fits")
        found_files = [t for t in found_files if t.find('remstriping') < 0]
        if (len(found_files) == 0):
            logger.error(f"Error! No input file \"{prefix}*/calibrated2_cals/{prefix}*_cal.fits\" is found!")
            sys.exit(-1)
        input_files.extend(found_files)
    
    
    
    input_files = sorted(input_files)
    
    
    # loop individual files, find unique 'obs_id' and 'target_name'
    info_dict = OrderedDict()
    info_dict['program'] = []
    info_dict['obs_id'] = []
    info_dict['target_group'] = []
    info_dict['target_name'] = []
    info_dict['target_RA'] = []
    info_dict['target_Dec'] = []
    info_dict['instrument'] = []
    info_dict['filter'] = []
    info_dict['pupil'] = []
    info_dict['file_path'] = []
    for input_filepath in input_files:
        # read fits header
        header = fits.getheader(input_filepath, 0)
        info_dict['program'].append(header['PROGRAM'].strip())
        info_dict['obs_id'].append(header['OBSERVTN'].strip())
        info_dict['target_group'].append(header['TARGPROP'].strip())
        info_dict['target_name'].append('"{}"'.format(header['TARGNAME'].strip()))
        info_dict['target_RA'].append(header['TARG_RA'])
        info_dict['target_Dec'].append(header['TARG_DEC'])
        info_dict['instrument'].append(header['INSTRUME'].strip())
        info_dict['filter'].append(header['FILTER'].strip())
        if 'PUPIL' in header:
            info_dict['pupil'].append(header['PUPIL'].strip())
        else:
            info_dict['pupil'].append('""')
        info_dict['file_path'].append(input_filepath)
    
    info_table = Table(info_dict)
    
    # sort by '{program}_{obs_id}_{instrument}_{filter}'
    info_table.sort(['instrument', 'filter', 'program', 'obs_id'])
    
    if not os.path.isdir('calibrated3_mosaics'):
        os.makedirs('calibrated3_mosaics')
    if os.path.isfile(f'calibrated3_mosaics/info_table{suffix_str}.txt'):
        shutil.move(f'calibrated3_mosaics/info_table{suffix_str}.txt', f'calibrated3_mosaics/info_table{suffix_str}.txt.backup')
    if os.path.isfile(f'calibrated3_mosaics/info_table{suffix_str}.csv'):
        shutil.move(f'calibrated3_mosaics/info_table{suffix_str}.csv', f'calibrated3_mosaics/info_table{suffix_str}.csv.backup')
    info_table.write(f'calibrated3_mosaics/info_table{suffix_str}.txt', format='ascii.fixed_width', delimiter=' ', bookend=True)
    with open(f'calibrated3_mosaics/info_table{suffix_str}.txt', 'r+') as fp:
        fp.seek(0)
        fp.write('#')
    info_table.write(f'calibrated3_mosaics/info_table{suffix_str}.csv', format='csv')
    
    
    # group by '{program}_{obs_id}_{instrument}_{filter}'
    unique_groups = info_table.group_by(['instrument', 'filter'])
    
    
    # make a string of all unique programs for output filename
    unique_programs = np.unique(info_table['program']).tolist()
    unique_programs = sorted(unique_programs)
    all_program_str = '+'.join(unique_programs)
    
    unique_program_obs_ids = []
    for i in range(len(info_table)):
        program_obs_id_str = '{}{}'.format(info_table['program'][i], info_table['obs_id'][i])
        if program_obs_id_str not in unique_program_obs_ids:
            unique_program_obs_ids.append(program_obs_id_str)
    unique_program_obs_ids = sorted(unique_program_obs_ids)
    all_program_obs_id_str = '+'.join(unique_program_obs_ids)
    
    
    # prepare association file to process all rate files into one single output file
    for subgroup_key, subgroup_table in zip(unique_groups.groups.keys, unique_groups.groups):
        
        input_files = subgroup_table['file_path'].data.tolist()
        output_name = 'combined_{}_{}_{}'.format(
            all_program_str, 
            subgroup_key['instrument'],
            subgroup_key['filter'])
        output_dir = os.path.join('calibrated3_mosaics', output_name) # directly output to "calibrated3_mosaics" under current directory.
        output_file = output_name+'_i2d.fits'
        output_filepath = os.path.join(output_dir, output_file)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        else:
            # check existing output file, skip the processing or backup and overwrite it. 
            if os.path.isfile(output_filepath):
                if not overwrite:
                    logger.info("Found processed {} -> {} and overwrite is set to False. We will skip the processing.".format(input_files, output_filepath))
                    continue
        
        # check if another script is processing this data
        if os.path.exists(output_dir+'.touch'):
            logger.info("Found another script processing {} -> {}. We will skip processing it.".format(input_files, output_filepath))
            continue
        
        # create touch file
        with open(output_dir+'.touch', 'a'):
            os.utime(output_dir+'.touch', None)
        
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
        asn_dict['program'] = all_program_str # TODO
        asn_dict['constraints'] = 'No constraints' # TODO
        asn_dict['asn_id'] = all_program_obs_id_str # TODO
        asn_dict['asn_pool'] = 'none'
        asn_dict['products'] = []
        product_dict = OrderedDict()
        product_dict['name'] = output_name
        product_dict['members'] = []
        asn_dict['products'].append(product_dict)
        for input_filepath in input_files:
            input_filename = os.path.basename(input_filepath)
            input_suffix = '_cal'
            product_dict['members'].append(
                {'expname': os.path.join('..', input_filepath),
                 'exptype': 'science'}
            )
        
        asn_file = os.path.join('calibrated3_mosaics', output_name+'_asn.json')
        
        if os.path.isfile(asn_file):
            shutil.move(asn_file, asn_file+'.backup')
        
        with open(asn_file, 'w') as fp:
            json.dump(asn_dict, fp, indent=4)
        
        
        # prepare a single output file 
        
        logger.info("Processing {} -> {}".format(input_files, output_filepath))
        
        
        # prepare to run
        pipeline_object = calwebb_image3.Image3Pipeline()

        # Set some parameters that pertain to the entire pipeline
        pipeline_object.output_dir = output_dir
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
        os.remove(output_dir+'.touch')
        
        # check
        assert os.path.isfile(output_filepath)
        
        # log
        logger.info("Processed {} -> {}".format(input_files, output_filepath))
    
    
    
    
    



