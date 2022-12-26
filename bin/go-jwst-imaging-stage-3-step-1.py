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
import os, sys, re, json, copy, datetime, time, glob, shutil
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
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# Import jwst package itself
import jwst

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
@click.option('--combine-program/--no-combine-program', is_flag=True, default=False, help='Combine all programs into one.')
@click.option('--combine-obsnum/--no-combine-obsnum', is_flag=True, default=False, help='Combine all obsnum into one.')
@click.option('--combine-filter/--no-combine-filter', is_flag=True, default=False, help='Combine all filters into one.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False, help='Overwrite?')
@click.option('--run-individual-steps/--no-run-individual-steps', is_flag=True, default=False, help='Run individual step of JWST stage3 pipeline? This is turned on if abs_refcat is provided!')
def main(
        input_cal_files, 
        output_dir, 
        select_program, 
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
        combine_program, 
        combine_obsnum, 
        combine_filter, 
        overwrite, 
        run_individual_steps,
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
            if info_dict['program'] not in select_program:
                continue
        if len(select_obsnum) > 0:
            if info_dict['obs_num'] not in select_obsnum:
                continue
        if len(select_instrument) > 0:
            if info_dict['instrument'] not in select_instrument:
                continue
        if len(select_filter) > 0:
            if info_dict['filter'] not in select_filter:
                continue
        info_list.append(info_dict)
    
    
    # check len(info_list)
    if len(info_list) == 0:
        if len(select_program) > 0 or \
           len(select_obsnum) > 0 or \
           len(select_instrument) > 0 or \
           len(select_filter) > 0:
            raise Exception('Error! No input file matching the specified ' + \
                'program {!r} obsnum {!r} instrument {!r} and filter {!r}'.format(
                select_program, select_obsnum, select_instrument, select_filter))
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
    group_by_columns = ['instrument', 'filter', 'program', 'obs_num']
    if combine_filter:
        if len(np.unique(info_table['filter'])) > 0:
            logger.warning('We are combining different filters!')
        group_by_columns.remove('filter')
    if combine_program:
        if len(np.unique(info_table['program'])) > 0:
            logger.warning('We are combining different programs!')
        group_by_columns.remove('program')
    if combine_obsnum:
        if len(np.unique(info_table['obs_num'])) > 0:
            logger.warning('We are combining different obs_num!')
        group_by_columns.remove('obs_num')
    logger.info('group_by_columns: {}'.format(group_by_columns))
    groupped_table = info_table.group_by(group_by_columns)
    
    
    # prepare association file to process all rate files into one single output file
    output_lookup_dict = OrderedDict()
    for subgroup_key, subgroup_table in zip(groupped_table.groups.keys, groupped_table.groups):
        
        unique_programs = np.unique(subgroup_table['program'])
        groupped_program_str = '+'.join(unique_programs)
        
        unique_obsnums = np.unique(subgroup_table['obs_num'])
        groupped_obsnum_str = '+'.join(unique_obsnums)
        
        unique_instruments = np.unique(subgroup_table['instrument'])
        groupped_instrument_str = '+'.join(unique_instruments)
        
        unique_filters = np.unique(subgroup_table['filter'])
        groupped_filter_str = '+'.join(unique_filters)
        
        #unique_obsnums = np.unique(subgroup_table['obs_num'])
        #groupped_obsnum_str = '+'.join(unique_obsnums)
        
        subgroup_files = subgroup_table['file_path'].data.tolist()
        
        output_name = 'jw{}_obs{}_{}_{}'.format(
            groupped_program_str, 
            groupped_obsnum_str, 
            groupped_instrument_str,
            groupped_filter_str,
        )
        
        output_subdir = os.path.join(output_dir, output_name)
        output_file = output_name + '_i2d.fits'
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
            logger.warning(f"WARNING WARNING WARNING! Found another script processing data dir {output_subdir!r}. " + 
                            "We will skip processing it for now!")
            continue
        
        # create touch file
        with open(output_subdir + '.touch', 'a'):
            os.utime(output_subdir + '.touch', None)
        
        # backup existing file
        if os.path.isfile(output_filepath):
            shutil.move(output_filepath, output_filepath+'.backup')
        
        # clean up existing cal.fits under output_subdir
        #for calfilename in os.listdir(output_subdir):
        #    if calfilename.endswith('_cal.fits'):
        #        if os.path.islink()
        
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
        catdict = OrderedDict() # see "jwst/tweakreg/tweakreg_step.py"
        for subgroup_file in subgroup_files:
            # get 'expname'
            # relative to asn file dir path, 
            # see "jwst/datamodels/container.py"
            # see also "jwst/tweakreg/tweakreg_step.py", 
            # ```
            # filename = member['expname']
            # member['expname'] = path.join(asn_dir, filename)
            # ```
            #expname = os.path.relpath(subgroup_file, os.path.dirname(output_subdir))
            # 
            # 20221203: now using symlink and run the processing from the output_subdir!
            expname = os.path.basename(subgroup_file)
            explink = os.path.join(output_subdir, expname)
            product_dict['members'].append(
                {'expname': expname, 
                 'exptype': 'science'
                }
            )
            if not os.path.exists(explink):
                if os.path.islink(explink):
                    os.remove(explink)
                os.symlink(os.path.relpath(subgroup_file, output_subdir), 
                           explink)
            # 
            # if there is a '{dataset_cal}_cat_for_tweakreg.csv' file, 
            # use it as catfile for tweakreg
            calcat_file = os.path.splitext(subgroup_file)[0]+'_cat_for_tweakreg.csv'
            if os.path.isfile(calcat_file):
                calcat_name = os.path.basename(calcat_file)
                calcat_link = os.path.join(output_subdir, calcat_name)
                if not os.path.exists(calcat_link):
                    if os.path.islink(calcat_link):
                        os.remove(calcat_link)
                    os.symlink(os.path.relpath(calcat_file, output_subdir), 
                               calcat_link)
                catkey = expname
                catdict[catkey] = calcat_name
        
        asn_filename = 'asn.json'
        asn_file = os.path.join(output_subdir, asn_filename)
        
        if os.path.isfile(asn_file):
            shutil.move(asn_file, asn_file+'.backup')
        
        with open(asn_file, 'w') as fp:
            json.dump(asn_dict, fp, indent=4)
        
        catfile = None
        if len(catdict) > 0:
            catfilename = 'catfile.txt'
            catfilepath = os.path.join(output_subdir, catfilename)
            if os.path.isfile(catfilepath):
                shutil.move(catfilepath, catfilepath+'.backup')
            with open(catfilepath, 'w') as fp:
                for catkey in catdict:
                    fp.write('{} {}\n'.format(catkey, catdict[catkey]))
            catfile = catfilename # we will run the process inside the output_subdir, so set only filename here
        
        
        if abs_refcat is not None and abs_refcat != '':
            abs_refcat = os.path.abspath(abs_refcat)
            run_individual_steps = True
        
        
        # Print progress
        logger.info("Processing {} -> {}".format(subgroup_files, output_filepath))
        
        
        # chdir
        current_dir = os.getcwd()
        logger.info("chdir {}".format(output_subdir))
        os.chdir(output_subdir)
        
        
        # prepare to run
        pipeline_object = calwebb_image3.Image3Pipeline()
        

        # Set some parameters that pertain to the entire pipeline
        #pipeline_object.input_dir = os.getcwd()
        #[20221203]#pipeline_object.output_dir = output_subdir
        pipeline_object.output_file = output_name # os.path.splitext(output_file)[0]
        pipeline_object.output_ext = ".fits" # default
        pipeline_object.save_results = True

        # Set some parameters that pertain to some of the individual steps
        # Turn on TweakRegStep
        #pipeline_object.tweakreg.skip = False
        #[20221203]#pipeline_object.tweakreg.output_dir = output_subdir
        pipeline_object.tweakreg.minobj = 7 # default is 15
        pipeline_object.tweakreg.searchrad = 1.0 # default is 2.0
        pipeline_object.tweakreg.separation = 1.0 # default is 0.1 arcsec
        pipeline_object.tweakreg.tolerance = 1.0 # default is 0.7 arcsec
        pipeline_object.tweakreg.save_catalogs = True # "./*_cal_cat.ecsv" # always current dir, see "tweakreg_step.py"
        pipeline_object.tweakreg.save_results = True # "{output_subdir}/*_cal_tweakreg.fits"
        pipeline_object.tweakreg.search_output_file = False # do not use output_file define in parent step
        pipeline_object.tweakreg.output_use_model = True # use DataModel.meta.filename as output_file
        # set catfile if found
        if catfile is not None:
            pipeline_object.tweakreg.use_custom_catalogs = True
            pipeline_object.tweakreg.catfile = catfile
        # set abs_refcat if user has input that
        if abs_refcat is not None and abs_refcat != '':
            pipeline_object.tweakreg.abs_refcat = abs_refcat
            #pipeline_object.abs_fitgeometry = 'rshift' # 'shift', 'rshift', 'rscale', 'general' (shift, rotation, and scale)
            pipeline_object.tweakreg.abs_minobj = 1 # default is 15
            pipeline_object.tweakreg.abs_searchrad = 1.0 # default is 6.0
            pipeline_object.tweakreg.abs_separation = 1.0 # default is 0.1 arcsec
            pipeline_object.tweakreg.abs_tolerance = 1.0 # default is 0.7 arcsec
            pipeline_object.tweakreg.abs_use2dhist = False # default is True, but ...
            #pipeline_object.tweakreg.abs_fitgeometry = 'shift' # 'rshift'
            pipeline_object.tweakreg.save_abs_catalog = True
            #pipeline_object.tweakreg.save_abs_catalog = True # only if abs_refcat `gaia_cat_name in SINGLE_GROUP_REFCAT`
        # do manual 2dhist for a better tweakreg
        if catfile is not None and abs_refcat is not None and abs_refcat != '':
            from util_match_cat_file_to_abs_refcat_with_2dhist import match_cat_file_to_abs_refcat_with_2dhist
            pipeline_object.tweakreg.catfile = match_cat_file_to_abs_refcat_with_2dhist(catfile, abs_refcat)
        
        # Turn on SkyMatchStep
        #pipeline_object.skymatch.skip = False
        pipeline_object.skymatch.subtract = True
        pipeline_object.skymatch.skymethod = "global+match"
        #pipeline_object.skymatch.lsigma = 2.0
        #pipeline_object.skymatch.usigma = 2.0
        #pipeline_object.skymatch.nclip = 10
        #pipeline_object.skymatch.upper = 1.0
        pipeline_object.skymatch.save_results = True
        
        # Set OutlierDetection
        #pipeline_object.outlier_detection.skip = False
        #[20221203]#pipeline_object.outlier_detection.output_dir = output_subdir
        pipeline_object.outlier_detection.pixfrac = pixfrac
        pipeline_object.outlier_detection.in_memory = False
        #pipeline_object.outlier_detection.suffix = 'crf'
        #pipeline_object.outlier_detection.save_results = pipeline_object.save_results

        
        # Set the ratio of input to output pixels to create an output mosaic 
        # on a 0.015"/pixel scale
        # see -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
        #pipeline_object.resample_data = True # True is the default
        #pipeline_object.resample.suffix = 'i2d'
        #pipeline_object.resample.save_results = pipeline_object.save_results
        pipeline_object.resample.kernel = kernel
        pipeline_object.resample.pixfrac = pixfrac
        #pipeline_object.resample.pixel_scale_ratio = 0.48
        pipeline_object.resample.pixel_scale_ratio = pixel_scale_ratio
        pipeline_object.resample.pixel_scale = pixel_scale
        
        # Set source_catalog
        #pipeline_object.source_catalog.save_results = pipeline_object.save_results
        
        
        # run
        if not run_individual_steps:
            
            pipeline_object.run(asn_filename)
        
        # we can also run individual steps following "jwst/pipeline/calwebb_image3.py" `process()`
        else: # TODO
            
            pipeline_object.log.info(
                f'Step {pipeline_object.name} running with asn file {asn_filename}.'
            )
            pipeline_object.log.info(
                f'Step {pipeline_object.name} parameters are: {pipeline_object.get_pars()}'
            )
            asn_exptypes = ['science']
            with datamodels.open(asn_filename, asn_exptypes=asn_exptypes) as input_models:
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
                pipeline_object.tweakreg.save_model(image_models, idx=None, suffix='tweakreg',
                                format=pipeline_object.tweakreg.name_format, force=True)
                # 
                # 2. skymatch bkgmatch
                pipeline_object.skymatch.output_file = output_name
                pipeline_object.skymatch.save_results = True
                pipeline_object.skymatch.search_output_file = False
                image_models = pipeline_object.skymatch(image_models)
                # 
                # 3. outlier detection
                # TODO: here we can split image_models into each obsnum group, 
                # then parallelize the outlier_detection process,
                # this also saves memory usage
                pipeline_object.outlier_detection.output_file = output_name
                pipeline_object.outlier_detection.save_results = True
                pipeline_object.outlier_detection.suffix = 'crf'
                pipeline_object.outlier_detection.search_output_file = False
                if len(unique_obsnums) == 1:
                    image_models = pipeline_object.outlier_detection(image_models)
                else:
                    for subsubgrup_obsnum in unique_obsnums:
                        subsubgroup_indices = np.argwhere(subgroup_table['obs_num'] == subsubgrup_obsnum).ravel()
                        subsubgroup_image_models = np.take(image_models, subsubgroup_indices)
                        pipeline_object.log.info(
                            'Running outlier_detection for models: {} (obsnum: {}, indices: {})'.format(
                                subsubgroup_image_models, subsubgrup_obsnum, subsubgroup_indices)
                        )
                        subsubgroup_image_models = pipeline_object.outlier_detection(subsubgroup_image_models)
                        for kk in range(len(subsubgroup_indices)):
                            image_models[subsubgroup_indices[kk]] = subsubgroup_image_models[kk]
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
        
        
        # move some output files to subdirectories
        # for example the "outlier_i2d.fits", so that we do not find them when doing an `ls *_i2d.fits`
        for i in range(len(subgroup_files)):
            subgroup_file = subgroup_files[i]
            dataset_name = subgroup_file.replace('_cal.fits', '')
            if os.path.isfile(f'{output_name}_{i}_outlier_i2d.fits'):
                if not os.path.isdir(f'{output_name}_{i}_outlier'):
                    os.makedirs(f'{output_name}_{i}_outlier')
                shutil.move(f'{output_name}_{i}_outlier_i2d.fits', 
                            f'{output_name}_{i}_outlier/{output_name}_{i}_outlier_i2d.fits')
        
        
        # chdir back
        logger.info("chdir {}".format(current_dir))
        os.chdir(current_dir)
        
        
        # save pars as asdf file
        asdf_filepath = os.path.splitext(output_filepath)[0] + '_calwebb_image3.asdf'
        if os.path.isfile(asdf_filepath):
            shutil.move(asdf_filepath, asdf_filepath+'.backup')
        asdf_object = AsdfFile(pipeline_object.get_pars())
        asdf_object.write_to(asdf_filepath)
        logger.info('Parameters are saved into {}'.format(asdf_filepath))
        
        
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



