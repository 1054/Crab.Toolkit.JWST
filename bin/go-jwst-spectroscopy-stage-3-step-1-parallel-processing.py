#!/usr/bin/env python
#
"""
Process a number of JWST cal data, output multi-exposure combined cal data. 

This script runs the JWST calwebb_spec3 pipeline.

Inputs
    
    jw*_cal.fits

Outputs

    jw*-o*_s*_*_cal.fits
    
Last update:
    
    2022-12-24 DZLIU. 

Notes:
    

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

# The entire calwebb_spec3 pipeline
from jwst import datamodels
from jwst.pipeline import calwebb_spec3
from jwst.pipeline import calwebb_image3
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# Packages used in calwebb_spec3.py
from collections import defaultdict
import os.path as op
from jwst.associations.lib.rules_level3_base import format_product
from jwst.exp_to_source import multislit_to_container

# Import jwst package itself
import jwst

# Import AsdfFile
from asdf import AsdfFile

# Setup logging
import logging

# Setup multiprocessing
import multiprocessing as mp

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

def process_asn_async(asn_filename, output_name, index_list, overwrite):
    # only for WFSS
    # following the code in calwebb_spec3.py
    pipeline_object = calwebb_spec3.Spec3Pipeline()
    pipeline_object.output_file = output_name
    pipeline_object.save_results = True
    pipeline_object.master_background.skip = True
    pipeline_object.mrs_imatch.skip = True
    pipeline_object.outlier_detection.skip = True
    pipeline_object.resample_spec.skip = True
    pipeline_object.cube_build.skip = True
    pipeline_object.extract_1d.suffix = 'x1d'
    pipeline_object.extract_1d.save_results = pipeline_object.save_results
    pipeline_object.combine_1d.suffix = 'c1d'
    pipeline_object.combine_1d.save_results = pipeline_object.save_results
    pipeline_object.spectral_leak.skip = True
    pipeline_object.pixel_replace.skip = True
    asn_exptypes = ['science', 'background']
    input_models = datamodels.open(asn_filename, asn_exptypes=asn_exptypes)
    for model in input_models:
        model.meta.asn.table_name = op.basename(input_models.asn_table_name)
    # 
    exptype = input_models[0].meta.exposure.type # 'NRC_WFSS'
    model_type = input_models[0].meta.model_type # 'MultiSlitModel'
    output_file = input_models.meta.asn_table.products[0].name # 'jw05893_obs017_NIRCAM_F444W'
    # 
    members_by_type = defaultdict(list)
    product = input_models.meta.asn_table.products[0].instance
    for member in product['members']:
        members_by_type[member['exptype'].lower()].append(member['expname'])
    # 
    source_models = input_models
    sources = [source_models]
    if model_type in calwebb_spec3.MULTISOURCE_MODELS:
        sources = [(name, model) for name, model in multislit_to_container(source_models).items()]
    # 
    for i in index_list:
        if i>=0 and i < len(sources):
            source = sources[i]
            if isinstance(source, tuple):
                source_id, result = source
                print(f'DEBUG: processing source_id {source_id}')

                # NIRSpec fixed-slit data
                if exptype == "NRS_FIXEDSLIT":
                    # Output file name is constructed using the source_id and the slit name
                    slit_name = pipeline_object._create_nrsfs_slit_name(result)
                    srcid = f's{source_id:>09s}'
                    pipeline_object.output_file = format_product(output_file, source_id=srcid, slit_name=slit_name)

                # NIRSpec MOS/MSA data
                elif exptype == "NRS_MSASPEC":
                    # Construct the specially formatted source_id to use in the output file
                    # name that separates source, background, and virtual slits
                    srcid = pipeline_object._create_nrsmos_source_id(result)
                    pipeline_object.output_file = format_product(output_file, source_id=srcid)
                    logger.debug(f"output_file = {pipeline_object.output_file}")

                else:
                    # All other types just use the source_id directly in the file name
                    srcid = f's{source_id:>09s}'
                    pipeline_object.output_file = format_product(output_file, source_id=srcid)

            else:
                result = source

            # 
            if not overwrite:
                print('Checking output {}'.format(pipeline_object.output_file))
                if os.path.isfile(pipeline_object.output_file):
                    print('Skipping output {}'.format(pipeline_object.output_file))
                    continue

            # 
            resample_complete = None

            # Call outlier detection and pixel replacement
            if exptype not in calwebb_spec3.SLITLESS_TYPES:
                # Update the asn table name to the level 3 instance so that
                # the downstream products have the correct table name since
                # the _cal files are not saved they will not be updated
                for cal_array in result:
                    cal_array.meta.asn.table_name = op.basename(input_models.asn_table_name)
                if exptype in calwebb_spec3.IFU_EXPTYPES:
                    pipeline_object.outlier_detection.mode = 'ifu'
                else:
                    pipeline_object.outlier_detection.mode = 'spec'
                result = pipeline_object.outlier_detection(result)

                # interpolate pixels that have a NaN value or are flagged
                # as DO_NOT_USE or NON_SCIENCE.
                result = pipeline_object.pixel_replace(result)
                # Resample time. Dependent on whether the data is IFU or not.
                resample_complete = None
                if exptype in calwebb_spec3.IFU_EXPTYPES:
                    result = pipeline_object.cube_build(result)
                    try:
                        resample_complete = result[0].meta.cal_step.cube_build
                    except AttributeError:
                        pass
                else:
                    result = pipeline_object.resample_spec(result)
                    try:
                        resample_complete = result.meta.cal_step.resample
                    except AttributeError:
                        pass

            # Do 1-D spectral extraction
            if exptype in calwebb_spec3.SLITLESS_TYPES:

                # interpolate pixels that have a NaN value or are flagged
                # as DO_NOT_USE or NON_SCIENCE
                result = pipeline_object.pixel_replace(result)

                # For slitless data, extract 1D spectra and then combine them
                if exptype in ['NIS_SOSS']:
                    # For NIRISS SOSS, don't save the extract_1d results,
                    # instead run photom on the extract_1d results and save
                    # those instead.

                    pipeline_object.extract_1d.save_results = False
                    result = pipeline_object.extract_1d(result)

                    # SOSS F277W may return None - don't bother with that.
                    if result is not None:
                        pipeline_object.photom.save_results = pipeline_object.save_results
                        pipeline_object.photom.suffix = 'x1d'
                        result = pipeline_object.photom(result)
                else:
                    result = pipeline_object.extract_1d(result)

                result = pipeline_object.combine_1d(result)

            elif resample_complete is not None and resample_complete.upper() == 'COMPLETE':

                # If 2D data were resampled and combined, just do a 1D extraction

                if exptype in calwebb_spec3.IFU_EXPTYPES:
                    pipeline_object.extract_1d.search_output_file = False
                    if exptype in ['MIR_MRS']:
                        if not pipeline_object.spectral_leak.skip:
                            pipeline_object.extract_1d.save_results = False
                            pipeline_object.spectral_leak.suffix = 'x1d'
                            pipeline_object.spectral_leak.search_output_file = False
                            pipeline_object.spectral_leak.save_results = pipeline_object.save_results

                result = pipeline_object.extract_1d(result)

                if exptype in ['MIR_MRS']:
                    result = pipeline_object.spectral_leak(result)
            
            # else:
            #     logger.warning(
            #         'Resampling was not completed. Skipping extract_1d.'
            #     )

        input_models.close()

        #self.log.info('Ending calwebb_spec3')






# Defaults
# see -- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/arguments.html
# DEFAULT_KERNEL = 'square'
# DEFAULT_PIXFRAC = 1.0 # 0.5
# DEFAULT_PIXEL_SCALE_RATIO = 0.48
# DEFAULT_PIXEL_SCALE = None



# Main
@click.command()
@click.argument('input_cal_files', nargs=-1, type=str)
@click.argument('output_dir', type=click.Path(exists=False))
@click.option('--nparallel', type=int, default=4, help='Run the processing (mostly spectral extraction) in N parallel processes. Default is 4.')
@click.option('--nsource-per-process', type=int, default=30, help='In each parallel process, how many sources to be processed. Default is 30.')
@click.option('--program', 'select_program', type=str, multiple=True, default=None, help='Specify a program.')
@click.option('--obsnum', 'select_obsnum', type=str, multiple=True, default=None, help='Specify a obsnum.')
@click.option('--instrument', 'select_instrument', type=str, multiple=True, default=None, help='Specify an instrument.')
@click.option('--filter', 'select_filter', type=str, multiple=True, default=None, help='Specify a filter.')
# @click.option('--kernel', type=click.Choice(['point', 'square', 'gaussian', 'tophat', 'turbo', 'lanczos2', 'lanczos3']), 
#                           default=DEFAULT_KERNEL, 
#                           help='Drizzle kernel.'+
#                                'The form of the kernel function used to distribute flux onto the output image. '+
#                                'Available kernels are square, gaussian, point, tophat, turbo, lanczos2, and lanczos3.')
# @click.option('--pixfrac', type=float, 
#                            default=DEFAULT_PIXFRAC, 
#                            help='Drizzle pixfrac.')
# @click.option('--pixel-scale-ratio', type=float, 
#                                      default=DEFAULT_PIXEL_SCALE_RATIO, 
#                                      help='Drizzle pixel scale ratio. Ratio of input to output pixel scale. '+
#                                           'A value of 0.5 means the output image would have 4 pixels sampling '+
#                                           'each input pixel. Ignored when --pixel-scale are provided.')
# @click.option('--pixel-scale', type=float, 
#                                default=DEFAULT_PIXEL_SCALE, 
#                                help='Drizzle pixel scale. '+
#                                     'Absolute pixel scale in arcsec. '+
#                                     'When provided, overrides pixel_scale_ratio.')
# @click.option('--abs-refcat', type=click.Path(exists=True), 
#                               default=None, 
#                               help='Absolute reference catalog, must contain `RA` and `DEC` columns, optionally `weight`. See `tweakwcs/imalign.py` `align_wcs`.')
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
#@click.option('--run-individual-steps/--no-run-individual-steps', is_flag=True, default=False, help='Run individual step of JWST stage3 pipeline? This is turned on if abs_refcat is provided!')
def main(
        input_cal_files, 
        output_dir, 
        nparallel, 
        nsource_per_process, 
        select_program, 
        select_obsnum, 
        select_instrument, 
        select_filter, 
        #kernel, 
        #pixfrac, 
        #pixel_scale_ratio, 
        #pixel_scale, 
        #abs_refcat, 
        save_info_table_dir, 
        save_info_table_name, 
        combine_program, 
        combine_obsnum, 
        combine_filter, 
        overwrite, 
        #run_individual_steps,
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
        
        unique_obsnums = np.unique(subgroup_table['obs_num'])
        groupped_obsnum_str = '+'.join(unique_obsnums)
        
        subgroup_files = subgroup_table['file_path'].data.tolist()
        
        output_name = 'jw{}_obs{}_{}_{}'.format(
            groupped_program_str, 
            groupped_obsnum_str, 
            groupped_instrument_str,
            groupped_filter_str,
        )
        
        output_subdir = os.path.join(output_dir, output_name)
        
        # check meta.exposure.type and set output file name
        with datamodels.open(subgroup_files[0]) as model:
            exp_type = model.meta.exposure.type
        if exp_type in ['NRS_MSATA', 'NRS_TACONFIRM']:
            target_acquisition_mode = True
            output_file = output_name + '_i2d.fits'
        else:
            target_acquisition_mode = False
            output_file = output_name + '_all_sources.csv' # here we output a source table
        
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
        # catdict = OrderedDict() # see "jwst/tweakreg/tweakreg_step.py"
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
            # # if there is a '{dataset_cal}_cat_for_tweakreg.csv' file, 
            # # use it as catfile for tweakreg
            # calcat_file = os.path.splitext(subgroup_file)[0]+'_cat_for_tweakreg.csv'
            # if os.path.isfile(calcat_file):
            #     calcat_name = os.path.basename(calcat_file)
            #     calcat_link = os.path.join(output_subdir, calcat_name)
            #     if not os.path.exists(calcat_link):
            #         if os.path.islink(calcat_link):
            #             os.remove(calcat_link)
            #         os.symlink(os.path.relpath(calcat_file, output_subdir), 
            #                    calcat_link)
            #     catkey = expname
            #     catdict[catkey] = calcat_name
        
        asn_filename = 'asn.json'
        asn_file = os.path.join(output_subdir, asn_filename)
        
        if os.path.isfile(asn_file):
            shutil.move(asn_file, asn_file+'.backup')
        
        with open(asn_file, 'w') as fp:
            json.dump(asn_dict, fp, indent=4)
        
        # # catfile = None
        # # if len(catdict) > 0:
        # #     catfilename = 'catfile.txt'
        # #     catfilepath = os.path.join(output_subdir, catfilename)
        # #     if os.path.isfile(catfilepath):
        # #         shutil.move(catfilepath, catfilepath+'.backup')
        # #     with open(catfilepath, 'w') as fp:
        # #         for catkey in catdict:
        # #             fp.write('{} {}\n'.format(catkey, catdict[catkey]))
        # #     catfile = catfilename # we will run the process inside the output_subdir, so set only filename here
        
        
        # if abs_refcat is not None and abs_refcat != '':
        #     abs_refcat = os.path.abspath(abs_refcat)
        #     run_individual_steps = True
        
        
        # Print progress
        logger.info("Processing {} -> {}".format(subgroup_files, output_filepath))
        
        
        # chdir
        current_dir = os.getcwd()
        logger.info("chdir {}".format(output_subdir))
        os.chdir(output_subdir)
        
        
        # check meta.exposure.type and prepare to run
        if target_acquisition_mode:
            pipeline_object = calwebb_image3.Image3Pipeline()
            pipeline_object.source_catalog.skip = True # 'apcorr' in pipeline_object.source_catalog.reference_file_types
        else:
            pipeline_object = calwebb_spec3.Spec3Pipeline()
        
        
        # Set some parameters that pertain to the entire pipeline
        #pipeline_object.input_dir = os.getcwd()
        pipeline_object.output_file = output_name # os.path.splitext(output_file)[0]
        #pipeline_object.output_ext = ".fits" # default
        pipeline_object.save_results = True

        
        # run
        if os.path.isfile(asn_filename+'.done'):
            logger.info('Found existing {!r}, will not re-process the asn file.')
        else:
            # get exptype
            run_in_parallel = False
            nsources = 0
            asn_exptypes = ['science', 'background']
            logger.info('Examining ASN file: {}'.format(asn_filename))
            with datamodels.open(asn_filename, asn_exptypes=asn_exptypes) as models:
                exptype = models[0].meta.exposure.type # 'NRC_WFSS'
                if exptype == 'NRC_WFSS':
                    sources = [(name, model) for name, model in multislit_to_container(models).items()]
                    nsources = len(sources)
                    if len(sources) > 1:
                        if len(sources) > nparallel:
                            nparallel = len(sources)
                        if len(sources) > nparallel*nsource_per_process:
                            nsource_per_process = int(np.ceil(len(sources)/nparallel))
                        run_in_parallel = True
                models.close()
            # run in parallel
            if run_in_parallel:
                logger.info('nparallel: {}'.format(nparallel))
                logger.info('nsource_per_process: {}'.format(nsource_per_process))
                with mp.Pool(nparallel) as pool:
                    proc_list = []
                    for i in range(0, nsources, nsource_per_process):
                        index_list = np.arange(i, i+nsource_per_process).tolist()
                        proc = pool.apply_async(process_asn_async, args=(asn_filename, output_name, index_list, overwrite))
                        proc_list.append(proc)
                        time.sleep(1.0)
                    pool.close()
                    pool.join()
                    for i in range(0, nsources, nsource_per_process):
                        proc_list[i].get()
            else:
                pipeline_object.run(asn_filename)
            
            with open(asn_filename+'.done', 'w') as ofp:
                ofp.write('\n')
        
        
        # here we output source table
        if not target_acquisition_mode:
            x1d_files = sorted(glob.glob(f'{output_name}_s*_x1d.fits'))
            x1d_table = OrderedDict()
            x1d_table['source_id'] = []
            x1d_table['source_ra'] = []
            x1d_table['source_dec'] = []
            x1d_table['source_xpos'] = []
            x1d_table['source_ypos'] = []
            x1d_table['source_name'] = []
            x1d_table['source_alias'] = []
            x1d_table['x1d_fits_file'] = []
            x1d_table['x1d_ascii_file'] = []
            #with datamodels.open(cal_file) as model:
            #    [slit.source_id for slit in model.slits]
            for x1d_file in x1d_files:
                with fits.open(x1d_file) as hdul:
                    x1d_table['source_id'].append(hdul[1].header['SOURCEID'])
                    x1d_table['source_ra'].append(hdul[1].header['SRCRA'])
                    x1d_table['source_dec'].append(hdul[1].header['SRCDEC'])
                    x1d_table['source_xpos'].append(hdul[1].header['SRCXPOS'])
                    x1d_table['source_ypos'].append(hdul[1].header['SRCYPOS'])
                    if 'SRCNAME' in hdul[1].header:
                        x1d_table['source_name'].append(hdul[1].header['SRCNAME'])
                    else:
                        x1d_table['source_name'].append(hdul[1].header['SLTNAME'])
                    if 'SRCALIAS' in hdul[1].header:
                        x1d_table['source_alias'].append(hdul[1].header['SRCALIAS'])
                    else:
                        x1d_table['source_alias'].append('NULL')
                    x1d_ascii_file = os.path.splitext(x1d_file)[0] + '.dat'
                    Table(hdul[1].data).write(x1d_ascii_file, format='ascii.fixed_width', 
                                              delimiter=' ', bookend=True, overwrite=True)
                    with open(x1d_ascii_file, 'r+') as ofp:
                        ofp.seek(0)
                        ofp.write('#')
                    x1d_table['x1d_fits_file'].append(x1d_file)
                    x1d_table['x1d_ascii_file'].append(x1d_ascii_file)
            Table(x1d_table).write(output_file, overwrite=True)
            if not output_file.endswith('.dat'):
                output_file2 = os.path.splitext(output_file)[0] + '.dat'
                Table(x1d_table).write(output_file2, format='ascii.fixed_width', 
                                       delimiter=' ', bookend=True, overwrite=True)
                with open(output_file2, 'r+') as ofp:
                    ofp.seek(0)
                    ofp.write('#')
        
        
        # chdir back
        logger.info("chdir {}".format(current_dir))
        os.chdir(current_dir)
        
        
        # save pars as asdf file
        if target_acquisition_mode:
            asdf_filepath = os.path.splitext(output_filepath)[0] + '_calwebb_image3.asdf'
        else:
            asdf_filepath = os.path.splitext(output_filepath)[0] + '_calwebb_spec3.asdf'
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



