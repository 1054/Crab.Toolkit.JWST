#!/usr/bin/env python
#
"""
Subtract background for a cal image using jwst.skymatch. 

The input cal image file will be updated inplace, unless --no-inplace is set!

Inputs:
    
    *_cal.fits

Outputs:
    
    *_cal.fits (inplace)
    *_cal_skymatchstep.fits

"""

# Packages that allow us to get information about objects:
import os, sys, re, json, copy, datetime, time, glob, shutil
import asdf
import click
from collections import OrderedDict
try:
    import packaging.version
    LooseVersion = packaging.version.Version
except:
    from distutils.version import LooseVersion
    # DeprecationWarning: distutils Version classes are deprecated. Use packaging.version instead.

# Set CRDS
if not ('CRDS_PATH' in os.environ):
    os.environ['CRDS_PATH'] = os.path.expanduser('~/jwst_crds_cache')
if not ('CRDS_SERVER_URL' in os.environ):
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'
#if not ('CRDS_CONTEXT' in os.environ):
#    os.environ['CRDS_CONTEXT'] = 

# Numpy library:
import numpy as np

# Astropy tools:
from astropy.io import fits

# Import JWST pipeline-related modules

# The entire calwebb_image2 pipeline
from jwst.pipeline import calwebb_image2

# Individual steps that make up calwebb_image2
from jwst.background import BackgroundStep
from jwst.assign_wcs import AssignWcsStep
from jwst.flatfield import FlatFieldStep
from jwst.photom import PhotomStep
from jwst.resample import ResampleStep
from jwst import datamodels
from jwst.datamodels import dqflags

# importing an individual pipeline step
from jwst.skymatch import SkyMatchStep

# Import jwst package itself
import jwst

# Import stdatamodels package
import stdatamodels

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
@click.argument('input_cal_file', type=click.Path(exists=True))
@click.option('--inplace/--no-inplace', is_flag=True, default=True)
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
def main(
        input_cal_file, 
        inplace, 
        overwrite, 
    ):
    
    # Add script dir to sys path
    if not (get_script_dir() in sys.path):
        sys.path.append(get_script_dir())
    
    # Setup logger
    logger = setup_logger()
    
    # Print JWST pipeline version
    logger.info('JWST pipeline version: {}'.format(jwst.__version__))
    
    # check CRDS 
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
    
    # Set input_filepath, skymatchstep_filepath
    input_filepath = input_cal_file
    input_filename = os.path.splitext(os.path.basename(input_filepath))[0] # no .fits suffix, has _cal suffix.
    work_dir = os.path.dirname(input_filepath)
    backup_filepath = os.path.join(work_dir, input_filename + '_before_skymatchstep.fits')
    skymatchstep_filepath = os.path.join(work_dir, input_filename + '_skymatchstep.fits')
    skymatchstep_filename = input_filename + '_skymatchstep'
    if not os.path.isfile(backup_filepath): #<TODO># 
        shutil.copy2(input_filepath, backup_filepath) #<TODO># 
    if inplace: 
        # Check input_filepath history
        with datamodels.open(input_filepath) as model:
            #if model.meta.instrument.name.upper() != 'NIRCAM': #20220914 do this also for MIRI?
            #    logger.warning('The input data is not a NIRCam imaging data. Will not run background subtraction.')
            #    return
            for entry in model.history:
                for k,v in entry.items():
                    if 'Done background subtraction; go-jwst-imaging-stage-2-step-2-do-bkgsub.py' in v: 
                        if os.path.isfile(skymatchstep_filepath):
                            logger.warning(f'{input_filepath} already had background subtraction. Skipping!')
                            return
                        else:
                            break
        if os.path.isfile(skymatchstep_filepath):
            shutil.move(skymatchstep_filepath, skymatchstep_filepath+'.backup')
    else:
        # Check skymatchstep_filepath existence
        if os.path.isfile(skymatchstep_filepath):
            if overwrite:
                shutil.move(skymatchstep_filepath, skymatchstep_filepath+'.backup')
            else:
                logger.info('Found existing output file {!r} and overwrite is False. Do nothing.'.format(
                            skymatchstep_filepath))
                return
    
    
    # Print progress
    logger.info("Processing {} -> {}".format(input_filepath, skymatchstep_filepath))
    
    
    # following CEERS, 
    # do sky subtraction at stage 2, 
    # this needs an association file
    
    
    # read program info
    header = fits.getheader(input_filepath, 0)
    program = header['PROGRAM'].strip()
    obs_id = header['OBSERVTN'].strip()
    if 'TARGPROP' in header:
        target_group = header['TARGPROP'].strip()
    else:
        target_group = 'TARGET'
    if 'TARGNAME' in header:
        target_name = header['TARGNAME'].strip()
    else:
        target_name = ''
    if target_name == '':
        target_name = target_group
    
    
    # prepare association content
    asn_dict = OrderedDict()
    asn_dict['asn_type'] = 'None'
    asn_dict['asn_rule'] = 'DMS_Level3_Base'
    asn_dict['version_id'] = None
    asn_dict['code_version'] = jwst.__version__
    asn_dict['degraded_status'] = 'No known degraded exposures in association.'
    asn_dict['program'] = program # 'noprogram'
    asn_dict['constraints'] = 'No constraints'
    asn_dict['asn_id'] = obs_id
    asn_dict['target'] = target_name
    asn_dict['asn_pool'] = 'none'
    asn_dict['products'] = []
    product_dict = OrderedDict()
    product_dict['name'] = skymatchstep_filename # input_filename
    product_dict['members'] = [
            {'expname': os.path.abspath(input_filepath), # os.path.relpath(input_filepath, work_dir), # note that work_dir is just dirname(input_filepath)
             'exptype': 'science',
            }
        ]
    asn_dict['products'].append(product_dict)
    
    asn_filepath = os.path.join(work_dir, skymatchstep_filename + '_asn.json')
    
    if os.path.isfile(asn_filepath):
        shutil.move(asn_filepath, asn_filepath+'.backup')
    
    with open(asn_filepath, 'w') as fp:
        json.dump(asn_dict, fp, indent=4)
    
    
    # Set parameters for SkyMatchStep
    skymatch = SkyMatchStep()
    skymatch.save_results = True
    skymatch.output_dir = work_dir
    skymatch.output_file = input_filename # SkyMatchStep will append 'skymatchstep' to the input filename if output_file is undefined.
    skymatch.output_ext = ".fits"
    skymatch.suffix = 'cal_skymatchstep' # so that the output is "*_cal_skymatchstep.fits". In default the code will remove "_cal".
    # set sky statistics parameters
    skymatch.skymethod = "local" # the default is global+match, doesn't matter as we're processing files individually
    skymatch.lsigma = 2.0
    skymatch.usigma = 2.0
    skymatch.nclip = 10
    #skymatch.upper = 1.0 # for NIRCam this is okay, but not for miri
    # set the 'subtract' parameter so the calculated sky value is removed from the image
    # (subtracting the calculated sky value from the image is off by default)
    skymatch.subtract = True
    
    
    # Run SkyMatchStep
    sky = skymatch.run(asn_filepath)
    
    
    # check bkgsub output
    print('Checking {!r}'.format(skymatchstep_filepath))
    if not os.path.isfile(skymatchstep_filepath):
        skymatchstep_filepath2 = skymatchstep_filepath.replace('cal_skymatchstep.fits', '0_cal_skymatchstep.fits') # 20241217
        if os.path.isfile(skymatchstep_filepath2):
            shutil.copy2(skymatchstep_filepath2, skymatchstep_filepath)
    assert os.path.isfile(skymatchstep_filepath)
    
    
    # add history entry
    with datamodels.open(skymatchstep_filepath) as model:
        timeobj = datetime.datetime.now()
        historystr = 'Done background subtraction; go-jwst-imaging-stage-2-step-2-do-bkgsub.py' + \
                     '; time = ' + timeobj.strftime('%Y-%m-%d %H:%M:%S')
        historyobj = stdatamodels.util.create_history_entry(historystr)
        model.history.append(historyobj)
        
        # for jwst.__version__ <= 1.6.2, 
        # the skymatch module is not actually saving the background-subtracted image, 
        # see "jwst/skymatch/skymatch_step.py" `process()`, which converts `img.models_grouped` into `images`
        # but `images` are not saved back into `img`. The `process()` returns `img`, not `images`, so
        # only FITS headers are updated, but the background-subtracted image are not saved to disk.
        #if LooseVersion(jwst.__version__) <= LooseVersion("1.6.2"):
        
        model.save(
            skymatchstep_filepath,
            overwrite=True,
        )
    
    
    # for jwst.__version__ <= 1.6.2, 
    # the skymatch module is not actually saving the background-subtracted image, 
    # see "jwst/skymatch/skymatch_step.py" `process()`, which converts `img.models_grouped` into `images`
    # but `images` are not saved back into `img`. The `process()` returns `img`, not `images`, so
    # only FITS headers are updated, but the background-subtracted image are not saved to disk.
    if LooseVersion(jwst.__version__) <= LooseVersion("1.6.2"):
        
        # output to "_cal.fits"
        with fits.open(skymatchstep_filepath) as hdul:
            sky = hdul[0].header['BKGLEVEL']
            assert (hdul[1].header['EXTNAME'] == 'SCI')
            hdul[0].header['HISTORY'] = 'Subtracted background level {} in the SCI image; {}'.format(
                sky, datetime.datetime.now().strftime("%Y:%m:%d %Hh%Mm%Ss") + time.tzname[time.daylight])
            hdul[1].data -= sky
            hdul.writeto(skymatchstep_filepath)
    
    
    # update "_cal.fits"
    # directly copy "_skymatchstep.fits" to "_cal.fits"
    if inplace:
        logger.info('Updating {} inplace'.format(input_filepath))
        shutil.copy2(skymatchstep_filepath, input_filepath)
    
    
    # log
    logger.info("Processed {} -> {} -> {}".format(input_filepath, skymatchstep_filepath, input_filepath))




# Entry
if __name__ == '__main__':
    main()



