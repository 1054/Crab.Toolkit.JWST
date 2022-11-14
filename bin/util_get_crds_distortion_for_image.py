#!/usr/bin/env python
# 
"""
Finding source emission in the image, obtain 2D background.

Usage: 
    ./util_make_seed_image_for_rate_image.py ngc0628_miri_lvl3_f770w_i2d_hackedbgr.fits

By Daizhong Liu @MPE. 

Last update: 2022-07-27.
Last update: 2022-10-03 rewritten based on "util_make_seed_image_for_rate_image.py".
"""
import os, sys, re, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.modeling import models as apy_models
from astropy.modeling import fitting as apy_fitting
from astropy.convolution import convolve as apy_convolve
from astropy.convolution import Gaussian2DKernel
from photutils.background import Background2D
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# mpl.rcParams['savefig.dpi'] = 300
# mpl.rcParams['figure.dpi'] = 300

if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

import jwst
from jwst import datamodels

import asdf

# logging
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# version
VERSION = '20221003'



def prep_crds_getreferences_kwargs(model):
    crds_getreferences_kwargs = {}
    
    import crds
    crds_context = None
    if 'CRDS_CONTEXT' in os.environ: 
        if os.environ['CRDS_CONTEXT'] != '':
            crds_context = os.environ['CRDS_CONTEXT']
    if crds_context is None:
        crds_context = crds.get_default_context()
    
    # CRDS query parameters
    # https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/reference_files.html
    instrument_name = model.meta.instrument.name
    if instrument_name.upper() in ['NIRCAM', 'NIRISS']: 
        crds_dict = {'INSTRUME':instrument_name, 
                     'DETECTOR':model.meta.instrument.detector, 
                     'CHANNEL':model.meta.instrument.channel, 
                     'FILTER':model.meta.instrument.filter, 
                     'PUPIL':model.meta.instrument.pupil, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    elif instrument_name.upper() == 'MIRI': 
        crds_dict = {'INSTRUME':instrument_name, 
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'BAND':model.meta.instrument.band, 
                     'READPATT':model.meta.exposure.readpatt, 
                     'SUBARRAY':model.meta.subarray.name, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    else:
        raise ValueError("Error! The input model (model.meta.instrument.name) is not NIRCAM or MIRI!")
    
    crds_getreferences_kwargs['parameters'] = crds_dict
    crds_getreferences_kwargs['reftypes'] = ['distortion']
    crds_getreferences_kwargs['context'] = crds_context
    
    return crds_getreferences_kwargs



def get_distortion_from_CRDS(model, raise_exception=True):
    logger.info('Querying flat from CRDS')
    import crds
    
    if isinstance(model, str):
        with datamodels.open(model) as model:
            crds_getreferences_kwargs = prep_crds_getreferences_kwargs(model)
    else:
        crds_getreferences_kwargs = prep_crds_getreferences_kwargs(model)
    logger.info('crds_getreferences_kwargs: {}'.format(crds_getreferences_kwargs))
    ret = crds.getreferences(**crds_getreferences_kwargs)
    
    # check if CRDS got the flat correctly
    try:
        retfile = ret['distortion']
    except KeyError:
        errmsg = 'Flat was not found in CRDS with the parameters: {}'.format(crds_dict)
        logger.error(errmsg)
        if raise_exception:
            raise Exception(errmsg)
        else:
            return None
    
    return retfile
    






@click.command()
@click.argument('fits_image', type=click.Path(exists=True))
def main(
        fits_image, 
    ):
    
    distortionfile = get_distortion_from_CRDS(
        fits_image, 
    )
    print('Got distortion file: {!r}'.format(distortionfile))

    with asdf.open(distortionfile) as dist_file:
        coord_transform = dist_file.tree['model']
    print('coord_transform', coord_transform)


    # from mirage "catalog_seed_image.py"
    # from mirage.utils import siaf_interface
    # siaf_inst = 'NIRCam'
    # instrument_siaf = siaf_interface.get_instance(siaf_inst)
    # self.siaf = instrument_siaf[self.params['Readout']['array_name']]
    # self.local_roll, self.attitude_matrix, self.ffsize, \
    #     self.subarray_bounds = siaf_interface.get_siaf_information(instrument_siaf,
    #                                                                self.params['Readout']['array_name'],
    #                                                                self.ra, self.dec,
    #                                                                self.params['Telescope']['rotation'])
    
    # from mirage.seed_image.catalog_seed_images import Catalog_seed
    # seed = Catalog_seed()




if __name__ == '__main__':
    
    main()



