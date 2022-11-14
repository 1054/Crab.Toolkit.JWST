#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = os.path.expanduser('~/jwst_mirage_data')
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

import click
import numpy as np
import crds
import jwst
import mirage
from astropy.table import Table
from mirage.apt import apt_inputs, read_apt_xml
from mirage.utils.constants import CRDS_FILE_TYPES
from mirage.reference_files import crds_tools

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirage-fetch-crds')


DEFAULT_INSTRUMENT = 'NIRCam'
DEFAULT_FILTERS = 'F115W,F200W,F277W,F444W'
DEFAULT_DETECTORS = 'NRCA1,NRCA2,NRCA3,NRCA4,NRCA5,NRCB1,NRCB2,NRCB3,NRCB4,NRCB5'
DEFAULT_DATES = '2023-01-01'
DEFAULT_READOUT = 'FAST'





####################
### MAIN PROGRAM ###
####################

@click.command()
@click.option('--instrument', type=str, default=DEFAULT_INSTRUMENT)
@click.option('--filters', type=str, default=DEFAULT_FILTERS)
@click.option('--detectors', type=str, default=DEFAULT_DETECTORS)
@click.option('--dates', type=str, default=DEFAULT_DATES)
@click.option('--readout', type=str, default=DEFAULT_READOUT)
@click.option('--verbose', is_flag=True, default=True)
def main(
        instrument, 
        filters, 
        detectors, 
        dates,
        readout,
        verbose, 
    ):
    
    if verbose:
        logger.info('crds version: {}'.format(crds.__version__))
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    
    instrument = instrument.upper()
    if instrument not in ['NIRCAM']:
        raise NotImplementedError()
    
    
    detector_list = detectors.split(',')
    filter_list = filters.split(',')
    for detector in detector_list:
        detector = detector.upper()
        for filter_name in filter_list:
            # following "mirage/yaml/yaml_generator.py" `add_crds_reffile_names`
            pupil_name = 'CLEAR'
            if instrument == 'NIRCAM':
                exptype = 'NRC_IMAGE'
            elif instrument == 'NIRISS':
                exptype = 'NIS_IMAGE'
            elif instrument == 'FGS':
                exptype = 'FGS_IMAGE'
            elif instrument == 'MIRI':
                exptype = 'MIR_IMAGE'
                detector = 'MIRIMAGE'
            status_dict = {
                'INSTRUME': instrument, 'DETECTOR': detector,
                'FILTER': filter_name, 'PUPIL': pupil_name,
                #'READPATT': readout, 
                'EXP_TYPE': exptype,
                'DATE-OBS': dates, 
                #'TIME-OBS': time,
                'SUBARRAY': 'FULL',
            }
            if instrument == 'NIRCAM':
                if detector in ['NRCA5', 'NRCB5', 'NRCALONG', 'NRCBLONG', 'A5', 'B5']:
                    status_dict['CHANNEL'] = 'LONG'
                else:
                    status_dict['CHANNEL'] = 'SHORT'
            files_no_transmission = list(CRDS_FILE_TYPES.values())
            files_no_transmission.remove('transmission')
            reffiles = crds_tools.get_reffiles(
                status_dict, 
                files_no_transmission,
                download = True,
            )
            # crds.getreferences(
            #     status_dict, 
            #     ['mask', 'distortion', 'ipc']
            # )




if __name__ == '__main__':
    main()



