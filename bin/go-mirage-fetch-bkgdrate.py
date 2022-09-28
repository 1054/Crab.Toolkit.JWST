#!/usr/bin/env python
# 
import os, sys, re, json, shutil
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = '/n17data/dzliu/Data/jwst_crds_cache'
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = '/n23data1/hjmcc/jwst/mirage/mirage_data'
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

import astropy.units as u
import click
import numpy as np
import crds
import jwst
import mirage
from astropy.coordinates import SkyCoord, FK5
from mirage.utils import siaf_interface
from mirage.utils.backgrounds import calculate_background, day_of_year_background_spectrum
from mirage.utils.constants import CRDS_FILE_TYPES, MEAN_GAIN_VALUES
from mirage.utils.utils import get_filter_throughput_file

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirage-fetch-bkgdrate')


DEFAULT_INSTRUMENT = 'NIRCam'
DEFAULT_FILTERS = 'F115W,F200W,F277W,F444W'
DEFAULT_DETECTORS = 'NRCA1,NRCA2,NRCA3,NRCA4,NRCA5,NRCB1,NRCB2,NRCB3,NRCB4,NRCB5'
DEFAULT_BKGLEVELS = 'low,medium,high'
DEFAULT_OUTPUT_DIR = './'
DEFAULT_OUTPUT_NAME = 'bkgdrate'





####################
### MAIN PROGRAM ###
####################

@click.command()
@click.argument('ra', type=str)
@click.argument('dec', type=str)
@click.option('--instrument', type=str, default=DEFAULT_INSTRUMENT)
@click.option('--filters', type=str, default=DEFAULT_FILTERS)
@click.option('--detectors', type=str, default=DEFAULT_DETECTORS)
#@click.option('--dates', type=str, default=None)
@click.option('--background-levels', 'bkglevels', type=str, default=DEFAULT_BKGLEVELS)
@click.option('--output-dir', type=click.Path(exists=False), default=DEFAULT_OUTPUT_DIR)
@click.option('--output-name', type=str, default=DEFAULT_OUTPUT_NAME)
@click.option('--update', is_flag=True, default=True)
@click.option('--requery', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        ra, 
        dec, 
        instrument, 
        filters, 
        detectors,
        #dates, 
        bkglevels, 
        output_dir,
        output_name,
        update, 
        requery, 
        verbose, 
    ):

    if verbose:
        logger.info('crds version: {}'.format(crds.__version__))
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    if output_name.endswith('.json'):
        output_name = re.sub(r'\.json$', r'', output_name)
    out_file = os.path.join(output_dir, output_name + '.json')
    if os.path.isfile(out_file):
        if not update:
            logger.info('Found existing file {!r} and update is set to False. Returning!'.format(out_file))
            return
    
    out_dict = {}
    if os.path.isfile(out_file):
        with open(out_file, 'r') as fp:
            try:
                out_dict = json.load(fp)
            except:
                pass
    
    
    if re.match(r'^[0-9.+-]+$', ra) and re.match(r'^[0-9.+-]+$', dec):
        ra, dec = float(ra), float(dec)
    else:
        scoord = SkyCoord(ra, dec, unit=(u.hour, u.deg), frame=FK5)
        ra, dec = scoord.ra.deg, scoord.dec.deg
    
    instrument = instrument.lower()
    if instrument not in ['nircam']:
        raise NotImplementedError()
    
    instrument_siaf = siaf_interface.get_instance('NIRCam')
    
    has_update = False
    
    detector_list = detectors.split(',')
    filter_list = filters.split(',')
    for detector in detector_list:
        detector = detector.lower()
        for filter_name in filter_list:
            module = detector[3]
            longwave = (detector[4]=='5')
            gain_value = MEAN_GAIN_VALUES[instrument][detector]
            filter_file = get_filter_throughput_file(
                instrument,
                filter_name,
                'CLEAR',
                fgs_detector=detector, 
                nircam_module=module,
            )
            
            siaf = instrument_siaf[detector.upper()+'_FULL'] # NRCA1_FULL
            
            # if dates is not None:
            #     bkgd_wave, bkgd_spec = day_of_year_background_spectrum(
            #         ra,
            #         dec,
            #         dates
            #     )
            #     bkgdrate = calculate_background(
            #         ra, 
            #         dec,
            #         filter_file, 
            #         True,
            #         gain_value, 
            #         siaf,
            #         back_wave = bkgd_wave,
            #         back_sig = bkgd_spec
            #     )
            #     out_dict[] = bkgdrate
            # else:
            for bkglevel in bkglevels.split(','):
                key = '{}_{}_{}_{:.8f}_{:.8f}_{}'.format(instrument, detector, filter_name, ra, dec, bkglevel)
                if key not in out_dict or requery:
                    if verbose:
                        logger.info('Querying bkgdrate for {}'.format(key))
                    bkgdrate = calculate_background(
                        ra, 
                        dec,
                        filter_file, 
                        False,
                        gain_value, 
                        siaf,
                        level = bkglevel,
                    )
                    out_dict[key] = bkgdrate
                    has_update = True
    
    if os.path.isfile(out_file):
        shutil.copy2(out_file, out_file+'.backup')
    with open(out_file, 'w') as fp:
        json.dump(out_dict, fp, indent=4)
    
    if verbose:
        if has_update:
            logger.info('Updated bkgdrate file {!r}'.format(out_file))
        else:
            logger.info('No update to the existing bkgdrate file {!r}'.format(out_file))




if __name__ == '__main__':
    main()



