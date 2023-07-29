#!/usr/bin/env python 
# 

"""
Tune the WCS for specific data sets. 

We need a json file under the `CODE_DATADIR` directory recording the correct CRPIX and CRVAL values. 
The json file name should start with the data set file name, e.g., 
'jw01837004020_02201_00001_nrca1_new_wcs.json'. 

"""

import os, sys, re, copy, json, shutil, glob, datetime
import click
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, FK5
from astropy.wcs import WCS
from astropy.modeling.models import Shift

# jwst
from jwst import datamodels
from jwst.datamodels import ImageModel
from jwst.assign_wcs.util import update_fits_wcsinfo
from stdatamodels import util


# code name and version
CODE_NAME = 'util_correct_astrometry_for_specific_data_sets.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230107'
CODE_HOMEPAGE = 'https://github.com/1054/Crab.Toolkit.JWST'
CODE_DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                            'util_correct_astrometry_for_specific_data_sets_wcs_files')

# logging
import logging
logging.basicConfig(level='INFO')
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--template-dir', type=click.Path(exists=True), default=CODE_DATADIR, help='Where to find the custom WCS json dictionary files.')
def main(
        input_file, 
        output_file, 
        template_dir, 
    ):
    """
    Input a cal image (which should already have a WCS assigned by the calibration pipeline). 
    """
    
    # find wcs json file(s)
    dataset_name = os.path.splitext(os.path.basename(input_file))[0]
    dataset_name = re.sub(r'^(.*)(_rate|_cal)$', r'\1', dataset_name)
    obsvisit_prefix = dataset_name.split('_')[0] # e.g. jw01837004020
    wcs_json_files = glob.glob(os.path.join(template_dir, dataset_name + '*.json'))
    if len(wcs_json_files) == 0:
        wcs_json_files = glob.glob(os.path.join(template_dir, '*', dataset_name + '*.json'))
    if len(wcs_json_files) == 0:
        wcs_json_files = glob.glob(os.path.join(template_dir, obsvisit_prefix + '*.json'))
    if len(wcs_json_files) == 0:
        wcs_json_files = glob.glob(os.path.join(template_dir, '*', obsvisit_prefix + '*.json'))
    
    if len(wcs_json_files) == 0:
        logger.warning('No WCS json file is found associating with the data set {!r} or obs-visit prefix {!r} in the directory {!r}. Skipping.'.format(
            dataset_name, obsvisit_prefix, template_dir, 
        ))
        return
    
    wcs_json_file = None
    if len(wcs_json_files) > 0:
        if len(wcs_json_files) > 1:
            wcs_json_files = sorted(wcs_json_files)
            wcs_json_file = wcs_json_files[-1]
            logger.warning('Found multiple WCS json file(s): {}'.format(wcs_json_files))
            logger.warning('Taking the last WCS json file: {}'.format(wcs_json_file))
        else:
            wcs_json_file = wcs_json_files[0]
            logger.info('Found WCS json file: {}'.format(wcs_json_file))
    
    
    # read WCS json dict
    new_wcs_params = {}
    with open(wcs_json_file, 'r') as fp:
        new_wcs_params = json.load(fp)
    
    
    # read input image
    logger.info('Reading input data: {!r}'.format(input_file))
    with datamodels.open(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Corrected astrometry'):
                    logger.info('{!r} already had astrometry corrected. Skipping!'.format(input_file))
                    return
        
    
    
    # output file
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_correcting_astrometry.fits'
        if not os.path.isfile(output_orig):
            shutil.copy2(output_file, output_orig) # do not overwrite
    else:
        if os.path.isfile(output_file):
            shutil.move(output_file, output_file+'.backup')
        elif os.path.dirname(output_file) != '' and not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        shutil.copy2(input_file, output_file) # copy input to output then update the fits header
    
    out_timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    old_wcs_params = {}
    with fits.open(output_file, mode='update') as hdul:
        if 'old_coordinate' in new_wcs_params and 'new_coordinate' in new_wcs_params:
            if isinstance(new_wcs_params['old_coordinate'], (list,tuple)):
                old_ra_str, old_dec_str = map(str, new_wcs_params['old_coordinate'])
            else:
                old_ra_str, old_dec_str = new_wcs_params['old_coordinate'].split()
            if re.match(r'^[0-9.]+$', old_ra_str) and re.match(r'^([+-]?)[0-9.]+$', old_dec_str):
                old_coordinate = SkyCoord(float(old_ra_str), float(old_dec_str), unit=(u.deg, u.deg), frame=FK5)
            else:
                old_coordinate = SkyCoord(old_ra_str+' '+old_dec_str, unit=(u.hour, u.deg), frame=FK5)
            # 
            if isinstance(new_wcs_params['new_coordinate'], (list,tuple)):
                new_ra_str, new_dec_str = map(str, new_wcs_params['new_coordinate'])
            else:
                new_ra_str, new_dec_str = new_wcs_params['new_coordinate'].split()
            if re.match(r'^[0-9.]+$', new_ra_str) and re.match(r'^([+-]?)[0-9.]+$', new_dec_str):
                new_coordinate = SkyCoord(float(new_ra_str), float(new_dec_str), unit=(u.deg, u.deg), frame=FK5)
            else:
                new_coordinate = SkyCoord(new_ra_str+' '+new_dec_str, unit=(u.hour, u.deg), frame=FK5)
            # 
            logger.info('Correcting astrometry with old_coordinate: {} -> new_coordinate: {}'.format(
                old_coordinate.to_string('hmsdms', sep=':', precision=3), 
                new_coordinate.to_string('hmsdms', sep=':', precision=3)))
            new_wcs_params['CRVAL1'] = hdul['SCI'].header['CRVAL1'] - old_coordinate.ra.deg + new_coordinate.ra.deg
            new_wcs_params['CRVAL2'] = hdul['SCI'].header['CRVAL2'] - old_coordinate.dec.deg + new_coordinate.dec.deg
        # 
        for key in new_wcs_params:
            if key in hdul['SCI'].header:
                old_wcs_params[key] = hdul['SCI'].header[key]
                hdul['SCI'].header[key] = new_wcs_params[key]
                logger.info('Updating header key: {}, value {} -> {}'.format(key, old_wcs_params[key], new_wcs_params[key]))
        hdul.flush()
    
    with ImageModel(output_file) as dm:
        
        #from jwst.assign_wcs.util import update_fits_wcsinfo
        #update_fits_wcsinfo(
        #    dm,
        #    max_pix_error=1., # default is 0.01 in "tweakreg_step.py"
        #    npoints=16
        #)
        
        dra = new_wcs_params['CRVAL1'] - old_wcs_params['CRVAL1'] # equiv. dm.meta.wcsinfo.crval1 - dm.meta.wcs.to_fits_sip()['CRVAL1']
        ddec = new_wcs_params['CRVAL2'] - old_wcs_params['CRVAL2'] # equiv. dm.meta.wcsinfo.crval2 - dm.meta.wcs.to_fits_sip()['CRVAL2']
        
        # dm.meta.wcs.available_frames : ['detector', 'v2v3', 'v2v3vacorr', 'world']
        # set_transform needs from_frame and to_frame in in sequence, so we only update the last sequence
        trans = dm.meta.wcs.get_transform('v2v3vacorr', 'world')
        trans2 = trans | (Shift(dra) & Shift(ddec))
        dm.meta.wcs.set_transform('v2v3vacorr', 'world', trans2)
        update_fits_wcsinfo(dm)
        
        # add history entry following CEERS
        stepdescription = f'Corrected astrometry with {CODE_NAME} ({out_timestamp})'
        software_dict = {'name':CODE_NAME,
                         'author':CODE_AUTHOR,
                         'version':CODE_VERSION,
                         'homepage':CODE_HOMEPAGE}
        substr = util.create_history_entry(stepdescription, software=software_dict)
        dm.history.append(substr)
        #print('dm.history', dm.history)
        dm.save(output_file)
        logger.info('Saved astrometry-corrected data into {!r}'.format(output_file))
    
    
    # save wcs param changes
    output_key_change_file = re.sub(r'\.fits$', r'', output_file) + '_correcting_astrometry_wcs_param_changes.json'
    if os.path.isfile(output_key_change_file):
        shutil.move(output_key_change_file, output_key_change_file+'.backup')
    with open(output_key_change_file, 'w') as fp:
        json.dump({'old':old_wcs_params, 'new':new_wcs_params}, fp, indent=4)
    
    



if __name__ == '__main__':
    main()



