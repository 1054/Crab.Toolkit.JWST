#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id, get science target uncal fits files. 

Rate, cal and i2d files are excluded unless specified.

The default download dir is 
./mastDownload/JWST/jwpppppooovvv_ggsaa_eeeee_detector/jwpppppooovvv_ggsaa_eeeee_detector_prodType.fits

"""

import os, sys, re, datetime
import click
from collections import namedtuple
from astroquery.mast import Observations

JWST_Dataset_Name = namedtuple('JWST_Dataset_Name', 
                               ['proposal_id', 'obs_num', 'visit_num', 
                                'visit_group', 'parallel', 'activity', 
                                'exposure', 'detector', 'prod_type'])

def parse_jwst_dataset_name(input_str, raise_exception=True):
    regex_format = r'^jw([0-9]{5})([0-9]{3})([0-9]{3})_([0-9]{2})([0-9]{1})([0-9]{2})_([0-9]{5})_([a-z0-9]+)(.*)$'
    regex_match = re.match(regex_format, input_str)
    if regex_match is not None:
        return JWST_Dataset_Name(*regex_match.groups())
    else:
        if raise_exception:
            raise Exception('Error! The input prefix does not seem to have the right format: {}'.format(regex_format))
        return None

def get_instrument_from_detector_name(detector_name):
    instrument_detector_name = {
        'NIRCam': ['nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcalong', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4', 'nrcblong'], 
        'MIRI': ['mirimage'],
    }
    for instrument in instrument_detector_name:
        if detector_name in instrument_detector_name[instrument]:
            return instrument
    return None


def get_proposal_id_from_survey_name(survey_name):
    survey_proposal_id_dict = {}
    survey_proposal_id_dict['CEERS'] = '01345'
    survey_proposal_id_dict['PRIMER'] = '01837'
    survey_proposal_id_dict['COSMOS'] = '01727'
    survey_proposal_id_dict['PHANGS'] = '02107'
    if survey_name in survey_proposal_id_dict:
        return survey_proposal_id_dict[survey_name]
    return None


# Setup logging
import logging
from astroquery import log as apy_log

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
    logger_streamhandler.setLevel(logging.INFO)

    log_file = get_script_name()
    log_time = datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    log_filepath = f"log_{log_file}_{log_time}.txt"
    logger_filehandler = logging.FileHandler(log_filepath, mode='a')
    logger_filehandler_formatter = logging.Formatter("[%(asctime)-15s] %(message)s", "%Y-%m-%d %H:%M:%S")
    logger_filehandler.setFormatter(logger_filehandler_formatter)

    logger = logging.getLogger()
    while len(logger.handlers) > 0:
        del logger.handlers[0]
    logger.addHandler(logger_streamhandler)
    logger.addHandler(logger_filehandler)
    logger.setLevel(logging.INFO)
    
    return logger, log_filepath



@click.command()
@click.argument('program', type=str)
@click.option('--obs-num', '--obsnum', '--obs', type=str, default=None, help='Observation number as listed in the program information webpage.')
@click.option('--calib-level', type=str, default='1', help='Selecting calib level. "1" for uncal data, "2" for rate data, "3" for drizzled image data. Can be multiple like "1,2,3".')
@click.option('--extension', type=str, default='fits', help='Selecting extension. Usually just "fits". Can be "fits,json,jpg" if you want to get those files too.')
@click.option('--preview', is_flag=True, default=False, help='Selecting preview files, e.g., i2d quicklook image in jpg format.')
@click.option('--auxiliary', is_flag=True, default=False, help='Selecting auxiliary files, e.g., guide star data.')
@click.option('--download/--no-download', is_flag=True, default=False)
@click.option('--download-dir', type=click.Path(exists=False), default='.')
def main(
        program, 
        obs_num, 
        calib_level, 
        extension, 
        preview, 
        auxiliary, 
        download, 
        download_dir, 
    ):

    # setup logger
    logger, log_file = setup_logger()
    logger.info('log file: {}'.format(log_file))
    
    # print user input program
    logger.info('program: {}'.format(program))
    
    # print user input download
    logger.info('download: {}'.format(download))
    
    # check if program is a proposal ID or a known survey name
    if re.match(r'^[0-9]+$', program):
        proposal_id = '{:05d}'.format(int(program))
    else:
        proposal_id = get_proposal_id_from_survey_name(program)
    
    if proposal_id is None:
        logger.error('Error! Cannot understand the input program ID. Example: 01345')
        sys.exit(255)
    
    # print proposal id to query
    logger.info('proposal_id: {}'.format(proposal_id))
    
    # check if user has specified an obs num
    if obs_num is not None:
        if obs_num.startswith('obs'):
            obs_num = re.sub(r'^obs', r'', obs_num)
        
        if re.match(r'^[0-9]+$', obs_num):
            obs_num = '{:03d}'.format(int(obs_num))
        elif re.match(r'^[0-9]+(,[0-9]+)+$', obs_num):
            obs_num = '{:03d}'.format(int(obs_num))
        else:
            logger.error('Error! Cannot understand the input obs_num. Example: 001')
            sys.exit(255)
    
        logger.info('obs_num: {}'.format(obs_num))
    
    # query the MAST databse
    logger.info('Running Observations.query_criteria')
    obs_list = Observations.query_criteria(
        obs_collection = "JWST",
        proposal_id = proposal_id,
        calib_level = [3], # cannot directly find level 1 products, must search level3 first then call get_product_list
    )
    logger.info('Returned {} rows'.format(len(obs_list)))
    
    # parse user input calib_level, extension, etc. for later filtering the product list
    if calib_level.find(',')>=0:
        calib_level = [int(t) for t in calib_level.replace(' ','').split(',')]
    else:
        calib_level = [int(calib_level)]
    logger.info('calib_level: {}'.format(calib_level))
    
    if extension.find(',')>=0:
        extension = extension.replace(' ','').split(',')
    else:
        extension = [extension]
    logger.info('extension: {}'.format(extension))
    
    productType = ["SCIENCE"]
    if preview:
        productType.append("PREVIEW")
        if 'jpg' not in extension:
            extension.append('jpg')
    if auxiliary:
        productType.append("AUXILIARY")
    logger.info('productType: {}'.format(productType))
    
    # loop obs list
    for iobs, obs in enumerate(obs_list):
        #logger.info('obs_id: {}'.format(obs['obs_id']))
        is_obs_num_matched = True
        if obs_num is not None:
            is_obs_num_matched = False
            if isinstance(obs_num, str):
                obs_num_list = [obs_num]
            else:
                obs_num_list = obs_num
            for obs_num_str in obs_num_list:
                check_obs_id = 'jw{}-o{}_t[0-9]+_.*_.*'.format(
                    proposal_id, obs_num_str)
                if re.match(check_obs_id, obs['obs_id']) is not None:
                    is_obs_num_matched = True
                    break
        if not is_obs_num_matched:
            continue
        logger.info('*** --- ({}/{}) "{}" --- ***'.format(iobs+1, len(obs_list), obs['obs_id']))
        product_list = Observations.get_product_list(obs)
        if len(product_list) > 0:
            #print('product_list', product_list)
            products = Observations.filter_products(
                product_list,
                calib_level = calib_level,
                productType = productType,
                extension = extension,
            )
            if len(products) > 0: 
                for iproduct, product in enumerate(products):
                    logger.info('        ({}/{}) {} {}'.format(
                        iproduct+1, len(products), 
                        product['obs_id'], 
                        product['productFilename']
                    ))
                # 
                if download:
                    manifest = Observations.download_products(
                        products, 
                        download_dir = download_dir, 
                    )
                    logger.info('        checking downloads ...')
                    for iproduct, product in enumerate(products):
                        downloaded_filepath = '{}/mastDownload/JWST/{}/{}'.format(
                            download_dir, 
                            product['obs_id'], 
                            product['productFilename']
                        )
                        if os.path.isfile(downloaded_filepath):
                            downloaded_str = 'Downloaded'
                        else:
                            downloaded_str = 'Failed to download'
                        logger.info('        ({}/{}) {} {}'.format(
                            iproduct+1, len(products), 
                            downloaded_str, 
                            downloaded_filepath
                        ))
                    
            else:
                logger.info('No product after filtering calib_level {} productType {} extension {}!'.format(
                    calib_level, productType, extension
                    ))
        else:
            logger.info('No product list in obs "{}"'.format(obs['obs_id']))



# 
if __name__ == '__main__':
    main()



