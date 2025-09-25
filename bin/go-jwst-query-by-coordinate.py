#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id, get science target uncal fits files. 

Rate, cal and i2d files are excluded unless specified.

The default download dir is 
./mastDownload/JWST/jwpppppooovvv_ggsaa_eeeee_detector/jwpppppooovvv_ggsaa_eeeee_detector_prodType.fits

"""

import os, sys, shutil, re, datetime
import numpy as np
import click
import astropy.units as u
from collections import namedtuple
from astropy.coordinates import SkyCoord
from astropy.table import Table
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

def check_guide_star_data(product):
    if re.match(r'^.*_gs-(fg|track|acq[0-9]|id)_.*$', product['productFilename']):
        return True
    return False



class SkyCoordParamType(click.ParamType):
    name = "SkyCoord"

    def convert(self, value, param, ctx):
        scoord = None
        if isinstance(value, SkyCoord):
            scoord = value
        
        elif isinstance(value, str):
            if value.find(':') >= 0 or value.find('h') >= 0:
                try:
                    scoord = SkyCoord(value, unit=(u.hourangle, u.deg))
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
            else:
                try:
                    scoord = SkyCoord(value, unit=(u.deg, u.deg))
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
        
        elif isinstance(value, (tuple, list)) and len(value) == 2:
            if str(value[0]).find(':') >= 0 or str(value[0]).find('h') >= 0:
                try:
                    scoord = SkyCoord(value[0], value[1], unit=(u.hourangle, u.deg))
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
            else:
                try:
                    scoord = SkyCoord(value[0], value[1], unit=(u.deg, u.deg))
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
        
        return scoord

SkyCoordParamType_ = SkyCoordParamType()


class SkyCoordArgument(click.Argument):
    
    def type_cast_value(self, ctx, value) -> list:
        #print('SkyCoordArgument type_cast_value')
        try:
            #print('SkyCoordArgument type_cast_value value {} type {}'.format(value, type(value)))
            assert len(value)%2==0
            ra_dec_list = []
            for ra, dec in zip(value[0::2], value[1::2]):
                ra_dec_list.append(SkyCoordParamType_.convert((ra, dec), param=self, ctx=ctx))
            #print('SkyCoordArgument type_cast_value ra_dec_list {}'.format(ra_dec_list))
            return ra_dec_list
        except Exception:
            raise click.BadParameter(value)


class SkyCoordOption(click.Option):
    
    def type_cast_value(self, ctx, value) -> list:
        if value is None:
            return None
        #print('SkyCoordOption type_cast_value')
        try:
            #print('SkyCoordOption type_cast_value value {} type {}'.format(value, type(value)))
            #assert len(value)%2==0
            ra_dec_list = []
            for ra_dec_pair in value:
                ra, dec = ra_dec_pair
                #print('SkyCoordOption type_cast_value ra dec {} {}'.format(ra, dec))
                ra_dec_list.append(SkyCoordParamType_.convert((ra, dec), param=self, ctx=ctx))
                #print('SkyCoordOption type_cast_value ra_dec_list {}'.format(ra_dec_list))
            return ra_dec_list
        except Exception:
            raise click.BadParameter(value)




@click.command()
@click.argument('skycoord', nargs=2, cls=SkyCoordArgument, type=SkyCoordParamType_)
@click.option('--calib-level', type=str, default='3', help='Selecting calib level. "1" for uncal data, "2" for rate data, "3" for drizzled image data. Can be multiple like "1,2,3".')
@click.option('--extension', type=str, default='fits', help='Selecting extension. Usually just "fits". Can be "fits,json,jpg,ecsv" if you want to get those files too.')
@click.option('--preview', is_flag=True, default=False, help='Selecting preview files, e.g., i2d quicklook image in jpg format.')
@click.option('--auxiliary', is_flag=True, default=False, help='Selecting auxiliary files, e.g., guide star data.')
@click.option('--with-image', is_flag=True, default=True)
@click.option('--with-nircam-imaging/--no-nircam-imaging', is_flag=True, default=True)
@click.option('--with-miri-imaging/--no-miri-imaging', is_flag=True, default=True)
@click.option('--with-niriss-imaging/--no-niriss-imaging', is_flag=True, default=True)
@click.option('--with-spectroscopy', is_flag=True, default=True)
@click.option('--with-nirspec-msa/--no-nirspec-msa', is_flag=True, default=True)
@click.option('--with-nirspec-slit/--no-nirspec-slit', is_flag=True, default=True)
@click.option('--with-nircam-wfss/--no-nircam-wfss', is_flag=True, default=True)
@click.option('--with-niriss-wfss/--no-niriss-wfss', is_flag=True, default=True)
@click.option('--with-niriss-soss/--no-niriss-soss', is_flag=True, default=True)
@click.option('--with-ifu/--no-ifu', is_flag=True, default=True)
@click.option('--with-nirspec-ifu/--no-nirspec-ifu', is_flag=True, default=True)
@click.option('--with-miri-ifu/--no-miri-ifu', is_flag=True, default=True)
@click.option('--download/--no-download', is_flag=True, default=False)
@click.option('--download-dir', type=click.Path(exists=False), default='.')
@click.option('--overwrite', is_flag=True, default=False)
def main(
        skycoord, 
        calib_level, 
        extension, 
        preview, 
        auxiliary, 
        with_image, 
        with_nircam_imaging, 
        with_miri_imaging, 
        with_niriss_imaging, 
        with_spectroscopy, 
        with_nirspec_msa, 
        with_nirspec_slit, 
        with_nircam_wfss, 
        with_niriss_wfss, 
        with_niriss_soss, 
        with_ifu, 
        with_nirspec_ifu, 
        with_miri_ifu, 
        download, 
        download_dir, 
        overwrite, 
    ):

    # setup logger
    logger, log_file = setup_logger()
    logger.info('log file: {}'.format(log_file))
    
    # print user input skycoord
    logger.info('skycoord: {}'.format(skycoord))
    
    # print user input download
    logger.info('download: {}'.format(download))

    # determine instrument_names
    # available names: https://outerspace.stsci.edu/display/MASTDOCS/JWST+Instrument+Names
    instrument_names = []
    if with_image:
        if with_nircam_imaging:
            instrument_names.append('NIRCAM/IMAGE')
        if with_miri_imaging:
            instrument_names.append('MIRI/IMAGE')
        if with_niriss_imaging:
            instrument_names.append('NIRISS/IMAGE')
    if with_spectroscopy:
        if with_nirspec_msa:
            instrument_names.append('NIRSPEC/MSA')
        if with_nirspec_slit:
            instrument_names.append('NIRSPEC/SLIT')
        if with_nircam_wfss:
            instrument_names.append('NIRCAM/WFSS')
        if with_niriss_wfss:
            instrument_names.append('NIRISS/WFSS')
        if with_niriss_soss:
            instrument_names.append('NIRISS/SOSS')
    if with_ifu:
        if with_nirspec_ifu:
            instrument_names.append('NIRSPEC/IFU')
        if with_miri_ifu:
            instrument_names.append('MIRI/IFU')
    
    # skycoord_str
    skycoord_str = skycoord.to_string('hmsdms', precision=1).replace(' ', '_')

    # query the MAST databse
    logger.info('Running Observations.query_criteria')
    obs_lists = {}
    product_lists = {}
    for instrument_name in instrument_names:
        instrument_name_str = instrument_name.lower().replace('/','_').replace('*','')
        obs_list_csv_file = os.path.join(download_dir, f'obs_list_{skycoord_str}_{instrument_name_str}.csv')
        if os.path.exists(obs_list_csv_file) and not overwrite:
            logger.info(f'Reading from existing query "{obs_list_csv_file}". (Use --overwrite to force re-query.)')
            obs_list = Table.read(obs_list_csv_file)
            obs_list['obsid'] = [str(t) for t in obs_list['obsid']]
        elif os.path.exists(obs_list_csv_file+'.nodata') and not overwrite:
            logger.info(f'According to "{obs_list_csv_file}.nodata", no data found for this query. Skipping. (Use --overwrite to force re-query.)')
            obs_list = None
        else:
            logger.info('Running query with Observations.query_criteria(instrument_name={!r})'.format(instrument_name))
            radius = "0.0 deg"
            if instrument_name.find('WFSS')>=0:
                radius = "5.0 arcsec"
            obs_list = Observations.query_criteria(
                obs_collection = "JWST",
                instrument_name = instrument_name,
                coordinates = skycoord,
                radius = radius,
                #calib_level = [3], # cannot directly find level 1 products, must search level3 first then call get_product_list
            )
            logger.info('Returned {} rows'.format(len(obs_list)))
            if obs_list is not None and len(obs_list) == 0:
                obs_list = None
                logger.info(f'Writing to "{obs_list_csv_file}.nodata".')
                with open(obs_list_csv_file+'.nodata', 'w') as fp:
                    fp.write('')
            else:
                logger.info(f'Writing to "{obs_list_csv_file}".')
                if os.path.exists(obs_list_csv_file):
                    shutil.move(obs_list_csv_file, obs_list_csv_file+'.backup')
                obs_list.write(obs_list_csv_file)
        
        obs_lists[instrument_name] = obs_list

        if obs_list is None:
            continue
    
        # get level-1 (uncal) product list
        product_lists[instrument_name] = {}
        for obs in obs_list:
            if obs['provenance_name'] != 'APT': # APT means planned observations, there are no real data.
                product_list = Observations.get_product_list(obs) # str(obs['obsID'])
                product_lists[instrument_name][obs['obs_id']] = product_list








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
        
        # check obs_num constraint if provided
        is_obs_num_matched = True
        if obs_num is not None:
            is_obs_num_matched = False
            if isinstance(obs_num, str):
                obs_num_list = [obs_num]
            else:
                obs_num_list = obs_num
            for obs_num_str in obs_num_list:
                # example: jw01345-o052_t022_nircam_clear-f200w
                #          jw01345-o061_s00410_nirspec_f170lp-g235m
                #          jw01727-o099_t093_miri_f770w
                check_obs_id = 'jw{}-o{}_(t|s)[0-9]+_.*_.*'.format(proposal_id, obs_num_str)
                if re.match(check_obs_id, obs['obs_id']) is not None:
                    is_obs_num_matched = True
                    break
                else:
                    # example: jw01837004001_02101_00002_mirimage
                    check_obs_id = 'jw{}{}[0-9]+_[x0-9]+_[0-9]+_.*'.format(proposal_id, obs_num_str)
                    if re.match(check_obs_id, obs['obs_id']) is not None:
                        is_obs_num_matched = True
                        break
        if not is_obs_num_matched:
            continue
        
        # check image type constraint if provided
        is_image_type_matched = True
        if image_type is not None:
            is_image_type_matched = False
            if image_type == 'nircam':
                check_ending1 = 'nircam'
                check_ending2 = 'nrc(a|b|)(1|2|3|4|5|long|)'
            elif image_type == 'miri':
                check_ending1 = 'miri'
                check_ending2 = 'mirimage'
            else:
                raise Exception('Unexpected image_type {!r}'.format(image_type))
            # example: jw01345-o052_t022_nircam_clear-f200w
            check_image_type = 'jw{}-o[0-9]+_(t|s)[0-9]+_{}_.*'.format(proposal_id, check_ending1)
            print("obs['obs_id']", obs['obs_id'], 'check_image_type', check_image_type)
            if re.match(check_image_type, obs['obs_id']) is not None:
                is_image_type_matched = True
            else:
                # example: jw01837004001_02101_00002_mirimage
                check_image_type = 'jw{}[0-9]+_[x0-9]+_[0-9]+_{}_.*'.format(proposal_id, check_ending1)
                print("obs['obs_id']", obs['obs_id'], 'check_image_type', check_image_type)
                if re.match(check_image_type, obs['obs_id']) is not None:
                    is_image_type_matched = True
        print("obs['obs_id']", obs['obs_id'], 'is_image_type_matched', is_image_type_matched)
        if not is_image_type_matched:
            continue
        
        # 
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
                # de-select guide star data
                if auxiliary and (not guide_star_data):
                    products = products[np.array([not check_guide_star_data(product) for product in products])]
                # 
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



