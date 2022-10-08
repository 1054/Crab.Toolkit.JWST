#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id, get science target uncal fits files. 

Rate, cal and i2d files are excluded unless specified.

The default download dir is 
./mastDownload/JWST/jwpppppooovvv_ggsaa_eeeee_detector/jwpppppooovvv_ggsaa_eeeee_detector_prodType.fits

"""

import re
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



@click.command()
@click.argument('program', type=str)
@click.option('--obs-num', '--obsnum', '--obs', type=str, default=None)
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
    
    # print user input program
    print('program: {}'.format(program))
    
    # check if program is a proposal ID or a known survey name
    if re.match(r'^[0-9]+$', program):
        proposal_id = '{:05d}'.format(int(program))
    else:
        proposal_id = get_proposal_id_from_survey_name(program)
    
    if proposal_id is None:
        print('Error! Cannot understand the input program ID. Example: 01345')
        sys.exit(255)
    
    # print proposal id to query
    print('proposal_id: {}'.format(proposal_id))
    
    # check if user has specified an obs num
    if obs_num is not None:
        if obs_num.startswith('obs'):
            obs_num = re.sub(r'^obs', r'', obs_num)
        
        if re.match(r'^[0-9]+$', obs_num):
            obs_num = '{:03d}'.format(int(obs_num))
        else:
            print('Error! Cannot understand the input obs_num. Example: 001')
            sys.exit(255)
    
        print('obs_num: {}'.format(obs_num))
    
    # query the MAST databse
    obs_list = Observations.query_criteria(
        obs_collection = "JWST",
        proposal_id = proposal_id,
        calib_level = [3], # cannot directly find level 1 products, must search level3 first then call get_product_list
    )
    
    # parse user input calib_level, extension, etc. for later filtering the product list
    if calib_level.find(',')>=0:
        calib_level = calib_level.replace(' ','').split(',')
    else:
        calib_level = [calib_level]
    
    if extension.find(',')>=0:
        extension = extension.replace(' ','').split(',')
    
    productType = ["SCIENCE"]
    if preview:
        productType.append("PREVIEW")
        if 'jpg' not in extension:
            extension.append('jpg')
    if auxiliary:
        productType.append("AUXILIARY")
    
    for iobs, obs in enumerate(obs_list):
        #print('obs_id: {}'.format(obs['obs_id']))
        if obs_num is not None:
            check_obs_id = 'jw{}-o{}_t[0-9]+_.*_.*'.format(
                proposal_id, obs_num)
            if re.match(check_obs_id, obs['obs_id']) is None:
                continue
        product_list = Observations.get_product_list(obs)
        if len(product_list) > 0:
            products = Observations.filter_products(
                product_list,
                calib_level = calib_level,
                productType = productType,
                extension = extension,
            )
            if len(products) > 0: 
                print('*** --- ({}/{}) --- ***'.format(iobs+1, len(obs_list)))
                print('obs_id: {}'.format(obs['obs_id']))
                #print(products.colnames)
                for product in products:
                    print('        {} {}'.format(product['obs_id'], product['productFilename']))
                    #print(product)
            # 
            if len(products) > 0: 
                if download:
                    manifest = Observations.download_products(
                        products, 
                        download_dir = download_dir, 
                    )



# 
if __name__ == '__main__':
    main()



