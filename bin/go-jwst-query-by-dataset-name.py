#!/usr/bin/env python
#
"""
Query astroquery.MAST by data set name, e.g., jwpppppooovvv_ggsaa_eeeee_detector_prodType.

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


@click.command()
@click.argument('dataset_names', nargs=-1, type=str)
@click.option('--download', is_flag=True, default=False)
@click.option('--download-dir', type=click.Path(exists=False), default='.')
def main(
        dataset_names, 
        download, 
        download_dir, 
    ):
    
    print('dataset_names: {}'.format(dataset_names))
    
    cached_obs_list = {}
    
    for dataset_name in dataset_names:
        
        jwst_dataset = parse_jwst_dataset_name(dataset_name)
        print('jwst_dataset: {}'.format(jwst_dataset))
        
        if jwst_dataset is None:
            raise Exception('Error! The input JWST dataset name {!r} seems incorrect!'.format(dataset_name))

        #calib_level = [3] # cannot directly find level -1 products, must search level3 first then call get_product_list
        extension = "fits"

        if jwst_dataset.proposal_id in cached_obs_list:
            obs_list = cached_obs_list[jwst_dataset.proposal_id]
        else:
            obs_list = Observations.query_criteria(
                obs_collection = "JWST",
                proposal_id = jwst_dataset.proposal_id,
                calib_level = [3], 
            )
        
        for obs in obs_list:
            print('obs_id: {}'.format(obs['obs_id']))
            check_obs_id = 'jw{}-o{}_t[0-9]+_{}_.*'.format(
                jwst_dataset.proposal_id, jwst_dataset.obs_num, get_instrument_from_detector_name(jwst_dataset.detector).lower())
            if re.match(check_obs_id, obs['obs_id']) is None:
                continue
            product_list = Observations.get_product_list(obs)
            if len(product_list) > 0:
                products = Observations.filter_products(
                    product_list,
                    productType = ["SCIENCE", "PREVIEW"],
                    calib_level = [1],
                    extension = extension,
                )
                if len(products) > 0: 
                    products = products[products['obs_id'] == dataset_name]
                # 
                if len(products) > 0: 
                    if download:
                        manifest = Observations.download_products(
                            products, 
                            download_dir = download_dir, 
                        )
                    else:
                        print('*** --- ***')
                        print(products)



# 
if __name__ == '__main__':
    main()



