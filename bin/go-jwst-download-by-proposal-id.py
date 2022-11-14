#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id.

"""

import os, sys, re, click
import numpy as np
from astroquery.mast import Observations

@click.command()
@click.argument('proposal_id', type=int)
@click.option('--download-dir', type=str, default=None)
@click.option('--calib-level', type=str, default='3', help='Calib level, e.g. "1" or "1,2,3".') # e.g., '1,2,3'
@click.option('--product-type', type=str, default='science,preview', help='Product type, science, preview and auxiliary.')
@click.option('--extension', type=str, default='fits', help='File extension, e.g., "fits" or "fits,jpg".') # e.g., 'fits,jpg'
@click.option('--dataset', type=str, multiple=True, default=[], help='Can select only one or multiple specific datasets, e.g., "jw01345002001_14201_00001_nrcb4".')
@click.option('--download/--no-download', is_flag=True, default=True)
@click.option('--guide-star/--no-guide-star', is_flag=True, default=False, help='Select guide star data or not. Default is no.')
def main(
        proposal_id, 
        download_dir, 
        calib_level, 
        product_type, 
        extension, 
        dataset, 
        download, 
        guide_star, 
    ):
    # 
    current_dir = os.getcwd()
    if download_dir is not None:
        if not os.path.isdir(download_dir):
            os.makedirs(download_dir)
        print('Current directory {!r}'.format(current_dir))
        print('Changing directory to {!r}'.format(download_dir))
        os.chdir(download_dir)
    # 
    calib_level_str = calib_level
    if re.match(r'\[(.*)\]', calib_level_str):
        calib_level_str = re.sub(r'\[(.*)\]', r'\1', calib_level_str)
    calib_level = []
    for calib_level_item in calib_level_str.split(','):
        calib_level.append(int(calib_level_item))
    # 
    if isinstance(dataset, str):
        dataset_list = [dataset]
    else:
        dataset_list = dataset
    if len(dataset_list) > 0:
        if 1 not in calib_level:
            calib_level.insert(0, 1)
    # 
    product_type_str = product_type
    if re.match(r'\[(.*)\]', product_type_str):
        product_type_str = re.sub(r'\[(.*)\]', r'\1', product_type_str)
    if product_type_str == '':
        product_type = ['SCIENCE', 'PREVIEW']
    else:
        product_type = []
        for product_type_item in product_type_str.split(','):
            product_type.append(product_type_item.upper())
    # 
    if guide_star:
        if 'AUXILIARY' not in product_type:
            product_type.append(product_type)
    # 
    extension_str = extension
    if re.match(r'\[(.*)\]', extension_str):
        extension_str = re.sub(r'\[(.*)\]', r'\1', extension_str)
    if extension_str == '':
        extension = None
    else:
        extension = []
        for extension_item in extension_str.split(','):
            extension.append(extension_item)
    # 
    print('------'*6)
    print('proposal_id: {}'.format(proposal_id))
    print('calib_level: {}'.format(calib_level))
    print('product_type: {}'.format(product_type))
    print('extension: {}'.format(extension))
    print('download: {}'.format(download))
    print('------'*6)
    proposal_id_str = '{:05d}'.format(proposal_id)
    obs_list = Observations.query_criteria(obs_collection="JWST",
                                           proposal_id=proposal_id_str)
    for obs in obs_list:
        product_list = Observations.get_product_list(obs)
        if len(product_list) > 0:
            if len(calib_level) > 0:
                products = Observations.filter_products(
                    product_list,
                    calib_level = calib_level,
                    productType = product_type,
                    extension = extension,
                )
                # see -- https://astroquery.readthedocs.io/en/latest/mast/mast.html
                # productType=["SCIENCE", "PREVIEW"],
                # productSubGroupDescription='DRZ', 
                # see -- https://masttest.stsci.edu/api/v0/_productsfields.html
                # 
                # 
                select = np.full(len(products), fill_value=True)
                for iproduct, product in enumerate(products):
                    if len(dataset_list) > 0:
                        for dataset_item in dataset_list:
                            if product['productFilename'].count(dataset_item) == 0:
                                select[iproduct] = False
                            else:
                                print("product['productFilename']: {} (select = {})".format(product['productFilename'], select[iproduct]))
                    #print("product['productFilename']: {} (select = {})".format(product['productFilename'], select[iproduct]))
                products = products[select]
                # 
                if len(products) > 0 and download:
                    manifest = Observations.download_products(products)
            elif download:
                manifest = Observations.download_products(product_list)



if __name__ == '__main__':
    main()


