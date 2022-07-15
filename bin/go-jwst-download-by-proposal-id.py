#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id.

"""

import os, sys, re, click
from astroquery.mast import Observations

#calib_level = []
#calib_level = [2,3] # [1,2,3]

@click.command()
@click.argument('proposal_id', type=int)
@click.option('--calib-level', type=list, default=[3]) # e.g., [1,2,3]
@click.option('--extension', type=str, default=None) # e.g., "fits"
@click.option('--download', type=bool, default=True) # e.g., "fits"
def main(proposal_id, calib_level, extension, download):
    calib_level_raw = calib_level
    calib_level = []
    for calib_level_item in calib_level_raw:
        if re.match('^[0-9]+$', str(calib_level_item)):
            calib_level.append(int(calib_level_item))
    print('------'*6)
    print('proposal_id: {}'.format(proposal_id))
    print('calib_level: {}'.format(calib_level))
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
                products = Observations.filter_products(product_list,
                                                        calib_level=calib_level,
                                                        extension=extension)
                # see -- https://astroquery.readthedocs.io/en/latest/mast/mast.html
                # productType=["SCIENCE", "PREVIEW"],
                # productSubGroupDescription='DRZ', 
                # 
                if len(products) > 0 and download:
                    manifest = Observations.download_products(products)
            elif download:
                manifest = Observations.download_products(product_list)



if __name__ == '__main__':
    main()


