#!/usr/bin/env python
#
"""
Query astroquery.MAST by program id.

"""

from astroquery.mast import Observations

#calib_level = []
calib_level = [1,2,3] # [3] # [-1,2,3]

download = True

extension = None # "fits" # None if calib_level 1 2

obs_list = Observations.query_criteria(obs_collection="JWST",
                                       proposal_id="01837")

for obs in obs_list:
    product_list = Observations.get_product_list(obs)
    if len(product_list) > 0:
        if len(calib_level) > 0:
            products = Observations.filter_products(product_list,
                                                    calib_level=calib_level,
                                                    extension=extension)
            if len(products) > 0 and download:
                manifest = Observations.download_products(products)
        elif download:
            manifest = Observations.download_products(product_list)



