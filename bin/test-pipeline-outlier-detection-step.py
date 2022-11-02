#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
os.environ["CRDS_PATH"] = '/n17data/dzliu/Data/jwst_crds_cache' # '/n23data1/hjmcc/jwst/mirage/crds_cache' #<DZLIU>#
os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu' # 'https://crds-serverless-mode.stsci.edu'
os.environ["CRDS_CONTEXT"] = 'jwst_1009.pmap' # 2022-10-26
import click
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area

import jwst

#from functools import partial
#from jwst.stpipe import Step
#from jwst.outlier_detection import outlier_detection, OutlierDetectionStep
#outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, outlier_detection.OutlierDetection)
#outlier_detection.OutlierDetection.output_ext = 'fits'
#outlier_detection.OutlierDetection.make_output_path = OutlierDetectionStep.make_output_path

from jwst.tweakreg import TweakRegStep
from jwst.skymatch import SkyMatchStep
from jwst.outlier_detection import OutlierDetectionStep
from jwst.resample import ResampleStep
from jwst.source_catalog import SourceCatalogStep
from jwst import datamodels
from jwst.datamodels import ImageModel
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('__name__')


print('jwst version: {}'.format(jwst.__version__))




visit_sca = 'combined'

visit_sca_imgfile_list = [
"../processing/jw01345004001_02201_00001_nrca1/calibrated2_cals/jw01345004001_02201_00001_nrca1_cal.fits",
"../processing/jw01345004001_02201_00001_nrca2/calibrated2_cals/jw01345004001_02201_00001_nrca2_cal.fits",
"../processing/jw01345004001_02201_00001_nrca3/calibrated2_cals/jw01345004001_02201_00001_nrca3_cal.fits",
"../processing/jw01345004001_02201_00001_nrca4/calibrated2_cals/jw01345004001_02201_00001_nrca4_cal.fits",
"../processing/jw01345004001_02201_00001_nrcb1/calibrated2_cals/jw01345004001_02201_00001_nrcb1_cal.fits",
"../processing/jw01345004001_02201_00001_nrcb2/calibrated2_cals/jw01345004001_02201_00001_nrcb2_cal.fits",
"../processing/jw01345004001_02201_00001_nrcb3/calibrated2_cals/jw01345004001_02201_00001_nrcb3_cal.fits",
"../processing/jw01345004001_02201_00001_nrcb4/calibrated2_cals/jw01345004001_02201_00001_nrcb4_cal.fits",
"../processing/jw01345004001_02201_00002_nrca1/calibrated2_cals/jw01345004001_02201_00002_nrca1_cal.fits",
"../processing/jw01345004001_02201_00002_nrca2/calibrated2_cals/jw01345004001_02201_00002_nrca2_cal.fits",
"../processing/jw01345004001_02201_00002_nrca3/calibrated2_cals/jw01345004001_02201_00002_nrca3_cal.fits",
"../processing/jw01345004001_02201_00002_nrca4/calibrated2_cals/jw01345004001_02201_00002_nrca4_cal.fits",
"../processing/jw01345004001_02201_00002_nrcb1/calibrated2_cals/jw01345004001_02201_00002_nrcb1_cal.fits",
"../processing/jw01345004001_02201_00002_nrcb2/calibrated2_cals/jw01345004001_02201_00002_nrcb2_cal.fits",
"../processing/jw01345004001_02201_00002_nrcb3/calibrated2_cals/jw01345004001_02201_00002_nrcb3_cal.fits",
"../processing/jw01345004001_02201_00002_nrcb4/calibrated2_cals/jw01345004001_02201_00002_nrcb4_cal.fits",
"../processing/jw01345004001_02201_00003_nrca1/calibrated2_cals/jw01345004001_02201_00003_nrca1_cal.fits",
"../processing/jw01345004001_02201_00003_nrca2/calibrated2_cals/jw01345004001_02201_00003_nrca2_cal.fits",
"../processing/jw01345004001_02201_00003_nrca3/calibrated2_cals/jw01345004001_02201_00003_nrca3_cal.fits",
"../processing/jw01345004001_02201_00003_nrca4/calibrated2_cals/jw01345004001_02201_00003_nrca4_cal.fits",
"../processing/jw01345004001_02201_00003_nrcb1/calibrated2_cals/jw01345004001_02201_00003_nrcb1_cal.fits",
"../processing/jw01345004001_02201_00003_nrcb2/calibrated2_cals/jw01345004001_02201_00003_nrcb2_cal.fits",
"../processing/jw01345004001_02201_00003_nrcb3/calibrated2_cals/jw01345004001_02201_00003_nrcb3_cal.fits",
"../processing/jw01345004001_02201_00003_nrcb4/calibrated2_cals/jw01345004001_02201_00003_nrcb4_cal.fits",
"../processing/jw01345004001_04201_00001_nrca1/calibrated2_cals/jw01345004001_04201_00001_nrca1_cal.fits",
"../processing/jw01345004001_04201_00001_nrca2/calibrated2_cals/jw01345004001_04201_00001_nrca2_cal.fits",
"../processing/jw01345004001_04201_00001_nrca3/calibrated2_cals/jw01345004001_04201_00001_nrca3_cal.fits",
"../processing/jw01345004001_04201_00001_nrca4/calibrated2_cals/jw01345004001_04201_00001_nrca4_cal.fits",
"../processing/jw01345004001_04201_00001_nrcb1/calibrated2_cals/jw01345004001_04201_00001_nrcb1_cal.fits",
"../processing/jw01345004001_04201_00001_nrcb2/calibrated2_cals/jw01345004001_04201_00001_nrcb2_cal.fits",
"../processing/jw01345004001_04201_00001_nrcb3/calibrated2_cals/jw01345004001_04201_00001_nrcb3_cal.fits",
"../processing/jw01345004001_04201_00001_nrcb4/calibrated2_cals/jw01345004001_04201_00001_nrcb4_cal.fits",
"../processing/jw01345004001_04201_00002_nrca1/calibrated2_cals/jw01345004001_04201_00002_nrca1_cal.fits",
"../processing/jw01345004001_04201_00002_nrca2/calibrated2_cals/jw01345004001_04201_00002_nrca2_cal.fits",
"../processing/jw01345004001_04201_00002_nrca3/calibrated2_cals/jw01345004001_04201_00002_nrca3_cal.fits",
"../processing/jw01345004001_04201_00002_nrca4/calibrated2_cals/jw01345004001_04201_00002_nrca4_cal.fits",
"../processing/jw01345004001_04201_00002_nrcb1/calibrated2_cals/jw01345004001_04201_00002_nrcb1_cal.fits",
"../processing/jw01345004001_04201_00002_nrcb2/calibrated2_cals/jw01345004001_04201_00002_nrcb2_cal.fits",
"../processing/jw01345004001_04201_00002_nrcb3/calibrated2_cals/jw01345004001_04201_00002_nrcb3_cal.fits",
"../processing/jw01345004001_04201_00002_nrcb4/calibrated2_cals/jw01345004001_04201_00002_nrcb4_cal.fits",
"../processing/jw01345004001_04201_00003_nrca1/calibrated2_cals/jw01345004001_04201_00003_nrca1_cal.fits",
"../processing/jw01345004001_04201_00003_nrca2/calibrated2_cals/jw01345004001_04201_00003_nrca2_cal.fits",
"../processing/jw01345004001_04201_00003_nrca3/calibrated2_cals/jw01345004001_04201_00003_nrca3_cal.fits",
"../processing/jw01345004001_04201_00003_nrca4/calibrated2_cals/jw01345004001_04201_00003_nrca4_cal.fits",
"../processing/jw01345004001_04201_00003_nrcb1/calibrated2_cals/jw01345004001_04201_00003_nrcb1_cal.fits",
"../processing/jw01345004001_04201_00003_nrcb2/calibrated2_cals/jw01345004001_04201_00003_nrcb2_cal.fits",
"../processing/jw01345004001_04201_00003_nrcb3/calibrated2_cals/jw01345004001_04201_00003_nrcb3_cal.fits",
"../processing/jw01345004001_04201_00003_nrcb4/calibrated2_cals/jw01345004001_04201_00003_nrcb4_cal.fits",
]



def test():
    
    # following Anton's script
    
    asn = asn_from_list.asn_from_list(visit_sca_imgfile_list, rule=DMS_Level3_Base, product_name=visit_sca)
    asn_file = 'tweakreg_asn.json'
    
    with open(asn_file, 'w') as outfile:
        name, serialized = asn.dump(format='json')
        outfile.write(serialized)
    
    tweakreg = TweakRegStep()
    tweakreg.save_catalogs = True
    tweakreg.save_results = True
    tweakreg.run(asn_file)
    
    
    
    # 
    
    imgfile_list = sorted(glob.glob('*_tweakreg.fits'))
    asn_file = 'outlier_detection_asn.json'
    
    with open(asn_file, 'w') as outfile:
        name, serialized = asn.dump(format='json')
        outfile.write(serialized)
    
    outlier_detection = OutlierDetectionStep()
    outlier_detection.resample_data = True
    outlier_detection.in_memory = False
    outlier_detection.output_dir = './'
    outlier_detection.run(asn_file)




if __name__ == '__main__':
    
    test()


