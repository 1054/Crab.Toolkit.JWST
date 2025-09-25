#!/usr/bin/env python                                                                                                                                                                                 
#
# see -- https://mirage-data-simulator.readthedocs.io/en/latest/reference_files.html
# 
import os
from mirage.reference_files import downloader
download_path = os.path.expanduser('~/jwst_mirage_data')
downloader.download_reffiles(download_path, instrument='all', dark_type='linearized', skip_darks=False, single_dark=False, skip_cosmic_rays=False, skip_psfs=False)
#dark_type: linearized, raw, both
#instrument: all, 'NIRCam,NIRSpec,NIRISS,MIRI'

