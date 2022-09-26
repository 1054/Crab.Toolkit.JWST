#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
os.environ["CRDS_PATH"]       = '/n17data/dzliu/Data/jwst_crds_cache' # '/n23data1/hjmcc/jwst/mirage/crds_cache' #<DZLIU>#
#os.environ["CRDS_DATA"]       = '/n17data/dzliu/Data/jwst_crds_cache' # '/n23data1/hjmcc/jwst/mirage/crds_cache' #<DZLIU># # '/n23data1/mfranco/crds_cache'
os.environ["MIRAGE_DATA"]     = '/n23data1/hjmcc/jwst/mirage/mirage_data' # '/n23data1/hjmcc/jwst/mirage/mirage_data' #<DZLIU># # '/n23data1/mfranco/mirage_data'
os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu' #<DZLIU># # 'https://crds-serverless-mode.stsci.edu'

import click
import numpy as np
#import pysiaf
import photutils # used by mirage.seed_image.fits_seed_image
#from photutils.psf.matching import TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow # see "mirage/seed_image/fits_seed_image.py", this is a bug in mirage!
for key in 'TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow'.split(','): # workaround for the above bug
    setattr(photutils, key.strip(), getattr(photutils.psf.matching, key.strip()))
import jwst
import mirage # pip install --upgrade git+https://github.com/spacetelescope/mirage.git
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from functools import partial
from jwst.stpipe import Step
from jwst.outlier_detection import outlier_detection, OutlierDetectionStep
#outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, outlier_detection.OutlierDetection)
#outlier_detection.OutlierDetection.output_ext = 'fits'
#outlier_detection.OutlierDetection.make_output_path = OutlierDetectionStep.make_output_path
from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.catalogs import catalog_generator
from mirage.catalogs.catalog_generator import ExtendedCatalog
from mirage.dark import dark_prep
from mirage.imaging_simulator import ImgSim
from mirage.ramp_generator import obs_generator
from mirage.seed_image import catalog_seed_image, blot_image
from mirage.seed_image.fits_seed_image import ImgSeed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger(__name__)


print('jwst version:  '  ,jwst.__version__)
print('mirage version:', mirage.__version__)



class DzliuMirageSimulator():
    """DzliuMirageSimulator
    
    For COSMOS JWST MIRAGE simulation.
    
    We need input APT files in three formats: "{apt_name}.xml", "{apt_name}.pointing" and "{apt_name}.aptx".
    
    """
    def __init__(
            self, 
            apt_dir = 'apt_files', 
            apt_name = 'cosmosweb_revised_jun2022_onlyDEC2022', 
            apt_xml_file = None, 
            pointing_file = None, 
        ):
        self.logger = logging.getLogger('DzliuMirageSimulator')
        self.apt_dir = apt_dir
        self.apt_name = apt_name
        self.yam = None
        if apt_xml_file is None:
            apt_xml_file = os.path.join(apt_dir, apt_name+'.xml')
        if not os.path.isfile(apt_xml_file):
            raise Exception('Error! File not found: {!r}'.format(apt_xml_file))
        if pointing_file is None:
            pointing_file = os.path.join(apt_dir, apt_name+'.pointing')
        if not os.path.isfile(pointing_file):
            raise Exception('Error! File not found: {!r}'.format(pointing_file))
        if apt_xml_file is not None and pointing_file is not None:
            self.load_apt_xml_and_pointing_files(apt_xml_file, pointing_file)
    
    def load_apt_xml_and_pointing_files(
            self, 
            apt_xml_file, 
            pointing_file, 
        ):
        self.logger.info('Loading APT and pointing files: {!r}, {!r}'.format(apt_xml_file, pointing_file))
        self.apt_xml_file = apt_xml_file
        self.pointing_file = pointing_file
        #create_catalog.create_basic_exposure_list()
        #create_catalog.for_proposal()
        cosmic_rays = {'library': 'SUNMAX', 'scale': 1.0}
        self.yam = yaml_generator.SimInput(
            input_xml = self.apt_xml_file, 
            pointing_file = self.pointing_file,
            catalogs = cat_dict, 
            cosmic_rays = cosmic_rays,
            background = background, 
            roll_angle = pav3,
            dates = dates, 
            reffile_defaults = reffile_defaults,
            add_ghosts = ghosts, 
            convolve_ghosts_with_psf = convolve_ghosts,
            verbose = True, 
            output_dir = os.path.join(path,output_dir),
            simdata_output_dir = os.path.join(path,simulation_dir),
            datatype = datatype, 
            use_JWST_pipeline = True, 
        )
    
    




# Update Yaml File
templateParamFile = 'jw01727043001_01101_00001_nrca1.yaml' # do not overwrite this file
paramfile = 'jw01727043001_01101_00001_nrca1_go.yaml'
with open(templateParamFile, 'r') as ifp:
    y = yaml.safe_load(ifp)
    print("y['simSignals']['galaxyListFile'] = {}".format(y['simSignals']['galaxyListFile'])) # e.g., /automnt/n17data/dzliu/Work/20220919-Max-script/catalogs/galaxies_DEC_2022_input.cat
    print("y['simSignals']['pointsource'] = {}".format(y['simSignals']['pointsource'])) # e.g., /automnt/n17data/dzliu/Work/20220919-Max-script/catalogs/ptsrc_all_filters_region_BEST_sw.cat
    print("y['simSignals']['psfpath'] = {}".format(y['simSignals']['psfpath'])) # e.g., psfpath: /n23data1/hjmcc/jwst/mirage/mirage_data/nircam/gridded_psf_library
    with open(paramfile, 'w') as ofp:
        yaml.dump(y, ofp)


# MIRAGE ImgSim Tutorial
# https://github.com/spacetelescope/mirage/blob/master/examples/Simulated_data_from_mosaic_image.ipynb


# Define function to get pixel scale of a FITS image
def get_pixel_scale(fits_file):
    data, header = fits.getdata(fits_file, header=True)
    wcs = WCS(header, naxis=2)
    pixsc = np.sqrt(proj_plane_pixel_area(wcs)) * 3600.0
    return pixsc


# Define function to create a single-pixel PSF
def prepare_mosaic_psf_file(psf_file, pixel_scale=None, overwrite=False):
    data = np.zeros([201, 201])
    data[data.shape[0]//2, data.shape[1]//2] = 1.0
    hdu = fits.PrimaryHDU(data)
    hdu.header['TELESCOP'] = 'N/A'
    hdu.header['INSTRUME'] = 'N/A'
    if pixel_scale is not None:
        hdu.header['CD1_1'] = -float(pixel_scale)/3600.0 # see "mirage/psf/tools.py" get_psf_metadata
        hdu.header['CD2_2'] = float(pixel_scale)/3600.0
    hdu.writeto(psf_file, overwrite=overwrite)
    return


# Define function to get image center ra dec
def get_image_center(fits_file):
    data, header = fits.getdata(fits_file, header=True)
    wcs = WCS(header, naxis=2)
    (ra,), (dec,) = wcs.wcs_pix2world([(header['NAXIS1']+1.0)/2.0,], [(header['NAXIS2']+1.0)/2.0,], 1)
    return ra, dec


# Define function to hack image center
def hack_image_center(fits_file, out_file, ra, dec):
    data, header = fits.getdata(fits_file, header=True)
    wcs = WCS(header, naxis=2)
    radec = np.array([[ra, dec]])
    mosaic_xy_at_center_radec = wcs.wcs_world2pix(radec, 1)


# Resample Stamp Image
if True:

    # Read observation RA Dec PA
    with open(paramfile, 'r') as fp:
        y = yaml.safe_load(fp)
    ra = y['Telescope']['ra']
    dec = y['Telescope']['dec']
    pav3 = y['Telescope']['rotation']

    # Crop from the mosaic and resample for the desired detector/aperture
    mosaic_file = 'hlsp_candels_hst_acs_gs-tot-sect23_f814w_v1.0_drz.fits'
    magnitude = None
    #mosaic_fwhm = 0.045 # None -- does not work # 0.09 -- ERROR: FWHM of the mosaic image is larger than that of the JWST PSF. Unable to create a matching PSF kernel using photutils.
    #mosaic_fwhm_units = 'arcsec'
    # -- TypeError: get_HST_PSF() takes 1 positional argument but 2 were given -- "mirage/seed_image/fits_seed_image.py", line 692, in prepare_psf
    # TODO: change the fits file meta data, let TELESCOP to some others and provide a psf_file (OR use gaussian_psf=True)
    # create a single pixel psf_file
    psf_file = mosaic_file.rstrip('.fits') + '.psf.fits'
    pixel_scale = get_pixel_scale(mosaic_file)
    prepare_mosaic_psf_file(psf_file, pixel_scale=pixel_scale, overwrite=True)
    mosaic_fwhm = pixel_scale # 1
    mosaic_fwhm_units = 'arcsec' # 'pixels' -- there is a bug in "mirage/seed_image/fits_seed_image.py", line 364, in crop_and_blot -- UnboundLocalError: local variable 'mosaic_fwhm_arcsec' referenced before assignment
    instrument = 'nircam'
    filter_name = 'F200W'
    stamp_name = 'stamp_'+instrument+'_'+filter_name
    stamp_catalog_file = stamp_name+'.cat'
    seed = ImgSeed(
        paramfile = paramfile, 
        mosaic_file = mosaic_file, 
        mosaic_fwhm = mosaic_fwhm,
        mosaic_fwhm_units = mosaic_fwhm_units, 
        cropped_file = stamp_name+'_cropped.fits',
        blotted_file = stamp_name+'_cropped_blotted.fits', 
        outdir = './',
        psf_file = psf_file,
        gaussian_psf = False,
        #save_intermediates = True,
    )
    ra, dec = get_image_center(mosaic_file)
    seed.crop_center_ra = ra # see "mirage/seed_image/fits_seed_image.py" read_param_file -- In default it will use the ra dec in the param file. Here we test to crop the input image center.
    seed.crop_center_dec = dec
    seed.blot_center_ra = ra
    seed.blot_center_dec = dec
    #seed.outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, seed.outlier_detection.OutlierDetection)
    #seed.blot_image.outlier_detection.OutlierDetection.search_output_file = True
    #seed.blot_image.outlier_detection.OutlierDetection.output_file = 
    #blot_image.outlier_detection.OutlierDetection.output_ext = 'fits'
    #outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, outlier_detection.OutlierDetection)
    #outlier_detection.OutlierDetection.output_ext = 'fits'
    seed.crop_and_blot()
    # BUG AGAIN -- NameError: name 'psf_filename' is not defined -- "mirage/seed_image/fits_seed_image.py", line 708, in prepare_psf
    # FIXING -- edit "mirage/seed_image/fits_seed_image.py", change all "psf_filename" to "self.psf_file"
    # BUG AGAIN -- ValueError: cannot convert float NaN to integer -- "/mirage/seed_image/crop_mosaic.py", line 139, in extract
    # FIXING -- setting "seed.center_ra = ra # TODO: hack the center" and "...dec..." above
    # BUG AGAIN -- AttributeError: None object has no attribute 'search_output_file' -- "stpipe/step.py", line 1071, in _make_output_path
    # FIXING -- 
    #   edit "/home/dzliu/Software/CONDA/miniconda3/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py", 
    #       find "model_path = self.make_output_path(", go to next line, CHANGE 
    #           "basename=blot_root," -> "basepath=blot_root,"
    #       then add following new line:
    #           "ext='fits',"
    #       #find "self.make_output_path = pars.get", got to second next line, CHANGE
    #       #    "partial(Step._make_output_path, None)" -> "partial(Step._make_output_path, self)"
    #   edit "/home/dzliu/Software/CONDA/miniconda3/lib/python3.9/site-packages/mirage/seed_image/blot_image.py", 
    #   must pass a 'make_output_path' object to 'pars' when calling outlier_detection.OutlierDetection:
    #       find "outlier_detection.OutlierDetection()", two lines above, add new line:
    #           from functools import partial
    #           from jwst.stpipe import Step
    #           stepx = OutlierDetectionStep()
    #           pars['make_output_path'] = stepx.make_output_path # partial(Step._make_output_path, ) #<DZLIU>#
    #       
    # 

    # Read observation RA Dec PA
    with open(paramfile, 'r') as fp:
        y = yaml.safe_load(fp)
    ra = y['Telescope']['ra']
    dec = y['Telescope']['dec']
    pav3 = y['Telescope']['rotation']

    # Prepare stamp catalog
    stamp_catalog = ExtendedCatalog(
        filenames = [stamp_name+'_cropped_blotted.fits'],
        ra = [ra],
        dec = [dec],
        position_angle = [pav3],
    )
    stamp_catalog.add_magnitude_column(
        [str(magnitude)],
        #instrument = instrument,
        #filter_name = filter_name, # comment out to add a generic one, see "mirage/catalogs/catalog_generator.py"
    )
    #stamp_catalog.add_magnitude_column(
    #    [str(magnitude)],
    #    instrument = instrument,
    #    filter_name = 'F444W',
    #)
    stamp_catalog.save(
        stamp_catalog_file
    )


# Update Yaml file to use the stamp
if True:
    with open(paramfile, 'r') as ifp:
        y = yaml.safe_load(ifp)
    y['simSignals']['extended'] = stamp_catalog_file
    y['simSignals']['galaxyListFile'] = None # "mirage/catalogs/utils.py" will check lower case str of this being 'none' or not
    y['simSignals']['pointsource'] = None
    with open(paramfile, 'w') as ofp:
        yaml.dump(y, ofp)


# 
if True:
    m = ImgSim() # override_dark=dark
    m.paramfile = paramfile
    m.create()





