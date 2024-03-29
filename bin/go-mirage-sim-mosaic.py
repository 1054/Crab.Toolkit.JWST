#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = os.path.expanduser('~/jwst_crds_cache')
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = os.path.expanduser('~/jwst_mirage_data')
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

import click
import astropy.units as u
import astropy.constants as const
import numpy as np
#import pysiaf
import photutils # used by mirage.seed_image.fits_seed_image
#from photutils.psf.matching import TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow # see "mirage/seed_image/fits_seed_image.py", this is a bug in mirage!
for key in 'TopHatWindow, TukeyWindow, CosineBellWindow, SplitCosineBellWindow, HanningWindow'.split(','): # workaround for the above bug
    setattr(photutils, key.strip(), getattr(photutils.psf.matching, key.strip()))
import jwst
import mirage # pip install --upgrade git+https://github.com/spacetelescope/mirage.git
from astropy.coordinates import SkyCoord, FK5
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict
from functools import partial
from jwst.stpipe import Step
from jwst.outlier_detection import outlier_detection, OutlierDetectionStep
#outlier_detection.OutlierDetection.make_output_path = partial(Step._make_output_path, outlier_detection.OutlierDetection)
#outlier_detection.OutlierDetection.output_ext = 'fits'
#outlier_detection.OutlierDetection.make_output_path = OutlierDetectionStep.make_output_path
from mirage import imaging_simulator
from mirage.catalogs import create_catalog
from mirage.catalogs import catalog_generator
from mirage.dark import dark_prep
from mirage.imaging_simulator import ImgSim
from mirage.ramp_generator import obs_generator
from mirage.seed_image import catalog_seed_image, blot_image
from mirage.seed_image.fits_seed_image import ImgSeed
from mirage.utils.utils import ensure_dir_exists
from mirage.yaml import yaml_generator

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirage-sim-mosaic')


DEFAULT_APT_XML_FILE = 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.xml'
DEFAULT_POINTING_FILE = 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.pointing'
DEFAULT_MOSAIC_FILE = None # 'input_mosaic_images/dust_opa_u_lensed_H_EUC.fits' # 'hlsp_candels_hst_acs_gs-tot-sect23_f814w_v1.0_drz.fits'
DEFAULT_MOSAIC_CENTER = None # The coordinate in the mosaic image WCS # ('10:00:28.6', '02:12:21.0')
DEFAULT_STAR_CATALOG = None # 'input_catalogs/ptsrc_pointings_BEST_sw_tot.cat'
DEFAULT_COSMIC_RAYS = {'library': 'SUNMAX', 'scale': 1.0}
DEFAULT_INSTRUMENT = 'NIRCam'
DEFAULT_FILTER = 'F277W'
DEFAULT_DATES = '2023-01-01'
DEFAULT_BACKGROUND = 'medium'
DEFAULT_PA_V3 = 293.09730273
DEFAULT_YAML_OUTPUT_DIR = 'yaml_files'
DEFAULT_SIM_OUTPUT_DIR = 'sim_data'
DEFAULT_OBSERVATION_LIST_FILE = 'observation_list.yaml'
DEFAULT_NEW_PSF_LIBRARY = None # '~/Work/JWST-MIRAGE-Simulation/20221008_make_mirage_psf_library/mirage_data/nircam/gridded_psf_library'
DEFAULT_NEW_PSF_WING_THRESHOLD = None # '~/Work/JWST-MIRAGE-Simulation/20221008_make_mirage_psf_library/custom_nircam_psf_wing_rate_thresholds.txt'



# # Update Yaml File
# templateParamFile = 'jw01727043001_01101_00001_nrca1.yaml' # do not overwrite this file
# paramfile = 'jw01727043001_01101_00001_nrca1_go.yaml'
# with open(templateParamFile, 'r') as ifp:
#     y = yaml.safe_load(ifp)
#     print("y['simSignals']['galaxyListFile'] = {}".format(y['simSignals']['galaxyListFile'])) # e.g., ~/Work/20220919*/catalogs/galaxies_DEC_2022_input.cat
#     print("y['simSignals']['pointsource'] = {}".format(y['simSignals']['pointsource'])) # e.g., ~/Work/20220919*/catalogs/ptsrc_all_filters_region_BEST_sw.cat
#     print("y['simSignals']['psfpath'] = {}".format(y['simSignals']['psfpath'])) # e.g., psfpath: ~/jwst_mirage_data/nircam/gridded_psf_library
#     with open(paramfile, 'w') as ofp:
#         yaml.dump(y, ofp)


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


# Define function to hack image center <TODO>
def hack_image_center(fits_file, out_file, ra, dec):
    data, header = fits.getdata(fits_file, header=True)
    wcs = WCS(header, naxis=2)
    radec = np.array([[ra, dec]])
    mosaic_xy_at_center_radec = wcs.wcs_world2pix(radec, 1)


# Define function to resample the input mosaic image to detector WCS and PSF (including convolution)
def resample_mosaic_image(
        yaml_file, 
        mosaic_file,
        mosaic_psf = None, # can be a float number or psf_file
        instrument = 'nircam',
        filter_name = 'F200W',
        output_name = None, 
        output_dir = './',
        recenter = False, 
        center_ra = None, 
        center_dec = None, 
        overwrite = False, 
        verbose = True, 
    ):
    """
    Inputs:
        User given mosaic image and psf. 
        If mosaic_psf is None, we will construct a single-pixel PSF.
    
    Outputs:
        FITS images for Mirage simulation, 
            output_dir+'/'+output_prefix+instrument+'_'+filter_name+'_cropped.fits'
            output_dir+'/'+output_prefix+instrument+'_'+filter_name+'_cropped_blotted.fits'
        An ASCII file as an extended source catalog for Mirage simulation,
            output_dir+'/'+output_prefix+instrument+'_'+filter_name+'.cat'
    
    Options:
        `recenter` means that we will make the cutout centered at the given 
        `center_ra`, `center_dec` coordinate, which are in the input `mosaic_file` WCS. 
    
    """
    # Crop from the mosaic and resample for the desired detector/aperture
    #mosaic_fwhm = 0.045 # None -- does not work # 0.09 -- ERROR: FWHM of the mosaic image is larger than that of the JWST PSF. Unable to create a matching PSF kernel using photutils.
    #mosaic_fwhm_units = 'arcsec'
    # -- TypeError: get_HST_PSF() takes 1 positional argument but 2 were given -- "mirage/seed_image/fits_seed_image.py", line 692, in prepare_psf
    # TODO: change the fits file meta data, let TELESCOP to some others and provide a psf_file (OR use gaussian_psf=True)
    
    # Check output file
    if output_name is None:
        output_name = 'resampled_'+instrument+'_'+filter_name
    elif output_name.endswith('.fits'):
        output_name = re.sub(r'\.fits$', r'', output_name)
    output_resampled_file = os.path.join(output_dir, output_name+'.fits')
    if os.path.isfile(output_resampled_file):
        if overwrite:
            shutil.move(output_resampled_file, output_resampled_file+'.backup')
        else:
            if verbose:
                logger.info('Found existing resampled image {!r} and overwrite is set to False.'.format(output_resampled_file))
            return output_resampled_file
    # 
    # Check output dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # 
    # Get pixel scale
    pixel_scale = get_pixel_scale(mosaic_file)
    # 
    # Prepare psf_file
    if mosaic_psf is None:
        # create a single pixel psf_file
        psf_file = re.sub(r'\.fits$', r'', mosaic_file) + '.psf.fits'
        prepare_mosaic_psf_file(psf_file, pixel_scale=pixel_scale, overwrite=True)
        mosaic_fwhm = pixel_scale / 10. # 1 #<TODO><20230502># try this
        mosaic_fwhm_units = 'arcsec' # 'pixels' -- there is a bug in "mirage/seed_image/fits_seed_image.py", line 364, in crop_and_blot -- UnboundLocalError: local variable 'mosaic_fwhm_arcsec' referenced before assignment
    else:
        #if re.match(r'^[0-9eE.+-]+$', mosaic_psf):
        try:
            float(mosaic_psf)
            psf_file = None # 
            mosaic_fwhm = float(mosaic_psf)
            mosaic_fwhm_units = 'arcsec'
        except ValueError:
            psf_file = mosaic_psf
            mosaic_fwhm = None # If None, a Gaussian2D model will be fit to the PSF to estimate the FWHM -- "mirage/seed_image/fits_seed_image.py"
            mosaic_fwhm_units = 'arcsec'
    # 
    # Create ImgSeed
    seed = ImgSeed(
        paramfile = yaml_file, 
        mosaic_file = mosaic_file, 
        mosaic_fwhm = mosaic_fwhm,
        mosaic_fwhm_units = mosaic_fwhm_units, 
        cropped_file = os.path.join(output_dir, output_name+'_intermediate_cropped.fits'), # bug/feature: we need to add the dir path here.
        blotted_file = output_name+'.fits', 
        outdir = output_dir,
        psf_file = psf_file,
        gaussian_psf = False,
        #save_intermediates = True,
    )
    # 
    # Check if need to recenter
    if recenter:
        if center_ra is None or center_dec is None:
            ra, dec = get_image_center(mosaic_file)
            if center_ra is None:
                center_ra = ra
            if center_dec is None:
                center_dec = dec
        if isinstance(center_ra, str) or isinstance(center_dec, str):
            if re.match(r'^[0-9.+-]+$', center_ra) and re.match(r'^[0-9.+-]+$', center_dec):
                center_ra, center_dec = float(center_ra), float(center_dec)
            else:
                center_scoord = SkyCoord(str(center_ra), str(center_dec), unit=(u.hour, u.deg), frame=FK5) #TODO frame FK5
                center_ra, center_dec = center_scoord.ra.deg, center_scoord.dec.deg
        seed.crop_center_ra = center_ra # see "mirage/seed_image/fits_seed_image.py" read_param_file -- In default it will use the ra dec in the param file. Here we test to crop the input image center.
        seed.crop_center_dec = center_dec
        seed.blot_center_ra = center_ra
        seed.blot_center_dec = center_dec
    # 
    # Run the crop_and_blot
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
    #   edit "~/Software/CONDA/miniconda3/lib/python3.9/site-packages/jwst/outlier_detection/outlier_detection.py", 
    #       find "model_path = self.make_output_path(", go to next line, CHANGE 
    #           "basename=blot_root," -> "basepath=blot_root,"
    #       then add following new line:
    #           "ext='fits',"
    #       #find "self.make_output_path = pars.get", got to second next line, CHANGE
    #       #    "partial(Step._make_output_path, None)" -> "partial(Step._make_output_path, self)"
    #   edit "~/Software/CONDA/miniconda3/lib/python3.9/site-packages/mirage/seed_image/blot_image.py", 
    #   must pass a 'make_output_path' object to 'pars' when calling outlier_detection.OutlierDetection:
    #       find "outlier_detection.OutlierDetection()", two lines above, add new line:
    #           from functools import partial
    #           from jwst.stpipe import Step
    #           stepx = OutlierDetectionStep()
    #           pars['make_output_path'] = stepx.make_output_path # partial(Step._make_output_path, )
    #       
    # BUG AGAIN -- AttributeError: None object has no attribute 'search_output_file' -- "stpipe/step.py", line 1071, in _make_output_path
    # FIXING -- 
    #   edit "~/Software/CONDA/miniconda3/lib/python3.9/site-packages/mirage/seed_image/save_seed.py", 
    #   comment out the following line:
    #       #<DZLIU># kw['PIXARMAP'] = parameters['Reffiles']['pixelAreaMap'] #<DZLIU># pixelAreaMap is not used
    # 
    if verbose:
        logger.info('Output to {!r}'.format(output_resampled_file))
    # 
    return output_resampled_file


# Define function to get RA Dec PAV3 from yamlfile
def get_ra_dec_pav3_from_yamlfile(yamlfile):
    # Read observation RA Dec PA
    with open(yamlfile, 'r') as fp:
        y = yaml.safe_load(fp)
    ra = y['Telescope']['ra']
    dec = y['Telescope']['dec']
    pav3 = y['Telescope']['rotation']
    return ra, dec, pav3


# Define function to get catalog row number
def get_catalog_row_number(catalog_file):
    tb = ascii.read(catalog_file)
    return len(tb)


# Define function to check yamlfile
def check_yamlfile_matches_the_filter(
        yamlfile, 
        filter_name = None, 
    ):
    check_okay = False
    with open(yamlfile, 'r') as fp:
        y = yaml.safe_load(fp)
    if 'Readout' in y:
        if 'filter' in y['Readout']:
            if filter_name is None:
                check_okay = True
            elif y['Readout']['filter'].lower() == filter_name.lower():
                check_okay = True
    return check_okay


# Define function to check yamlfile
def check_yamlfile_has_extended_catalog(
        yamlfile, 
        extended_catalog_file = None, 
    ):
    check_okay = False
    with open(yamlfile, 'r') as fp:
        y = yaml.safe_load(fp)
    if 'simSignals' in y:
        if 'extended' in y['simSignals']:
            if extended_catalog_file is None:
                check_okay = True
            elif y['simSignals']['extended'] == extended_catalog_file:
                check_okay = True
    return check_okay


# Define function to update yamlfile
def update_yamlfile_with_extended_catalog(
        yamlfile, 
        extended_catalog_file, 
        extended_scale = 1.0, 
        output_yamlfile = None,  # None for in-place update
        verbose = True, 
    ):
    with open(yamlfile, 'r') as fp:
        y = yaml.safe_load(fp)
    y['simSignals']['extended'] = extended_catalog_file
    y['simSignals']['extendedscale'] = extended_scale
    y['simSignals']['PSFConvolveExtended'] = False # we have already convolved the PSF (blotted) in resample_mosaic_image()!
    #y['simSignals']['galaxyListFile'] = None # "mirage/catalogs/utils.py" will check lower case str of this being 'none' or not
    #y['simSignals']['pointsource'] = None
    if output_yamlfile is None:
        output_yamlfile = yamlfile
    if os.path.isfile(output_yamlfile):
        shutil.copy2(output_yamlfile, output_yamlfile+'.backup')
    with open(output_yamlfile, 'w') as fp:
        yaml.dump(y, fp)
    if verbose:
        logger.info('Updated {!r} with extended catalog {!r}'.format(output_yamlfile, extended_catalog_file))
    return


# Define function to parse SkyCoord
class SkyCoordParamType(click.ParamType):
    name = "SkyCoord"

    def convert(self, value, param, ctx):
        scoord = None
        if isinstance(value, SkyCoord):
            scoord = value
        
        elif isinstance(value, str):
            if value.find(':') >= 0 or value.find('h') >= 0:
                try:
                    scoord = SkyCoord(value, unit=(u.hourangle, u.deg), frame=FK5)
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
            else:
                try:
                    scoord = SkyCoord(value, unit=(u.deg, u.deg), frame=FK5)
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
        
        elif isinstance(value, (tuple, list)) and len(value) == 2:
            if str(value[0]).find(':') >= 0 or str(value[0]).find('h') >= 0:
                try:
                    scoord = SkyCoord(value[0], value[1], unit=(u.hourangle, u.deg), frame=FK5)
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
            else:
                try:
                    scoord = SkyCoord(value[0], value[1], unit=(u.deg, u.deg), frame=FK5)
                except ValueError:
                    self.fail(f"{value!r} is not a valid SkyCoord", param, ctx)
        
        return scoord

SkyCoordParamType_ = SkyCoordParamType()


class SkyCoordArgument(click.Argument):
    
    def type_cast_value(self, ctx, value) -> list:
        try:
            ra, dec = value
            ra_dec_skycoord = SkyCoordParamType_.convert((ra, dec), param=self, ctx=ctx)
            return ra_dec_skycoord
        except Exception:
            raise click.BadParameter(value)


class SkyCoordOption(click.Option):
    
    def type_cast_value(self, ctx, value) -> list:
        if value is None:
            return None
        try:
            ra, dec = value
            ra_dec_skycoord = SkyCoordParamType_.convert((ra, dec), param=self, ctx=ctx)
            return ra_dec_skycoord
        except Exception:
            raise click.BadParameter(value)



####################
### MAIN PROGRAM ###
####################

@click.command()
@click.option('--xml-file', type=click.Path(exists=True), default=DEFAULT_APT_XML_FILE)
@click.option('--pointing-file', type=click.Path(exists=True), default=DEFAULT_POINTING_FILE)
@click.option('--mosaic-file', type=click.Path(exists=True), default=DEFAULT_MOSAIC_FILE)
@click.option('--mosaic-center', nargs=2, cls=SkyCoordOption, type=SkyCoordParamType_, default=DEFAULT_MOSAIC_CENTER, help='Center coordinate in the input mosaic image.')
@click.option('--mosaic-psf', type=str, default=None, help='If the input mosaic image has a PSF, please provide a PSF FWHM in arcsec or a PSF file. In the latter case, a Gaussian2D model will be fit to the PSF to estimate the FWHM, see "mirage/seed_image/fits_seed_image.py".')
@click.option('--star-catalog', type=click.Path(exists=True), default=DEFAULT_STAR_CATALOG)
@click.option('--galaxy-catalog', type=click.Path(exists=True), default=None)
@click.option('--instrument', type=str, default=DEFAULT_INSTRUMENT)
@click.option('--filter', 'filter_name', type=str, default=DEFAULT_FILTER)
@click.option('--dates', type=str, default=DEFAULT_DATES)
@click.option('--background', type=str, default=DEFAULT_BACKGROUND)
@click.option('--pa-v3', type=float, default=DEFAULT_PA_V3)
@click.option('--yaml-output-dir', type=click.Path(exists=False), default=DEFAULT_YAML_OUTPUT_DIR)
@click.option('--sim-output-dir', type=click.Path(exists=False), default=DEFAULT_SIM_OUTPUT_DIR)
@click.option('--observation-list-file', type=click.Path(exists=False), default=DEFAULT_OBSERVATION_LIST_FILE)
@click.option('--only-dataset', type=str, multiple=True, default=[], help='Can specify which dataset(s) to process.')
@click.option('--new-psf-library', type=click.Path(exists=True), default=DEFAULT_NEW_PSF_LIBRARY, help='New psf library directory.')
@click.option('--new-psf-wing-threshold', type=click.Path(exists=True), default=DEFAULT_NEW_PSF_LIBRARY, help='New psf library directory.')
@click.option('--match-zeropoints/--no-match-zeropoints', is_flag=True, default=True, help='Create a temporary flux_cal file with zeropoints matched to the CRDS photom reference file.')
@click.option('--overwrite-siminput', is_flag=True, default=False)
@click.option('--overwrite-simprep', is_flag=True, default=False)
@click.option('--overwrite-simdata', is_flag=True, default=False)
@click.option('--skip-simdata', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        xml_file, 
        pointing_file, 
        mosaic_file, 
        mosaic_center,
        mosaic_psf,
        star_catalog, 
        galaxy_catalog, 
        instrument, 
        filter_name, 
        dates, 
        background, 
        pa_v3,
        yaml_output_dir,
        sim_output_dir,
        observation_list_file,
        only_dataset,
        new_psf_library,
        new_psf_wing_threshold, 
        match_zeropoints, 
        overwrite_siminput, 
        overwrite_simprep, 
        overwrite_simdata, 
        skip_simdata, 
        verbose, 
    ):
    
    if verbose:
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    if str(mosaic_file).lower() in ['none', '']:
        mosaic_file = None
    
    if str(star_catalog).lower() in ['none', '']:
        star_catalog = None
    
    if str(galaxy_catalog).lower() in ['none', '']:
        galaxy_catalog = None
    
    if mosaic_file is None and star_catalog is None and galaxy_catalog is None:
        errmsg = 'No input source to simulate! (mosaic_file is None and star_catalog is None and galaxy_catalog is None)'
        logger.error(errmsg)
        raise Exception(errmsg)
        #sys.exit(255)
    
    if mosaic_file is not None:
        mosaic_metadata = mirage.psf.tools.get_psf_metadata(mosaic_file) # make sure this runs
    
    if mosaic_psf == '':
        mosaic_psf = None
    
    # Check input catalogs
    catalogs = {}
    if star_catalog is not None:
        catalogs['point_source'] = star_catalog
    if galaxy_catalog is not None:
        catalogs['galaxy'] = galaxy_catalog
        
    if len(catalogs) == 0:
        catalogs = None
    
    if observation_list_file.find(os.sep) < 0:
        observation_list_file = os.path.join(yaml_output_dir, observation_list_file)
    
    # Check observation list file, the output of yaml_generator.
    # It is `Table(yam.info)` after calling yam.create_inputs(), see "mirage/yaml/yaml_generator.py"
    observation_list_final_file = os.path.join(
        yaml_output_dir,
        'Observation_table_for_' + os.path.basename(xml_file) +
        '_with_yaml_parameters.csv' # following the final_file variable construction in "mirage/yaml/yaml_generator.py"
    )
    
    # Check observation list file
    if os.path.isfile(observation_list_final_file): 
        if overwrite_siminput:
            shutil.move(observation_list_final_file, observation_list_final_file+'.backup')
        else:
            if verbose:
                logger.info('Found existing observation list file {!r} and overwrite is set to False.'.format(observation_list_final_file))
    
    # Prepare YAML file using the input APT xml file and pointing file
    # -- https://mirage-data-simulator.readthedocs.io/en/latest/yaml_generator.html
    # -- https://github.com/spacetelescope/mirage/blob/master/mirage/yaml/yaml_generator.py
    if not os.path.isfile(observation_list_final_file): 
        # Run yaml_generator.SimInput to create observation table csv file
        yam = yaml_generator.SimInput(
            input_xml = xml_file, 
            pointing_file = pointing_file,
            catalogs = catalogs, 
            cosmic_rays = DEFAULT_COSMIC_RAYS,
            background = background, 
            roll_angle = pa_v3,
            dates = dates, 
            reffile_defaults = 'crds_full_name',
            add_ghosts = False, 
            convolve_ghosts_with_psf = False,
            verbose = True, 
            output_dir = os.path.abspath(yaml_output_dir),
            simdata_output_dir = os.path.abspath(sim_output_dir),
            observation_list_file = os.path.abspath(observation_list_file), 
            datatype = 'linear, raw', # output both
            use_JWST_pipeline = True, 
            #offline = True, # cannot be offline otherwise linearized_darkfile will not be set
        )

        yam.use_linearized_darks = True
        
        #yam.add_catalogs()
        
        # Create updated observation table csv file and many yaml files
        yam.create_inputs()
    
    
    # Load observation_table
    observation_table = Table.read(observation_list_final_file)
    
    # Loop each simulation entry (single detector exposure)
    for i in range(len(observation_table)):
        sim_data_filename = observation_table['outputfits'][i]
        sim_data_filepath = os.path.join(sim_output_dir, sim_data_filename)
        yaml_filename = observation_table['yamlfile'][i]
        yaml_file = os.path.join(yaml_output_dir, yaml_filename)
        yaml_name = re.sub(r'\.yaml$', r'', yaml_filename)
        
        # Check if a filter has been specified
        if filter_name is not None:
            if not check_yamlfile_matches_the_filter(yaml_file, filter_name):
                if verbose:
                    logger.info('*** Skipping observation {!r} ({}/{}) because of non-matched filter'.format(
                        yaml_name, i+1, len(observation_table)).ljust(96) + ' ***')
                continue
        
        # Check if use has set an only_dataset
        if only_dataset is None:
            only_dataset = []
        if isinstance(only_dataset, str):
            only_dataset = [only_dataset]
        if len(only_dataset) > 0:
            check_okay = False
            for dataset_name in only_dataset:
                if dataset_name in os.path.basename(yaml_file):
                    check_okay = True
                    break
            if not check_okay:
                if verbose:
                    logger.info('*** Skipping observation {!r} ({}/{}) because of non-matched dataset name'.format(
                        yaml_name, i+1, len(observation_table)).ljust(96) + ' ***')
                continue
        
        # Check if add mosaic image as extended image
        if mosaic_file is not None:
        
            if verbose:
                logger.info('*'*100)
                logger.info('*** Processing observation {!r} ({}/{})'.format(
                    yaml_name, i+1, len(observation_table)).ljust(96) + ' ***')
                logger.info('*'*100)
            
            # Check yaml file with user input mosaic image as extended image
            yaml_ext_file = os.path.join(yaml_output_dir, yaml_name+'_ext.yaml')
            if os.path.isfile(yaml_ext_file):
                if overwrite_simprep:
                    shutil.move(yaml_ext_file, yaml_ext_file+'.backup')
                else:
                    if verbose:
                        logger.info('Found existing yaml file {!r} and overwrite is set to False.'.format(yaml_ext_file))
            
            # Check ext img file, if missing, recreate "*_ext.yaml" anyway
            extended_img_file = os.path.join(sim_output_dir, yaml_name+'_ext.fits')
            if not os.path.isfile(extended_img_file) and os.path.isfile(yaml_ext_file):
                shutil.move(yaml_ext_file, yaml_ext_file+'.backup')
            
            # Check ext cat file, if missing, recreate "*_ext.yaml" anyway
            extended_cat_file = os.path.join(sim_output_dir, yaml_name+'_ext.cat')
            if not os.path.isfile(extended_cat_file) and os.path.isfile(yaml_ext_file):
                shutil.move(yaml_ext_file, yaml_ext_file+'.backup')
            
            # Process yaml file to use mosaic image as extended image
            if not os.path.isfile(yaml_ext_file):
                # instrument = observation_table['Instrument'][i]
                # if instrument == 'NIRCAM':
                #     aperture = observation_table['aperture'][i]
                #     if aperture.startswith('NRCA5') or aperture.startswith('NRCB5'):
                #         filter_name = observation_table['LongFilter'][i]
                #     else:
                #         filter_name = observation_table['ShortFilter'][i]
                # else:
                #     raise NotImplementedError()
                
                # recenter?
                if mosaic_center is not None:
                    recenter = True
                    center_ra = mosaic_center.ra.deg
                    center_dec = mosaic_center.dec.deg
                else:
                    recenter = False
                    center_ra = None
                    center_dec = None
            
                # Resample mosaic image
                extended_img = resample_mosaic_image(
                    yaml_file = yaml_file, 
                    mosaic_file = mosaic_file, 
                    mosaic_psf = mosaic_psf, 
                    output_dir = sim_output_dir, 
                    output_name = yaml_name+'_ext', 
                    recenter = recenter, 
                    center_ra = center_ra, 
                    center_dec = center_dec, 
                    overwrite = overwrite_simprep, 
                )
                
                # Get RA Dec PAV3
                ra, dec, pav3 = get_ra_dec_pav3_from_yamlfile(yaml_file)
                
                # Set extendedscale if needed (TODO)
                extended_scale = 1.0
                
                # Set starting_index
                starting_index = 1
                if star_catalog is not None:
                    starting_index += get_catalog_row_number(star_catalog)
                if galaxy_catalog is not None:
                    starting_index += get_catalog_row_number(galaxy_catalog)

                # Prepare extended source catalog
                extended_catalog = catalog_generator.ExtendedCatalog(
                    filenames = [extended_img],
                    ra = [ra],
                    dec = [dec],
                    #position_angle = [pav3],
                    position_angle = [-pav3], # the rotation is already corrected to detector frame by resample_mosaic_image()!
                    # to cancel out the rotation in "mirage/seed_image/catalog_seed_image.py" `calc_x_position_angle`, 
                    # we need to set -pav3 here!
                    starting_index = starting_index, 
                )
                extended_catalog.add_magnitude_column(
                    ['None'], # a generic magnitude
                    #instrument = instrument,
                    #filter_name = filter_name, # comment out to add a generic one, see "mirage/catalogs/catalog_generator.py"
                )
                #extended_catalog.add_magnitude_column( # we can add multiple
                #    [str(magnitude)],
                #    instrument = instrument,
                #    filter_name = 'F444W',
                #)
                
                # Write extended source catalog to simdata dir
                extended_catalog_file = os.path.join(sim_output_dir, yaml_name+'_ext.cat')
                if os.path.isfile(extended_catalog_file):
                    shutil.move(extended_catalog_file, extended_catalog_file+'.backup')
                extended_catalog.save(extended_catalog_file)
                
                # Update Yaml file to use the extended image
                update_yamlfile_with_extended_catalog(
                    yaml_file, 
                    extended_catalog_file, 
                    extended_scale, 
                    yaml_ext_file, 
                )
            
            yaml_file = yaml_ext_file
        
        
        # match zeropoints by updating flux_cal file, see below.
        # see also "go-mirage-sim-one-yaml-file.py". 
        
        # Run actual simulation to generate uncal files
        if not skip_simdata:
            
            # Check sim data output file
            if os.path.isfile(sim_data_filepath):
                if overwrite_simdata:
                    shutil.move(sim_data_filepath, sim_data_filepath+'.backup')
                else:
                    if verbose:
                        logger.info('Found existing data file {!r} and overwrite is set to False. Skipping!'.format(sim_data_filepath))
                        continue
            
            # Check sim data output dir
            if not os.path.isdir(sim_output_dir):
                os.makedirs(sim_output_dir)
            
            # load yaml_file
            with open(yaml_file, 'r') as fp:
                yamldict = yaml.safe_load(fp)
            
            # Update yaml file in the output directory, we will use that instead of the input yaml file for the simulation
            original_yaml_file = yaml_file
            new_yaml_file = os.path.join(sim_output_dir, os.path.basename(yaml_file))
            
            # Update yaml dict
            #yamldict['Output']['directory'] = sim_output_dir
            #yamldict['Output']['file'] = output_file
            
            if new_psf_library is not None:
                # defaults 
                # psfpath: $MIRAGE_DATA/nircam/gridded_psf_library/
                # psf_wing_threshold_file: $CONDA_PREFIX/lib/python3.9/site-packages/mirage/config/nircam_psf_wing_rate_thresholds.txt
                yamldict['simSignals']['psfpath'] = new_psf_library
                if verbose:
                    logger.info('Using new PSF library {!r}'.format(new_psf_library))
            
            if new_psf_wing_threshold is not None:
                yamldict['simSignals']['psf_wing_threshold_file'] = new_psf_wing_threshold
                if verbose:
                    logger.info('Using new PSF wing threshold file {!r}'.format(new_psf_wing_threshold))
            
            # match zeropoints by creating a new flux_cal that matches the CRDS photom reference file
            if match_zeropoints:
                instrument_name = yamldict['Inst']['instrument']
                filter_name = yamldict['Readout']['filter']
                pupil_name = yamldict['Readout']['pupil']
                array_name = yamldict['Readout']['array_name']
                photom_file = yamldict['Reffiles']['photom']
                flux_cal_file = yamldict['Reffiles']['flux_cal']
                # read photom_file to get pixar_sr
                pixar_sr = None
                photmjsr = None
                with fits.open(photom_file) as photom_hdul:
                    pixar_sr = photom_hdul[0].header['PIXAR_SR'] # u.sr
                    for irow in range(len(photom_hdul[1].data)):
                        photom_filter, photom_pupil, photom_photmjsr, photom_uncertainty = photom_hdul[1].data[irow]
                        if photom_filter == filter_name and photom_pupil == pupil_name:
                            photmjsr = photom_photmjsr
                            break
                if photmjsr is None:
                    logger.error('Error! Could not find filter {} pupil {} in photom file {}'.format(
                        filter_name, pupil_name, photom_file))
                    raise Exception('Error! Could not find filter {} pupil {} in photom file {}'.format(
                        filter_name, pupil_name, photom_file))
                
                # read flux_cal file to get Pivot_wave
                pivot_wave = None
                flux_cal_headers = None
                flux_cal_dict = None
                with open(flux_cal_file, 'r') as fp:
                    line_number = 0
                    for line_text in fp:
                        if line_text.strip() == '':
                            continue
                        line_number += 1
                        if line_number == 1:
                            flux_cal_headers = line_text.split()
                        else:
                            line_split = line_text.split()
                            if len(line_split) == len(flux_cal_headers):
                                flux_cal_dict = OrderedDict(zip(flux_cal_headers, line_split))
                                if flux_cal_dict['Filter'] == filter_name and \
                                    flux_cal_dict['Pupil'] == pupil_name and \
                                    flux_cal_dict['Module'] == array_name[3] and \
                                    flux_cal_dict['Detector'] == array_name[0:5]:
                                    pivot_wave = float(flux_cal_dict['Pivot_wave'])
                                    break
                if pivot_wave is None:
                    logger.error('Error! Could not find filter {} pupil {} module {} detector {} in flux_cal file {}'.format(
                        filter_name, pupil_name, array_name[3], array_name[0:5], flux_cal_file))
                    raise Exception('Error! Could not find filter {} pupil {} module {} detector {} in flux_cal file {}'.format(
                        filter_name, pupil_name, array_name[3], array_name[0:5], flux_cal_file))
                flux_cal_dict['VEGAMAG'] = float(flux_cal_dict['VEGAMAG'])
                flux_cal_dict['ABMAG'] = float(flux_cal_dict['ABMAG'])
                flux_cal_dict['STMAG'] = float(flux_cal_dict['STMAG'])
                flux_cal_dict['PHOTFLAM'] = float(flux_cal_dict['PHOTFLAM'])
                flux_cal_dict['PHOTFNU'] = float(flux_cal_dict['PHOTFNU'])
                flux_cal_dict['Pivot_wave'] = float(flux_cal_dict['Pivot_wave'])
                
                # prepare new flux_cal file
                old_flux_cal_file = os.path.join(sim_output_dir, os.path.splitext(os.path.basename(yaml_file))[0] + '_flux_cal_mirage.txt')
                with open(old_flux_cal_file, 'w') as fp:
                    fp.write(' '.join(flux_cal_headers)+'\n')
                    fp.write(' '.join([str(t) for t in flux_cal_dict.values()])+'\n')
                    logger.info('Extracted MIRAGE flux_cal file {!r} to temporary file {!r}'.format(flux_cal_file, old_flux_cal_file))
                
                # Filter Pupil Module Detector VEGAMAG ABMAG STMAG PHOTFLAM PHOTFNU Pivot_wave
                ABMAG = ((photmjsr * u.MJy/u.sr) * (pixar_sr * u.sr)).to(u.ABmag)
                logger.info('ABMAG: {} -> {}'.format(flux_cal_dict['ABMAG'], ABMAG))
                PHOTFNU = ABMAG.to(u.erg/u.s/u.cm**2/u.Hz)
                logger.info('PHOTFNU: {} -> {}'.format(flux_cal_dict['PHOTFNU'], PHOTFNU))
                Pivot_wave = pivot_wave * u.um
                logger.info('Pivot_wave: {}'.format(Pivot_wave))
                PHOTFLAM = PHOTFNU * (const.c.cgs/Pivot_wave.cgs).to(u.Hz) / (Pivot_wave).to(u.AA)
                logger.info('PHOTFLAM: {} -> {}'.format(flux_cal_dict['PHOTFLAM'], PHOTFLAM))
                STMAG = ABMAG.to(u.STmag, u.spectral_density(Pivot_wave)) #TODO error?
                logger.info('STMAG: {} -> {}'.format(flux_cal_dict['STMAG'], STMAG))
                VEGAMAG = float(flux_cal_dict['VEGAMAG']) + (ABMAG.value - float(flux_cal_dict['ABMAG']))
                logger.info('VEGAMAG: {} -> {}'.format(flux_cal_dict['VEGAMAG'], VEGAMAG))
                flux_cal_dict['VEGAMAG'] = VEGAMAG
                flux_cal_dict['ABMAG'] = ABMAG.value
                flux_cal_dict['STMAG'] = STMAG.value
                flux_cal_dict['PHOTFLAM'] = PHOTFLAM.value
                flux_cal_dict['PHOTFNU'] = PHOTFNU.value
                
                new_flux_cal_file = os.path.join(sim_output_dir, os.path.splitext(os.path.basename(yaml_file))[0] + '_flux_cal_crds.txt')
                with open(new_flux_cal_file, 'w') as fp:
                    fp.write(' '.join(flux_cal_headers)+'\n')
                    fp.write(' '.join([str(t) for t in flux_cal_dict.values()])+'\n')
                    logger.info('Converted CRDS photom file {!r} to new flux_cal file {!r}'.format(photom_file, new_flux_cal_file))
            
                yamldict['Reffiles']['flux_cal'] = new_flux_cal_file
                if verbose:
                    logger.info('Using new CRDS-based flux_cal file {!r}'.format(new_flux_cal_file))
            
            
            # Write yaml dict to disk
            if os.path.isfile(new_yaml_file):
                shutil.move(new_yaml_file, new_yaml_file+'.backup')
            with open(new_yaml_file, 'w') as fp:
                yaml.dump(yamldict, fp)
            if verbose:
                logger.info('Updated {!r}'.format(new_yaml_file))
            
            
            # Run simulation to generate data file
            if not os.path.isfile(sim_data_filepath):
                if verbose:
                    logger.info('*'*100)
                    logger.info('*** Running Actual Simulation for {!r}'.format(
                        new_yaml_file).ljust(96) + ' ***')
                    logger.info('*'*100)
                m = ImgSim() # override_dark=dark
                m.paramfile = new_yaml_file
                m.create()
    
    
    # return observation_table
    return observation_table




if __name__ == '__main__':
    main()



