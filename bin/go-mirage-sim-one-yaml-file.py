#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = '/n17data/dzliu/Data/jwst_crds_cache'
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = '/n23data1/hjmcc/jwst/mirage/mirage_data'
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


DEFAULT_NEW_PSF_LIBRARY = None # '/n17data/dzliu/Work/JWST-MIRAGE-Simulation/20221008_make_mirage_psf_library/mirage_data/nircam/gridded_psf_library'
DEFAULT_NEW_PSF_WING_THRESHOLD = None # '/n17data/dzliu/Work/JWST-MIRAGE-Simulation/20221008_make_mirage_psf_library/custom_nircam_psf_wing_rate_thresholds.txt'




####################
### MAIN PROGRAM ###
####################

@click.command()
@click.argument('yaml_file', type=click.Path(exists=True))
@click.option('--output-dir', type=click.Path(exists=False), default=None, help='Can specifiy a new output directory.')
@click.option('--output-file', type=str, default=None, help='Can specifiy a new output filename.')
@click.option('--new-psf-library', type=click.Path(exists=True), default=DEFAULT_NEW_PSF_LIBRARY, help='New psf library directory.')
@click.option('--new-psf-wing-threshold', type=click.Path(exists=True), default=DEFAULT_NEW_PSF_LIBRARY, help='New psf library directory.')
@click.option('--overwrite', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        yaml_file, 
        output_dir, 
        output_file, 
        new_psf_library,
        new_psf_wing_threshold, 
        overwrite, 
        verbose, 
    ):
    
    if verbose:
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    # get sim_data_filepath
    with open(yaml_file, 'r') as fp:
        yamldict = yaml.safe_load(fp)
        if output_dir is None:
            output_dir = yamldict['Output']['directory']
        if output_file is None:
            output_file = yamldict['Output']['file']
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    sim_data_filepath = os.path.join(output_dir, output_file)

    # Check sim data output file
    if os.path.isfile(sim_data_filepath):
        if overwrite_simdata:
            shutil.move(sim_data_filepath, sim_data_filepath+'.backup')
        else:
            if verbose:
                logger.info('Found existing data file {!r} and overwrite is set to False. Skipping.'.format(sim_data_filepath))
                return
    
    # Update yaml file in the output directory, we will use that instead of the input yaml file for the simulation
    original_yaml_file = yaml_file
    new_yaml_file = os.path.join(output_dir, os.path.basename(yaml_file))
    
    # Update yaml dict
    yamldict['Output']['directory'] = output_dir
    yamldict['Output']['file'] = output_file
    
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
    
    # what about flux_cal
    # defaults
    # flux_cal: $CONDA_PREFIX/lib/python3.9/site-packages/mirage/config/NIRCam_zeropoints.list
    #   Filter Pupil Module Detector VEGAMAG ABMAG STMAG PHOTFLAM PHOTFNU Pivot_wave
    #   F277W CLEAR A NRCA5 25.00545910447429 27.312679347975052 30.817225766862414 1.710447281981006e-21 4.314669729926466e-31 2.7499784072587885
    # Mirage 
    #   mirage.seed_image.catalog_seed_image.get_point_source_list()
    #     countrate = utils.magnitude_to_countrate(..., magsys, mag, photfnu=self.photfnu, ...)
    #   mirage.utils.utils.magnitude_to_countrate()
    #     if magsys.lower() == 'abmag':
    #       return count_scale * (10**((mag + 48.599934378) / -2.5) / photfnu)
    #   so if ABMAG = 27.312679347975052, PHOTFNU = 4.314669729926466e-31,
    #     (10**((27.312679347975052 + 48.599934378) / -2.5) / 4.314669729926466e-31) = 1.000
    #     PHOTFNU = 10**((ABMAG + 48.599934378) / -2.5)
    #     PHOTFNU = (ABMAG * u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz)
    #     Pivot_wave = 2.7499784072587885 * u.um
    #     PHOTFLAM = PHOTFNU * (const.c.cgs/Pivot_wave.cgs).to(u.Hz) / Pivot_wave.to(u.AA)
    #     STMAG = (ABMAG * u.ABmag).to(u.STmag, u.spectral_density(Pivot_wave))
    #     VEGAMAG = -2.5*np.log10(PHOTFNU.value) - 48.6
    #     VEGAMAG = PHOTFNU.to(u.ABmag).value
    #   
    # Max: 
    #   mag_zp_MJy = (mag_zp* u.ABmag).to(u.MJy)/pixar_sr
    #   pixar_sr   =  2.31E-14 for the sw
    #   pixar_sr   =  9.31E-14 for the lw
    #   --> (27.312679347975052 * u.ABmag).to(u.MJy) / (9.31E-14 * u.sr) = 0.46342406 MJy/sr
    # CRDS: 
    #   [2022-11-10 19:48:45] Prefetch for PHOTOM reference file is '/home/dzliu/jwst_crds_cache/references/jwst/nircam/jwst_nircam_photom_0111.fits'.
    #   hdul[0].header
    #   PIXAR_SR = 9.31E-14 / Nominal pixel area in steradians
    #   hdul[1].data
    #   [('filter', 'S12'), ('pupil', 'S12'), ('photmjsr', '>f4'), ('uncertainty', '>f4')]
    #   ('F277W', 'CLEAR',  0.4908   , 0.09816   ),
    #   ('F277W', 'MASKRND',  2.7386296, 0.5477259 ),
    #   ('F277W', 'MASKBAR',  2.7386296, 0.5477259 ),
    #   --> ((0.4908 * u.MJy/u.sr) * (9.31E-14 * u.sr)).to(u.ABmag)   --> 27.25036441 mag(AB)
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
                        pivot_wave = flux_cal_dict['Pivot_wave']
                        break
    if pivot_wave is None:
        logger.error('Error! Could not find filter {} pupil {} module {} detector {} in flux_cal file {}'.format(
            filter_name, pupil_name, array_name[3], array_name[0:5], flux_cal_file))
        raise Exception('Error! Could not find filter {} pupil {} module {} detector {} in flux_cal file {}'.format(
            filter_name, pupil_name, array_name[3], array_name[0:5], flux_cal_file))
    
    # prepare new flux_cal file
    old_flux_cal_file = os.path.join(output_dir, os.path.splitext(os.path.basename(yaml_file))[0] + '_flux_cal_old.txt')
    with open(old_flux_cal_file, 'w') as fp:
        fp.write(' '.join(flux_cal_headers)+'\n')
        fp.write(' '.join([str(t) for t in flux_cal_dict.values()])+'\n')
    
    # Filter Pupil Module Detector VEGAMAG ABMAG STMAG PHOTFLAM PHOTFNU Pivot_wave
    ABMAG = ((photmjsr * u.MJy/u.sr) * (9.31E-14 * u.sr)).to(u.ABmag)
    PHOTFNU = ABMAG.to(u.erg/u.s/u.cm**2/u.Hz)
    Pivot_wave = pivot_wave * u.um
    PHOTFLAM = PHOTFNU * (const.c.cgs/Pivot_wave.cgs).to(u.Hz) / (Pivot_wave).to(u.AA)
    STMAG = ABMAG.to(u.STmag, u.spectral_density(Pivot_wave))
    #VEGAMAG = -2.5*np.log10(PHOTFNU.value) - 48.6
    VEGAMAG = PHOTFNU.to(u.ABmag).value
    flux_cal_dict['VEGAMAG'] = VEGAMAG
    flux_cal_dict['ABMAG'] = ABMAG.value
    flux_cal_dict['STMAG'] = STMAG.value
    flux_cal_dict['PHOTFLAM'] = PHOTFLAM.value
    flux_cal_dict['PHOTFNU'] = PHOTFNU.value
    
    new_flux_cal_file = os.path.join(output_dir, os.path.splitext(os.path.basename(yaml_file))[0] + '_flux_cal.txt')
    with open(new_flux_cal_file, 'w') as fp:
        fp.write(' '.join(flux_cal_headers)+'\n')
        fp.write(' '.join([str(t) for t in flux_cal_dict.values()])+'\n')
    
    yamldict['Reffiles']['flux_cal'] = new_flux_cal_file
    if verbose:
        logger.info('Using new flux_cal file {!r}'.format(new_flux_cal_file))
    
    
    # Write yaml dict to disk
    if os.path.isfile(new_yaml_file):
        shutil.move(new_yaml_file, new_yaml_file+'.backup')
    with open(new_yaml_file, 'w') as fp:
        yaml.dump(yamldict, fp)
    if verbose:
        logger.info('Updated {!r}'.format(new_yaml_file))
    
    
    # Run simulation to generate data file
    if verbose:
        logger.info('*'*100)
        logger.info('*** Running Actual Simulation for {!r}'.format(
            new_yaml_file).ljust(96) + ' ***')
        logger.info('*'*100)
    m = ImgSim() # override_dark=dark
    m.paramfile = new_yaml_file
    m.create()
    
    if verbose:
        logger.info('*'*100)
        logger.info('*** Finished Simulation for {!r}'.format(
            new_yaml_file).ljust(96) + ' ***')
        logger.info('*'*100)




if __name__ == '__main__':
    main()



