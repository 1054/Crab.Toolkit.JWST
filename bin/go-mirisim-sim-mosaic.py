#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, time, json, yaml, asdf
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = os.path.expanduser('~/Data/JWST-CRDS') # '/n17data/dzliu/Data/jwst_crds_cache'
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = os.path.expanduser('~/Data/JWST-MIRAGE') # '/n23data1/hjmcc/jwst/mirage/mirage_data'
if "MIRISIM_ROOT" not in os.environ:
    os.environ["MIRISIM_ROOT"] = os.path.expanduser('~/Data/JWST-MIRISIM/mirisim') # os.path.expanduser('~/mirisim')
if "PYSYN_CDBS" not in os.environ:
    os.environ["PYSYN_CDBS"] = os.path.expanduser('~/Data/JWST-MIRISIM/mirisim/cdbs/') # os.path.expanduser('~/mirisim/cdbs')
if "CDP_DIR" not in os.environ:
    os.environ["CDP_DIR"] = os.path.expanduser('~/Data/mirisim_data/CDP')

sys.path.insert(0, os.path.expanduser('~/Cloud/Github/mirage'))
sys.path.insert(0, os.path.expanduser('~/Cloud/Github/mirisim_iap_fr'))

import click
import astropy.units as u
import numpy as np

from astropy.coordinates import SkyCoord, FK5
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from collections import OrderedDict

import pysiaf

import jwst
import jwst.assign_wcs
from jwst.assign_wcs.pointing import compute_roll_ref

import mirage # pip install --upgrade git+https://github.com/spacetelescope/mirage.git
from mirage.utils import siaf_interface
from mirage.utils.set_telescope_pointing_separated import compute_local_roll

import configobj # needed by mirisim
import pysynphot # needed by mirisim
import importlib_resources # needed by mirisim
import pysftp # needed by mirisim
import pySpecSim # needed by mirisim, copied from server, cannot be installed via pip
import mirisim
from mirisim import MiriSimulation
from mirisim.config_parser import SimConfig, SceneConfig, SimulatorConfig
from mirisim.imsim.ima import get_ima_subarray_to_skysim_transform
from mirisim.imsim.mirimImager import MirimImager
from mirisim.obssim import ObservationSimulation
from mirisim.obssim.pointing import Pointing
from mirisim.skysim import builder
from mirisim.skysim import scenemaker
from mirisim.skysim.catalogues import Catalogue
from mirisim.skysim.scenes import SkyScene, CompositeSkyScene
from mirisim.skysim.astrosources import Galaxy, Star
from mirisim.skysim.disks import SersicDisk
from mirisim.skysim.externalsources import Skycube, Skyimage

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirisim-sim-mosaic')


MIRI_SIAF = pysiaf.siaf.Siaf('MIRI')
MIRIM_ILLUM_SIAF = MIRI_SIAF.apertures['MIRIM_ILLUM']
MIRIM_CORONLYOT_SIAF = MIRI_SIAF.apertures['MIRIM_CORONLYOT']
MIRIM_ILLUM_V2REF = MIRIM_ILLUM_SIAF.V2Ref
MIRIM_ILLUM_V3REF = MIRIM_ILLUM_SIAF.V3Ref
MIRIM_ILLUM_V3IYANG = MIRIM_ILLUM_SIAF.V3IdlYAngle
MIRIM_ILLUM_VPARITY = MIRIM_ILLUM_SIAF.VIdlParity
logger.info("MIRIM_ILLUM_SIAF.V2Ref, V3Ref, V3IdlYAngle, VIdlParity = {}, {}, {}, {}".format( 
      MIRIM_ILLUM_SIAF.V2Ref, 
      MIRIM_ILLUM_SIAF.V3Ref, 
      MIRIM_ILLUM_SIAF.V3IdlYAngle, 
      MIRIM_ILLUM_VPARITY
))

DEFAULT_APT_XML_FILE = 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.xml'
DEFAULT_POINTING_FILE = 'apt_files/cosmosweb_revised_jun2022_onlyDEC2022.pointing'
#DEFAULT_MOSAIC_FILE = None # 'input_mosaic_images/dust_opa_u_lensed_H_EUC.fits' # 'hlsp_candels_hst_acs_gs-tot-sect23_f814w_v1.0_drz.fits'
#DEFAULT_STAR_CATALOG = 'input_catalogs/ptsrc_pointings_BEST_sw_tot.cat'
DEFAULT_FILTER = 'F770W'
DEFAULT_DATES = '2023-01-01'
DEFAULT_PA_V3 = 293.09730273
DEFAULT_SIM_OUTPUT_DIR = 'mirisim_output'




# Define A Catalogue Class

class MyCatalog(Catalogue):
    """ A subclass of mirisim.skysim.catalogues.Catalogue.
    """
    def __init__(self, filename, pointing:mirisim.obssim.pointing.Pointing):
        """ The input ra0, dec0 should be the MIRI instrument center coordinate.
        """
        self.catalogue = None
        v2_ref, v3_ref = pointing.get_v2v3_ref() # MIRIM_ILLUM center's (v2,v3) coordinate
        v2_off, v3_off = pointing.get_v2v3_offset_actual() # for this particular dither, arcsec
        self.roll = pointing.pa # this pa is the roll angle of MIRIM_ILLUM Ideal+Y to North
        self.v3_pa = -pointing.pa+MIRIM_ILLUM_V3IYANG # this telescope V3 PA
        self.v1_ra, self.v1_dec = pointing.get_radec_v1()
        self.dec0 = self.v1_dec + (v3_ref - v3_off)/3600.0 # MIRIM_ILLUM center's RA Dec
        self.ra0 = self.v1_ra + (v2_ref - v2_off)/3600.0/np.cos(np.deg2rad(self.v1_dec)) # MIRIM_ILLUM center's RA Dec
        super(MyCatalog, self).__init__(filename)
    
    def readcatalogue(self, filename, 
            col_id = 'ID', 
            col_ra = 'RA', 
            col_dec = 'Dec', 
            col_n = 'n', 
            col_re = 're', 
            col_q = 'q', 
            col_pa = 'pa', 
            **kwargs
        ):
        """Reads the catalogue.
        """
        self.col_id = col_id
        self.col_ra = col_ra
        self.col_dec = col_dec
        self.col_n = col_n
        self.col_re = col_re
        self.col_q = col_q
        self.col_pa = col_pa
        self.catalogue = Table.read(filename, **kwargs)
        assert self.col_id in self.catalogue.colnames
        assert self.col_ra in self.catalogue.colnames
        assert self.col_dec in self.catalogue.colnames
        assert self.col_n in self.catalogue.colnames
        assert self.col_re in self.catalogue.colnames
        assert self.col_q in self.catalogue.colnames
        assert self.col_pa in self.catalogue.colnames
    
    def parsecatalogue(self,):
        """returns CompositeScene object.
        """
        for i in range(len(self.catalogue)):
            kwargs = {}
            ra = self.catalogue[self.col_ra][i]
            dec = self.catalogue[self.col_dec][i]
            kwargs['Cen'] = [-(ra-self.ra0)*np.cos(np.deg2rad(self.dec0))*3600.0, 
                             (dec-self.dec0)*3600.0] # offset to detector MIRIM_ILLUM center, in arcsec
            kwargs['n'] = self.catalogue[self.col_n][i]
            kwargs['re'] = self.catalogue[self.col_re][i]
            kwargs['q'] = self.catalogue[self.col_q][i] # see mirisim.skysim.disks.Disk.__init__
            kwargs['pa'] = self.catalogue[self.col_pa][i]
            newscene = Galaxy(**kwargs)
            newscene.name = self.catalogue[self.col_id][i]
            if i == 0:
                scene = CompositeSkyScene(sources=[newscene], operations=[])
            else:
                scene += CompositeSkyScene(self.catalogue)
        
        return scene



# define function to read pointing table
def read_pointing_table(pointing_file):
    header = None
    table_data = []
    logger.info('Reading pointing file: {!r}'.format(pointing_file))
    with open(pointing_file, 'r') as fp:
        next_is_header_line = False
        next_is_data_line = False
        for line in fp:
            line = line.strip()
            if line == '' or line.startswith('#') or line.startswith('='):
                next_is_header_line = False
                next_is_data_line = False
            elif line.startswith('* Observation '):
                this_obsnum = int(line.replace('* Observation ', ''))
            elif line.startswith('** Visit '):
                this_visitnum = int(line.replace('** Visit ', '').split(':')[1])
                next_is_header_line = True
            elif next_is_header_line:
                if header is None:
                    line2 = line
                    line2 = line2.replace('Aperture Name', 'Aperture_Name')
                    line2 = line2.replace('Target', 'Target_ID Target')
                    line2 = 'obsnum' + ' ' + 'visitnum' + ' ' + line2 + '  is_base'
                    header = line2.split()
                next_is_header_line = False
                next_is_data_line = True
            elif next_is_data_line:
                line2 = line
                if line2.find(' (base)') >= 0:
                    line2 = line2.replace(' (base)', '  1')
                line2 = str(this_obsnum) + ' ' + str(this_visitnum) + ' ' + line2 + '  0  0'
                data = line2.split()
                print(line2)
                print(len(header), len(data))
                table_data.append(data[0:len(header)])
    table_dict = OrderedDict()
    for i in range(len(header)):
        table_dict[header[i]] = [t[i] for t in table_data]
    table = Table(table_dict)
    print('Pointing Table:\n', table)
    return table





####################
### MAIN PROGRAM ###
####################

@click.command()
@click.option('--xml-file', type=click.Path(exists=True), default=DEFAULT_APT_XML_FILE)
@click.option('--pointing-file', type=click.Path(exists=True), default=DEFAULT_POINTING_FILE)
@click.option('--mosaic-file', type=click.Path(exists=True), default=None)
@click.option('--star-catalog', type=click.Path(exists=True), default=None)
@click.option('--galaxy-catalog', type=click.Path(exists=True), default=None)
@click.option('--filter', 'filter_name', type=str, default=DEFAULT_FILTER)
@click.option('--dates', type=str, default=DEFAULT_DATES)
@click.option('--pa-v3', type=float, default=DEFAULT_PA_V3)
@click.option('--sim-output-dir', type=click.Path(exists=False), default=DEFAULT_SIM_OUTPUT_DIR)
@click.option('--overwrite', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        xml_file, 
        pointing_file, 
        mosaic_file, 
        star_catalog, 
        galaxy_catalog, 
        filter_name, 
        dates, 
        pa_v3,
        sim_output_dir,
        overwrite, 
        verbose, 
    ):
    
    if verbose:
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirisim version: {}'.format(mirisim.__version__))
    
    # 
    
    # #MiriSimulation.generate_configfiles()
    # #mysim = MiriSimulation.from_configfiles('simulation.ini')
    
    pointing_table = read_pointing_table(pointing_file)
    
    raise NotImplementedError()
    
    
    # Load config
    sim_config = SimConfig('ima_simulation.ini')
    scene_config = SceneConfig('scene.ini')
    simulator_config = SimulatorConfig('simulator.ini')
    
    # mysim = MiriSimulation(sim_config, scene_config, simulator_config)
    # mysim.run()
    
    
    if not os.path.isdir(sim_output_dir):
        os.makedirs(sim_output_dir)
    
    path_out = sim_output_dir
    path_cdp = os.environ['CDP_DIR']
    
    
    # Create simulation object
    simulation = ObservationSimulation(
        sim_config, scene_config, simulator_config,
        path_out, path_cdp,
    )
    
    
    # Assign the correct V1_RA V1_Dec V3_PA to each ExposureEvent
    telescope_v1_ra = +149.92961
    telescope_v1_dec = +2.43173
    telescope_v3_pa = 293.09730273
    telescope_v2_ref = +22.442
    telescope_v3_ref = -495.969
    
    simulation.ra = telescope_v1_ra # will go into pointing.radec_ref
    simulation.dec = telescope_v1_dec # will go into pointing.radec_ref
    simulation.pa = telescope_v3_pa
    
    for i in range(len(simulation.events)):
        if isinstance(simulation.events[i], mirisim.obssim.event.ExposureEvent):
            # # following mirage.utils.set_telescope_pointing_separated.add_wcs(), we compute local_roll
            # local_roll = compute_local_roll(
            #     telescope_v3_pa, 
            #     telescope_v1_ra, telescope_v1_dec, 
            #     telescope_v2_ref, telescope_v3_ref, 
            #     #MIRIM_ILLUM_V2REF, MIRIM_ILLUM_V3REF
            # )
            # following mirage.apt.apt_inputs(), use siaf_interface to get attitude_matrix
            local_roll, attitude_matrix, fullframesize, subarray_boundaries = \
                siaf_interface.get_siaf_information(
                    MIRI_SIAF, 'MIRIM_ILLUM', 
                    telescope_v1_ra, telescope_v1_dec, telescope_v3_pa,
                    v2_arcsec = telescope_v2_ref, v3_arcsec = telescope_v3_ref, # v2 v3 ref for the ra dec
                )
            # set attitude_matrix
            MIRIM_ILLUM_SIAF.set_attitude_matrix(attitude_matrix)
            # compute MIRIM_ILLUM center RA Dec
            ra_ref, dec_ref = MIRIM_ILLUM_SIAF.tel_to_sky(MIRIM_ILLUM_SIAF.V2Ref, MIRIM_ILLUM_SIAF.V3Ref)
            logger.info('simulation.events[{}], telescope ra, dec: ({}, {}), v2 v3: ({}, {}), pa: {}'.format(
                i, telescope_v1_ra, telescope_v1_dec, telescope_v2_ref, telescope_v3_ref, telescope_v3_pa, 
            ))
            logger.info('simulation.events[{}], MIRIM_ILLUM ra, dec: ({}, {}), v2 v3: ({}, {}), pa: {}'.format(
                i, ra_ref, dec_ref, MIRIM_ILLUM_SIAF.V2Ref, MIRIM_ILLUM_SIAF.V3Ref, local_roll,
            ))
            # pointing PA is the Rotation2D(pa) in mirisim.wcs.v2v3SkySimTransform(), 
            # it should be the angle from MIRIM_ILLUM +Y to North, counterclockwise.
            # local_roll is the angle from North to V3, counterclockwise
            # V3IdlYAngle is the angle from V3 to MIRIM_ILLUM +Y, counterclockwise
            simulation.events[i].pointing.pa = local_roll # roll angle of MIRIM_ILLUM
            simulation.events[i].pointing.set_radec_v1(telescope_v1_ra, telescope_v1_dec) # stored as pointing.radec_v1 for later fixing wcs
            #simulation.events[i].pointing.set_radec_pointing_ref(ra_ref, dec_ref) # stored as pointing.radec_ref
            #simulation.events[i].pointing.set_v2v3_ref(telescope_v2_ref, telescope_v3_ref) # stored as pointing.v2v3_ref
            #simulation.events[i].pointing.set_v2v3_pointing_ref(MIRIM_ILLUM_SIAF.V2Ref, MIRIM_ILLUM_SIAF.V3Ref) # stored as pointing.v2v3_pref
    
    
    # Run the simulation
    simulation.run()
    
    
    # Post-process illum_models, add WCS.
    for i in range(len(simulation.events)):
        if isinstance(simulation.events[i], mirisim.obssim.event.ExposureEvent):
            # DZLIU fixing WCS of illum_models
            file_base = os.path.join(
                simulation.events[i].path_out, 
                simulation.events[i].exposures['IMA'].path_illum_models, 
                simulation.events[i].exposures['IMA'].fn_illum_model,
            )
            # following mirisim.imsim.miriImager.MirimImager()
            #pa = simulation.events[i].pointing.pa # the same as pointing_v3_PA
            subarray = simulation.events[i].exposures['IMA'].subarray
            filterName = simulation.events[i].exposures['IMA'].filtername
            #subarray = sim_config['Integration_and_patterns']['IMA_configuration']['ReadDetect']
            #filterName = sim_config['Integration_and_patterns']['IMA_configuration']['filter']
            v2_ref, v3_ref = simulation.events[i].pointing.get_v2v3_ref() # see mirisim.wcs.v2v3SkySimTransform
            v2_off, v3_off = simulation.events[i].pointing.get_v2v3_offset_actual() # see above
            v1_ra, v1_dec = simulation.events[i].pointing.get_radec_v1()
            dec = v1_dec + (v3_ref - v3_off) # MIRIM_ILLUM center's RA Dec
            ra = v1_ra + (v2_ref - v2_off)/3600.0/np.cos(np.deg2rad(v1_dec)) # MIRIM_ILLUM center's RA Dec
            pa = simulation.events[i].pointing.pa # roll angle from MIRIM_ILLUM +Y to North, counterclockwise
            v3_y_angle = MIRIM_ILLUM_V3IYANG # angle from V3 to MIRIM_ILLUM +Y, counterclockwise
            v3_pa = -pa + MIRIM_ILLUM_V3IYANG # V3 PA, angle from North to V3, counterclockwise
            roll_ref = compute_roll_ref(
                0., 0., v3_pa, ra, dec, v2_ref, v3_ref,
            ) # PA of V3 at V2_REF, V3_REF; equaling V3 PA if (V2_REF, V3_REF) == (0.0, 0.0)
            all_transforms = get_ima_subarray_to_skysim_transform(
                subarray, filterName,
                v2_ref, v3_ref, v2_off, v3_off, pa,
                simulator_config = simulator_config,
            )
            # following mirisim.imsim.miriImager.get_points_from_scene()
            sky_to_colrow = all_transforms.get_transform('skysim_radec', 'fullarray_colrow')
            col, row = sky_to_colrow(0.0, 0.0)
            file_in = file_base + '_MIRIMAGE_' + filterName + '.fits'
            file_out = file_base + '_MIRIMAGE_' + filterName + '_fixedwcs.fits'
            with fits.open(file_in) as hdul:
                hdul[1].header['CRVAL1'] = ra
                hdul[1].header['CRVAL2'] = dec
                hdul[1].header['CRPIX1'] = col
                hdul[1].header['CRPIX2'] = row
                hdul[1].header['PC1_1'] = -np.cos(np.deg2rad(-roll_ref))
                hdul[1].header['PC1_2'] = -np.sin(np.deg2rad(-roll_ref))
                hdul[1].header['PC2_1'] = -np.sin(np.deg2rad(-roll_ref))
                hdul[1].header['PC2_2'] = np.cos(np.deg2rad(-roll_ref))
                hdul[1].header['V3YANGLE'] = v3_y_angle
                hdul[1].header['ROLL_REF'] = roll_ref
                hdul[1].header['RA_V1'] = ra
                hdul[1].header['DEC_V1'] = dec
                hdul[1].header['PA_V3'] = v3_pa
                hdul.writeto(file_out, overwrite=True)
                print('Output to {!r}'.format(file_out))




if __name__ == '__main__':
    main()
