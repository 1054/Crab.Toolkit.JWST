#!/usr/bin/env python
# 
import os, sys, re, json, shutil
if "CRDS_PATH" not in os.environ:
    os.environ["CRDS_PATH"] = '/n17data/dzliu/Data/jwst_crds_cache'
if "MIRAGE_DATA" not in os.environ:
    os.environ["MIRAGE_DATA"] = '/n23data1/hjmcc/jwst/mirage/mirage_data'
if "CRDS_SERVER_URL" not in os.environ:
    os.environ["CRDS_SERVER_URL"] = 'https://jwst-crds.stsci.edu'

import astropy.units as u
import click
import numpy as np
import crds
import jwst
import mirage
from mirage import imaging_simulator
from mirage.seed_image import catalog_seed_image

import logging
logging.basicConfig(level='DEBUG')
logger = logging.getLogger('go-mirage-fetch-bkgdrate')





####################
### MAIN PROGRAM ###
####################

@click.command()
@click.argument('yaml_dir', type=click.Path(exists=True))
@click.option('--requery', is_flag=True, default=False)
@click.option('--verbose', is_flag=True, default=True)
def main(
        yaml_dir, 
        requery, 
        verbose, 
    ):

    if verbose:
        logger.info('crds version: {}'.format(crds.__version__))
        logger.info('jwst version: {}'.format(jwst.__version__))
        logger.info('mirage version: {}'.format(mirage.__version__))
    
    if yaml_dir.endswith('/'):
        yaml_dir = yaml_dir.rstrip('/')
    
    list_files = glob.glob(yaml_dir+'/**/jw*.yaml')
    if len(list_files) == 0:
        raise Exception('Error! No file found: ' + yaml_dir+'/**/jw*.yaml')
    list_files = sorted(list_files)
    
    cache_bkgdrates = {}
    if os.path.isfile(yaml_dir+'/bkgdrates.json'):
        with open(yaml_dir+'/bkgdrates.json', 'r') as fp:
            cache_bkgdrates = json.load(fp)
        shutil.copy2(
            yaml_dir+'/bkgdrates.json', 
            yaml_dir+'/bkgdrates.json.backup'
        )

    icount = 0
    for paramfile in list_files:
        
        if paramfile in cache_bkgdrates:
            if verbose:
                logger.info('Skipping ' + paramfile)
            continue
        
        with open(paramfile, 'r') as fp:
            yamlcontent = yaml.safe_load(fp)
            
        if np.isreal(yamlcontent['simSignals']['bkgdrate']):
            bkgdrate = yamlcontent['simSignals']['bkgdrate']
            if paramfile not in cache_bkgdrates:
                cache_bkgdrates[paramfile] = bkgdrate
            if verbose:
                logger.info('Skipping ' + paramfile)
            continue
        
        m = imaging_simulator.ImgSim()
        m.paramfile = paramfile

        cat = catalog_seed_image.Catalog_seed(offline=m.offline)
        cat.paramfile = m.paramfile
        cat.readParameterFile()
        cat.check_params()
        #cat.params['simSignals']['bkgdrate']
        
        bkgdrate = float(cat.params['simSignals']['bkgdrate'])
        
        yamlcontent['simSignals']['bkgdrate'] = bkgdrate
        
        cache_bkgdrates[paramfile] = bkgdrate
        
        if verbose:
            logger.info('Updating ' + paramfile + ' with bkgdrate ' + str(yamlcontent['simSignals']['bkgdrate']))
        
        shutil.copy2(paramfile, paramfile + '.backup')
        
        with open(paramfile, 'w') as fp:
            yaml.dump(yamlcontent, fp)
        
        if icount > 0 and icount % 10 == 0:
            with open(yaml_dir+'/bkgdrates.json', 'w') as fp:
                json.dump(cache_bkgdrates, fp, indent=4)
                
        icount += 1


    # save the cached bkgrates into a json file

    with open(yaml_dir+'/bkgdrates.json', 'w') as fp:
        json.dump(cache_bkgdrates, fp, indent=4)

    if verbose:
        logger.info('Output to ' + yaml_dir+'/bkgdrates.json')

    if verbose:
        logger.info('All fixed.')




if __name__ == '__main__':
    main()



