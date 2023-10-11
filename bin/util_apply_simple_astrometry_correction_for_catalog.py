#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os, sys, re, shutil
import glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import click


# code name and version
CODE_NAME = 'util_apply_simple_astrometry_correction_for_catalog.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20231010'
CODE_HOMEPAGE = ''


# main
@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('output_catalog_file', type=click.Path(exists=False))
@click.option('--old-ra-dec', nargs=2, type=str, required=True, help='Input the old RA Dec for an anchor point, in degrees.')
@click.option('--new-ra-dec', nargs=2, type=str, required=True, help='Input the new RA Dec for an anchor point, in degrees.')
@click.option('--ra-column', type=str, default=None)
@click.option('--dec-column', type=str, default=None)
def main(
        input_catalog_file, 
        output_catalog_file, 
        old_ra_dec, 
        new_ra_dec,
        ra_column,
        dec_column,
    ):
    """
    This script will apply a simple astrometry correction to all RA Dec in the input catalog and save as the output catalog. 
    
    """
    
    print('Reading catalog %r'%(input_catalog_file))
    if re.match(r'.*\.(txt|dat|lis)$', input_catalog_file): 
        table = Table.read(input_catalog_file, format='ascii')
    elif re.match(r'.*\.ecsv$', input_catalog_file):
        table = Table.read(input_catalog_file, format='ascii.ecsv')
    else:
        table = Table.read(input_catalog_file)
    
    print('Reading columns ...')
    
    colRA = None
    colRA_list = ['RA', 'ALPHA_J2000', 'ra', 'sky_centroid.ra']
    if ra_column is not None:
        colRA_list = [ra_column]
    for colname in colRA_list:
        if colname in table.colnames:
            colRA = colname
            break
        elif colname.find('.')>=0 and colname.split('.')[0] in table.colnames:
            colRA = colname.split('.')
            break
    if colRA is None:
        print('Error! Could not find RA column ({}) in the input table ({})'.format(colRA_list, table.colnames))
        sys.exit(255)
    
    colDec = None
    colDec_list = ['DEC', 'Dec', 'DELTA_J2000', 'dec', 'sky_centroid.dec']
    if dec_column is not None:
        colDec_list = [dec_column]
    for colname in colDec_list:
        if colname in table.colnames:
            colDec = colname
            break
        elif colname.find('.')>=0 and colname.split('.')[0] in table.colnames:
            colDec = colname.split('.')
            break
    if colDec is None:
        print('Error! Could not find DEC column ({}) in the input table ({})'.format(colDec_list, table.colnames))
        sys.exit(255)
    
    
    # parsing input old and new RA Dec
    
    old_ra, old_dec = old_ra_dec
    new_ra, new_dec = new_ra_dec
    try:
        old_ra_dec_scoord = SkyCoord(float(old_ra)*u.deg, float(old_dec)*u.deg)
    except:
        old_ra_dec_scoord = SkyCoord(old_ra, old_dec, unit=(u.hour, u.deg))
    try:
        new_ra_dec_scoord = SkyCoord(float(new_ra)*u.deg, float(new_dec)*u.deg)
    except:
        new_ra_dec_scoord = SkyCoord(new_ra, new_dec, unit=(u.hour, u.deg))
    
    d_ra = new_ra_dec_scoord.ra.deg - old_ra_dec_scoord.ra.deg
    d_dec = new_ra_dec_scoord.dec.deg - old_ra_dec_scoord.dec.deg
    
    
    # update RA Dec in the table
    
    if isinstance(colRA, (list, tuple)):
        if isinstance(table[colRA[0]], SkyCoord):
            scoord = table[colRA[0]]
            table[colRA[0]] = SkyCoord(scoord.ra + (d_ra * u.deg), scoord.dec)
        else:
            dataRA = getattr(table[colRA[0]], colRA[1])
            if isinstance(dataRA, u.Quantity):
                setattr(table[colRA[0]], colRA[1], dataRA + (d_ra * u.deg))
            else:
                setattr(table[colRA[0]], colRA[1], dataRA + (d_ra))
    else:
        table[colRA] += d_ra
    
    if isinstance(colDec, (list, tuple)):
        if isinstance(table[colDec[0]], SkyCoord):
            scoord = table[colDec[0]]
            table[colDec[0]] = SkyCoord(scoord.ra, scoord.dec + (d_dec * u.deg))
        else:
            dataDec = getattr(table[colDec[0]], colDec[1])
            if isinstance(dataDec, u.Quantity):
                setattr(table[colDec[0]], colDec[1], dataDec + (d_dec * u.deg))
            else:
                setattr(table[colDec[0]], colDec[1], dataDec + (d_dec))
    else:
        table[colDec] += d_dec
    
    

    if os.path.isfile(output_catalog_file):
        print('Found existing output file, backing up as %r'%(output_catalog_file+'.backup'))
        shutil.move(output_catalog_file, output_catalog_file+'.backup')
    if output_catalog_file.find(os.sep)>=0:
        if not os.path.isdir(os.path.dirname(output_catalog_file)):
            os.makedirs(os.path.dirname(output_catalog_file))
    print('Writing catalog %r'%(output_catalog_file))
    if re.match(r'.*\.(txt|dat|lis)$', input_catalog_file): 
        table.write(output_catalog_file, format='ascii.fixed_width', delimiter=' ', bookend=True)
        with open(output_catalog_file, 'r+') as fp:
            fp.seek(0)
            fp.write('#')
    elif re.match(r'.*\.ecsv$', input_catalog_file):
        table.write(output_catalog_file, format='ascii.ecsv')
    else:
        table.write(output_catalog_file)



########
# MAIN #
########

if __name__ == '__main__':
    
    main()



