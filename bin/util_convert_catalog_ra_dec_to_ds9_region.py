#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.table import Table
from matplotlib import cm
from matplotlib import colors as mpl_colors
from tqdm import tqdm

# code name and version
CODE_NAME = 'util_convert_catalog_ra_dec_to_ds9_region.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230108'
CODE_HOMEPAGE = ''

# logging
#import logging
#logging.basicConfig(level=logging.INFO)
#logger = logging.getLogger(CODE_NAME)
#logger.setLevel(logging.DEBUG)



@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('output_region_file', type=click.Path(exists=False))
@click.option('--color', type=str, default='green')
@click.option('--radius', type=float, default=0.5)
@click.option('--ra-column', type=str, default=None)
@click.option('--dec-column', type=str, default=None)
@click.option('--redshift-column', type=str, default=None)
@click.option('--id-column', type=str, default=None)
@click.option('--color-by-redshift', is_flag=True, default=False)
def main(
        input_catalog_file, 
        output_region_file, 
        color, 
        radius, 
        ra_column,
        dec_column,
        redshift_column,
        id_column, 
        color_by_redshift, 
    ):
    
    print('Reading catalog file {!r}'.format(input_catalog_file))
    if re.match(r'.*\.(txt|dat|lis)$', input_catalog_file): 
        table = Table.read(input_catalog_file, format='ascii')
    else:
        table = Table.read(input_catalog_file)
    
    print('Reading columns ...')
    
    colRA = None
    colRA_list = ['RA', 'ALPHA_J2000', 'ra']
    if ra_column is not None:
        colRA_list = [ra_column]
    for colname in colRA_list:
        if colname in table.colnames:
            colRA = colname
            break
    if colRA is None:
        print('Error! Could not find RA column ({}) in the input table ({})'.format(colRA_list, table.colnames))
        sys.exit(255)
    
    colDec = None
    colDec_list = ['DEC', 'Dec', 'DELTA_J2000', 'dec']
    if dec_column is not None:
        colDec_list = [dec_column]
    for colname in colDec_list:
        if colname in table.colnames:
            colDec = colname
            break
    if colDec is None:
        print('Error! Could not find DEC column ({}) in the input table ({})'.format(colDec_list, table.colnames))
        sys.exit(255)
    
    colRedshift = None
    if color_by_redshift:
        colRedshift_list = ['z', 'PHOTOZ', 'z_phot', 'ez_z_phot']
        if redshift_column is not None:
            colRedshift_list = [redshift_column]
        for colname in colRedshift_list:
            if colname in table.colnames:
                colRedshift = colname
                break
        if colRedshift is None:
            print('Error! Could not find Redshift column ({}) in the input table ({})'.format(colRedshift_list, table.colnames))
            sys.exit(255)
        color_mapper = cm.ScalarMappable(norm=mpl_colors.Normalize(vmin=-1, vmax=8), cmap=cm.gist_rainbow)
    
    colID = None
    colID_list = ['ID', 'id']
    if id_column is not None:
        colID_list = [id_column]
    for colname in colID_list:
        if colname in table.colnames:
            colID = colname
            break
    #if color_by_ID and colID is None:
    #    print('Error! Could not find ID column ({}) in the input table ({})'.format(colID_list, table.colnames))
    #    sys.exit(255)
    
    print('Writing region file {!r} ...'.format(output_region_file))
    
    if os.path.isfile(output_region_file):
        shutil.move(output_region_file, output_region_file+'.backup')
    
    with open(output_region_file, 'w') as fp:
        fp.write('# DS9 region file\n')
        fp.write('global color={}\n'.format(color))
        fp.write('fk5\n')
        for i in tqdm(range(len(table))):
            line_str = 'circle({},{},{}")'.format(table[colRA][i], table[colDec][i], radius)
            text_str = ''
            color_str = ''
            if colID is not None:
                text_str += '{}'.format(table[colID][i])
            if colRedshift is not None:
                if text_str != '':
                    text_str += ','
                redshift_value = table[colRedshift][i]
                if np.isfinite(redshift_value) and redshift_value>0.0:
                    text_str += 'z={:.3f}'.format(redshift_value)
                    if color_by_redshift:
                        color_str = mpl_colors.to_hex(color_mapper.to_rgba(redshift_value))
            if text_str != '' and color_str != '': 
                line_str += ' # text={{{}}} color={{{}}}'.format(text_str, color_str)
            elif text_str != '': 
                line_str += ' # text={{{}}}'.format(text_str)
            elif color_str != '': 
                line_str += ' # color={{{}}}'.format(color_str)
            fp.write(line_str+'\n')
    
    print('Output to {!r}'.format(output_region_file))



if __name__ == '__main__':
    
    main()



