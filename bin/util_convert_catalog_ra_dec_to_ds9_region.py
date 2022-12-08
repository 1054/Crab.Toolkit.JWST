#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
from astropy.table import Table

@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('output_region_file', type=click.Path(exists=False))
@click.option('--color', type=str, default='green')
@click.option('--radius', type=float, default=0.5)
def main(
        input_catalog_file, 
        output_region_file, 
        color, 
        radius, 
    ):
    
    table = Table.read(input_catalog_file)
    
    colRA = None
    colRA_list = ['RA', 'ALPHA_J2000', 'ra']
    for colname in colRA_list:
        if colname in table.colnames:
            colRA = colname
            break
    if colRA is None:
        print('Error! Could not find RA column ({}) in the input table ({})'.format(colRA_list, table.colnames))
        sys.exit(255)
    colDec = None
    colDec_list = ['DEC', 'Dec', 'DELTA_J2000', 'dec']
    for colname in colDec_list:
        if colname in table.colnames:
            colDec = colname
            break
    if colDec is None:
        print('Error! Could not find DEC column ({}) in the input table ({})'.format(colDec_list, table.colnames))
        sys.exit(255)
    
    with open(output_region_file, 'w') as fp:
        fp.write('# DS9 region file\n')
        fp.write('global color={}\n'.format(color))
        fp.write('fk5\n')
        for i in range(len(table)):
            fp.write('circle({},{},{}")\n'.format(
                table[colRA][i],
                table[colDec][i],
                radius,
            ))
    print('Output to {!r}'.format(output_region_file))



if __name__ == '__main__':
    
    main()



