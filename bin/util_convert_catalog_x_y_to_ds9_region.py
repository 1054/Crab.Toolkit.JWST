#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.table import Table

@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('output_region_file', type=click.Path(exists=False))
@click.option('--x-column', type=str, default='x', help='The column name for x coordinate.')
@click.option('--y-column', type=str, default='y', help='The column name for y coordinate.')
@click.option('--radius-column', type=str, default='', help='The column name for radius in pixels.')
@click.option('--color', type=str, default='green')
@click.option('--radius', type=float, default=5.0)
@click.option('--index-base', type=click.Choice(['0', '1']), default='0', help='x y index base, 0 or 1. Deafult is 0.')
def main(
        input_catalog_file, 
        output_region_file, 
        x_column, 
        y_column, 
        radius_column, 
        color, 
        radius, 
        index_base,
    ):
    
    table = Table.read(input_catalog_file)
    assert (x_column in table.colnames)
    assert (y_column in table.colnames)
    x_array = table[x_column]
    y_array = table[y_column]
    if radius_column != '':
        radius_array = table[radius_column]
    else:
        radius_array = np.full(x_array.shape, fill_value=radius)
    to_ds9_pixcoord = 0
    if int(index_base) == 0:
        to_ds9_pixcoord = 1
    with open(output_region_file, 'w') as fp:
        fp.write('# DS9 region file\n')
        fp.write('global color={}\n'.format(color))
        fp.write('image\n')
        for i in range(len(table)):
            fp.write('circle({},{},{})\n'.format(
                x_array[i] + to_ds9_pixcoord,
                y_array[i] + to_ds9_pixcoord,
                radius_array[i],
            ))
    print('Output to {!r}'.format(output_region_file))



if __name__ == '__main__':
    
    main()



