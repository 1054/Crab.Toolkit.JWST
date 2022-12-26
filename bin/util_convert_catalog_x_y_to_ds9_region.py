#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
from astropy.table import Table

@click.command()
@click.argument('input_catalog_file', type=click.Path(exists=True))
@click.argument('output_region_file', type=click.Path(exists=False))
@click.option('--color', type=str, default='green')
@click.option('--radius', type=float, default=5.0)
@click.option('--index-base', type=click.Choice(['0', '1']), default='0', help='x y index base, 0 or 1.')
def main(
        input_catalog_file, 
        output_region_file, 
        color, 
        radius, 
        index_base,
    ):
    
    table = Table.read(input_catalog_file)
    assert ('x' in table.colnames)
    assert ('y' in table.colnames)
    to_ds9_pixcoord = 0
    if int(index_base) == 0:
        to_ds9_pixcoord = 1
    with open(output_region_file, 'w') as fp:
        fp.write('# DS9 region file\n')
        fp.write('global color={}\n'.format(color))
        fp.write('image\n')
        for i in range(len(table)):
            fp.write('circle({},{},{})\n'.format(
                table['x'][i] + to_ds9_pixcoord,
                table['y'][i] + to_ds9_pixcoord,
                radius,
            ))
    print('Output to {!r}'.format(output_region_file))



if __name__ == '__main__':
    
    main()



