#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.table import Table

# code name and version
CODE_NAME = 'util_convert_bintable_to_fits_image.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230626'
CODE_HOMEPAGE = 'https://github.com/1054/Crab.Toolkit.JWST'

# logging
#import logging
#logging.basicConfig(level=logging.INFO)
#logger = logging.getLogger(CODE_NAME)
#logger.setLevel(logging.DEBUG)



@click.command()
@click.argument('input_bintable_file', type=click.Path(exists=True))
@click.argument('output_fits_file', type=click.Path(exists=False))
def main(
        input_bintable_file, 
        output_fits_file, 
    ):
    
    table = None
    header = fits.Header()
    header['NAXIS'] = 2
    header['NAXIS1'] = 0
    header['NAXIS2'] = 0
    with fits.open(input_bintable_file) as hdul:
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU):
                table = Table(hdu.data)
                for key in hdu.header:
                    if key not in ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'PCOUNT', 'GCOUNT']:
                        header[key] = hdu.header[key]
                if len(table.colnames) == 1:
                    data = table[table.colnames[0]]
                    while len(data.shape) > 2 and data.shape[0] == 1:
                        data = data[0]
    
    if table is None:
        raise Exception('No BinTable is found in the input file {!r}'.format(input_bintable_file))
    
    if os.path.isfile(output_fits_file):
        shutil.move(output_fits_file, output_fits_file+'.backup')
    fits.PrimaryHDU(data=data).writeto(output_fits_file)
    
    print('Output to {!r}'.format(output_fits_file))



if __name__ == '__main__':
    
    main()



