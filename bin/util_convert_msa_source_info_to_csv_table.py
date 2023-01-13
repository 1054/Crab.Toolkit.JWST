#!/usr/bin/env python
# 

import os, sys, re, copy, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.table import Table

# code name and version
CODE_NAME = 'util_convert_msa_source_info_to_csv_table.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20230113'
CODE_HOMEPAGE = 'https://github.com/1054/Crab.Toolkit.JWST'

# logging
#import logging
#logging.basicConfig(level=logging.INFO)
#logger = logging.getLogger(CODE_NAME)
#logger.setLevel(logging.DEBUG)



@click.command()
@click.argument('input_msa_meta_file', type=click.Path(exists=True))
@click.argument('output_source_info_table', type=click.Path(exists=False))
def main(
        input_msa_meta_file, 
        output_source_info_table, 
    ):
    
    with fits.open(input_msa_meta_file) as hdul:
        table = Table(hdul['SOURCE_INFO'].data)
        if os.path.isfile(output_source_info_table):
            shutil.move(output_source_info_table, output_source_info_table+'.backup')
        if re.match(r'^.*\.(txt|dat)$', output_source_info_table):
            for colname in ['alias', 'source_name']:
                if np.any([t.find(' ')>=0 for t in table[colname]]):
                    table[colname] = ['"{}"'.format(t) for t in table[colname]] # add quotes
            table.write(output_source_info_table, format='ascii.fixed_width', delimiter=' ', bookend=True)
            with open(output_source_info_table, 'r+') as fp:
                fp.seek(0)
                fp.write('#')
        else:
            table.write(output_source_info_table)
    
    print('Output to {!r}'.format(output_source_info_table))



if __name__ == '__main__':
    
    main()



