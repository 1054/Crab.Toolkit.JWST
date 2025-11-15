#!/usr/bin/env python
#
# Find SLTNAME 
# 
import os, sys, re, json, copy, glob, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from pprint import pprint
from tqdm import tqdm


@click.command()
@click.argument('sltnames', type=str, nargs=-1)
@click.option('--output-dir', type=click.Path(), default=None, help='If specify an output directory, then we extract the fits extensions.')
def main(sltnames, output_dir):
    print('Input slit name(s): {}  (type {})'.format(sltnames, type(sltnames)))
    cal_files = glob.glob('jw*/calibrated2_cals/jw*_cal.fits')
    cal_files = list(sorted(cal_files))
    print('We are going to check {} cal_files'.format(len(cal_files)))

    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    found_matches = []
    for ical in tqdm(range(len(cal_files)), leave=False):
        cal_file = cal_files[ical]
        header0 = None
        with fits.open(cal_file) as hdul:
            for iext in tqdm(range(1, len(hdul), 7), desc='Cal file: {}'.format(cal_file)):
                if iext == 0:
                    header0 = hdul[0].header
                    scatfile = header0['SCATFILE']
                else:
                    hdu = hdul[iext]
                    header = hdu.header
                    if 'SLTNAME' in header and header['SLTNAME'] in sltnames:
                        sltname = header['SLTNAME']
                        ra = header['SRCRA']
                        dec = header['SRCDEC']
                        found_matches.append({
                            'sltname': sltname, 
                            'ra': ra, 
                            'dec': dec, 
                            'iext': iext, 
                            'cal_file': cal_file, 
                            'ical': ical, 
                        })
                        pprint(found_matches[-1])
                        if output_dir is not None:
                            output_hdul = fits.HDUList([
                                hdul[0], hdul[iext], 
                                hdul[iext+1], hdul[iext+2], hdul[iext+3], 
                                hdul[iext+4], hdul[iext+5], hdul[iext+6]
                            ])
                            #radecstr = '{:.7f}_{:.7f}'.format(ra, dec)
                            radecstr = SkyCoord(ra, dec, unit=('deg', 'deg')).to_string('hmsdms', precision=3)
                            radecstr = radecstr.replace(' ', '_')
                            cal_name = os.path.basename(cal_file).replace('.fits', '')
                            output_name = 'slit_{}_ra_dec_{}_{}.fits'.format(sltname, radecstr, cal_name)
                            output_file = os.path.join(output_dir, output_name)
                            output_hdul.writeto(output_file, overwrite=True)
                            print('Output to {!r}'.format(output_file))

    print('------')
    pprint(found_matches)


if __name__ == '__main__':
    main()


