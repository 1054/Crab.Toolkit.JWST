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
from jwst import datamodels


# class VariadicOption(click.Option):
#     def __init__(self, *args, **kwargs):
#         kwargs.setdefault("multiple", True)
#         super().__init__(*args, **kwargs)

#     def handle_parse_result(self, ctx, opts, args):
#         # Check if the flag was actually provided in the command line
#         if self.name in opts:
#             # 'args' contains the remaining tokens the parser didn't know what to do with
#             # We grab everything from 'args' until we hit another option (starting with -)
#             while args and not args[0].startswith('-'):
#                 # Note: Click's 'multiple' returns a tuple, so we extend it
#                 opts[self.name] = opts[self.name] + (args.pop(0),)
        
#         return super().handle_parse_result(ctx, opts, args)



@click.command()
@click.argument('input_slit_names_and_cal_files', type=str, nargs=-1, required=True)
#@click.argument('slt_names', type=str, nargs=-1, required=True)
#@click.option('--cal-files', 'input_cal_files', type=str, cls=VariadicOption) # nargs=-1
@click.option('--output-dir', '--output', 'output_dir', type=click.Path(), default=None, help='If specify an output directory, then we extract the fits extensions.')
def main(input_slit_names_and_cal_files, output_dir):
    slt_names = []
    cal_files = []
    for input_str in input_slit_names_and_cal_files:
        if input_str.endswith('.fits'):
            cal_files.append(input_str)
        elif input_str.find('*')>=0:
            cal_files.extend(glob.glob(input_str))
        else:
            slt_names.append(input_str)
    
    if len(cal_files) == 0:
        cal_files = glob.glob('jw*/calibrated2_cals/jw*_cal.fits')
    
    print('Input slit name(s): {}  (type {})'.format(slt_names, type(slt_names)))
    print('Input cal file(s): {}  (type {})'.format(cal_files, type(cal_files)))
    #raise NotImplementedError()

    if len(cal_files) == 0:
        raise Exception('Error! No cal file found in "jw*/calibrated2_cals/jw*_cal.fits" or provided via --cal-files.')
    cal_files = list(sorted(cal_files))
    print('We are going to check {} cal_files'.format(len(cal_files)))

    if output_dir is not None:
        print('We will output to {!r}'.format(output_dir))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    found_matches = []
    matched_dict = {}
    for ical in tqdm(range(len(cal_files)), leave=False):
        cal_file = cal_files[ical]
        header0 = None
        with fits.open(cal_file) as hdul:
            # cal file contains N*7 FITS extensions
            for iext in tqdm(range(1, len(hdul), 7), desc='Cal file: {}'.format(cal_file)):
                if iext == 0:
                    header0 = hdul[0].header
                    scatfile = header0['SCATFILE']
                else:
                    hdu = hdul[iext]
                    header = hdu.header
                    if 'SLTNAME' in header and header['SLTNAME'] in slt_names:
                        sltname = header['SLTNAME']
                        ra = header['SRCRA']
                        dec = header['SRCDEC']
                        radecstr = SkyCoord(ra, dec, unit=('deg', 'deg')).to_string('hmsdms', precision=3)
                        radecstr = radecstr.replace(' ', '_')
                        cal_name = os.path.basename(cal_file).replace('.fits', '')
                        idx = int(iext/7)
                        savename = 'slit_{}_ra_dec_{}_{}'.format(sltname, radecstr, cal_name)
                        
                        if cal_file not in matched_dict:
                            matched_dict[cal_file] = []
                        matched_dict[cal_file].append(
                            (idx, savename)
                        )

                        found_matches.append({
                            'sltname': sltname, 
                            'ra': ra, 
                            'dec': dec, 
                            'iext': iext, 
                            'cal_file': cal_file, 
                            'ical': ical, 
                        })
                        pprint(found_matches[-1])
                        # if output_dir is not None:
                        #     output_hdul = [
                        #         hdul[0], hdul[iext], 
                        #         hdul[iext+1], hdul[iext+2], hdul[iext+3], 
                        #         hdul[iext+4], hdul[iext+5], hdul[iext+6],
                        #     ]
                        #     if (hdul[-1].header['EXTNAME'] == 'ASDF'):
                        #         output_hdul.append(hdul[-1])
                        #     output_hdul = fits.HDUList(output_hdul)
                        #     #radecstr = '{:.7f}_{:.7f}'.format(ra, dec)
                        #     radecstr = SkyCoord(ra, dec, unit=('deg', 'deg')).to_string('hmsdms', precision=3)
                        #     radecstr = radecstr.replace(' ', '_')
                        #     cal_name = os.path.basename(cal_file).replace('.fits', '')
                        #     output_name = 'slit_{}_ra_dec_{}_{}.fits'.format(sltname, radecstr, cal_name)
                        #     output_file = os.path.join(output_dir, output_name)
                        #     output_hdul.writeto(output_file, overwrite=True)
                        #     print('Output to {!r}'.format(output_file))
        
    # use jwst.datamodels to extract the slit data
    if output_dir is not None:
        for cal_file in matched_dict:
            with datamodels.open(cal_file) as dm:
                #indicies = np.argwhere([slit.name in slt_names for slit in dm.slits]).ravel()
                for idx, savename in matched_dict[cal_file]:
                    #dm.slits[idx]
                    #see - site-packages/stdatamodels/properties.py:475-476
                    from stdatamodels.jwst.datamodels.multislit import MultiSlitModel
                    from stdatamodels.jwst.datamodels.slit import SlitDataModel
                    # see - site-packages/stdatamodels/jwst/datamodels/multislit.py:36
                    dm2 = MultiSlitModel(SlitModel(**dict(dm.slits[idx].items())))
                    output_name = savename + '.fits'
                    output_file = os.path.join(output_dir, output_name)
                    if os.path.exists(output_file):
                        shutil.move(output_file, output_file+'.backup')
                    dm2.save(output_file)

    print('------')
    pprint(found_matches)


if __name__ == '__main__':
    main()


