#!/usr/bin/env python
#
# Extract a source by SLTNAME, or RA DEC from stage 2 NIRCam WFSS cal.fits files.
# 
import os, sys, re, json, copy, glob, shutil
import click
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from click.exceptions import UsageError
from pprint import pprint
from tqdm import tqdm
from jwst import datamodels as dm
from stdatamodels.jwst.datamodels.multislit import MultiSlitModel
from stdatamodels.jwst.datamodels.slit import SlitModel, SlitDataModel


def radec2str(radec):
    return radec.to_string('hmsdms', precision=3).replace(' ', '_')


@click.command()
@click.argument('cal_files', type=str, nargs=-1, required=True)
@click.option('--slit-name', '--sltname', 'sltname', type=str, default=None)
@click.option('--ra-dec', 'radec', nargs=2, type=str, default=None)
@click.option('--matching-radius', 'matching_radius', type=float, default=2.0, help='Matching radius in arcsec, default is 2 arcsec.')
@click.option('--output-dir', 'output_dir', type=click.Path(), default='extract_wfss_slit', help='Output directory, default is "extract_wfss_slit".')
def main(cal_files, sltname, radec, matching_radius, output_dir):
    
    if sltname is None and radec is None:
        raise UsageError('Please input either --slit-name or --ra-dec.')
    elif radec is not None:
        ra_str, dec_str = radec
        if ra_str.find(':')>=0 or ra_str.find('m')>=0:
            radec = SkyCoord(ra_str, dec_str, unit=('hour', 'deg'))
        else:
            radec = SkyCoord(float(ra_str), float(dec_str), unit=('deg', 'deg'))
    
    if sltname is not None:
        print('Input slit name: {}'.format(sltname))
    else:
        print('Input RA Dec: {}, matching radius {} arcsec'.format(radec2str(radec), matching_radius))
    
    print('Input cal file(s): {}  (type {})'.format(cal_files, type(cal_files)))
    #raise NotImplementedError()

    cal_files = list(sorted(cal_files))
    print('We are going to search {} cal files'.format(len(cal_files)))

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
                    if 'SLTNAME' not in header:
                        continue
                    sltname2 = header['SLTNAME']
                    ra = header['SRCRA']
                    dec = header['SRCDEC']
                    radec2 = SkyCoord(ra, dec, unit=('deg', 'deg'))
                    if sltname is not None:
                        if sltname2 != sltname:
                            continue
                    else:
                        sep = radec2.separation(radec).to('arcsec').value
                        # if sltname2 == '672': # DEBUG
                        #     print(f'Checking sltname 672, sep {sep}') # DEBUG
                        if sep > matching_radius:
                            continue
                    # 
                    idx = int(iext/7)
                    basename = os.path.basename(cal_file).replace('.fits', '')
                    savename = 'slit_{}_ra_dec_{}_{}'.format(sltname2, radec2str(radec2), basename)
                    # 
                    if cal_file not in matched_dict:
                        matched_dict[cal_file] = []
                    matched_dict[cal_file].append(
                        (idx, savename)
                    )
                    # 
                    found_matches.append({
                        'sltname': sltname, 
                        'ra': ra, 
                        'dec': dec, 
                        'iext': iext, 
                        'cal_file': cal_file, 
                        'ical': ical, 
                    })
                    #pprint(found_matches[-1])
    
    if len(found_matches) == 0:
        print('No matched source is found!')
    else:
        print('We will output extracted slit data to {!r}'.format(output_dir))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # use jwst.datamodels to extract the slit data
        ical = 0
        ncals = len(matched_dict)
        for cal_file in matched_dict:
            print('Extracting wfss slit(s) from {!r} ({}/{})'.format(cal_file, ical+1, ncals))
            with dm.open(cal_file) as model:
                #indicies = np.argwhere([slit.name in slt_names for slit in model.slits]).ravel()
                for idx, savename in matched_dict[cal_file]:
                    #model.slits[idx]
                    # see - site-packages/stdatamodels/properties.py:475-476
                    # see - site-packages/stdatamodels/jwst/datamodels/multislit.py:36
                    slit = model.slits[idx]
                    slit2 = SlitModel(**dict(slit.items()))
                    slit2.meta = copy.deepcopy(slit.meta) # DZLIU: copy meta.wcs and meta.wcsinfo and meta.coordinates, see /home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1477/lib/python3.12/site-packages/jwst/extract_2d/grisms.py:531
                    #slit2.update(slit) # Error -- logs = d.pop("cal_logs", {}) -- AttributeError: No attribute 'pop'
                    model2 = MultiSlitModel(slit2)
                    model2.update(model)
                    # 
                    # DZLIU: copy WCS v3yangle cd matrix etc. see /home/dzliu/Software/CONDA/miniconda3/envs/jwstpmap1477/lib/python3.12/site-packages/jwst/extract_2d/grisms.py:663
                    for kk in model.meta.wcsinfo:
                        if not hasattr(model2.meta.wcsinfo, kk):
                            setattr(model2.meta.wcsinfo, kk, getattr(model.meta.wcsinfo, kk))
                    # 
                    output_name = savename + '.fits'
                    output_file = os.path.join(output_dir, output_name)
                    if os.path.exists(output_file):
                        shutil.move(output_file, output_file+'.backup')
                    model2.save(output_file)
            ical += 1



if __name__ == '__main__':
    main()


