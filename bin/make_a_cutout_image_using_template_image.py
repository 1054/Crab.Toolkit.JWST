#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A Python command line script by Daizhong Liu to reproject a FITS image/cube 
to the WCS of a template FITS image cube. 
"""
# 
# see also "alma_phangs_fits_cube_reproject_to_input_ra_dec_fov_pixsc.py"
# 
from __future__ import print_function
import os, sys, re, shutil, copy, time, datetime
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
from astropy.wcs import WCS
#from spectral_cube import SpectralCube, VaryingResolutionSpectralCube
#from radio_beam import Beam
from reproject import reproject_interp
import warnings
warnings.filterwarnings('ignore')


def usage():
    print('Usage: ')
    print('  make_a_cutout_image_using_reproject.py \\')
    print('    INPUT_FITS_IMAGE.fits \\')
    print('    TEMPLATE_FITS_IMAGE.fits \\')
    print('    OUTPUT_FITS_IMAGE.fits')
    print('')



def write_fits_cube_to_file(data, header, output_file, set_output_bitpix_32=False):
    #print('Writing to "%s"'%(output_file))
    # 
    if set_output_bitpix_32:
        header['BITPIX'] = -32
        data = data.astype(np.float32)
    # 
    hdu = fits.PrimaryHDU(data=data, header=header)
    # 
    if output_file.find(os.sep)>=0:
        if not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
    if os.path.isfile(output_file):
        shutil.move(output_file, output_file+'.backup')
    hdu.writeto(output_file)
    print('Output to "%s"'%(output_file))




# 
# main
# 
if __name__ == '__main__':
    
    # Check user input
    input_file = ''
    template_file = ''
    output_file = ''
    scale_by_error_propagation = False
    iarg = 1
    #argmode = ''
    while iarg <= len(sys.argv)-1:
        arg_str = sys.argv[iarg]
        if re.match(r'^[-]+[a-zA-Z].*', arg_str):
            arg_str = re.sub(r'^[-]+', '-', arg_str).lower()
            if arg_str == '-image' or arg_str == '-cube' or arg_str == '-input':
                if iarg+1 <= len(sys.argv)-1:
                    iarg += 1
                    input_file = sys.argv[iarg]
                    print('input_file = %r'%(input_file))
            elif arg_str == '-template':
                if iarg+1 <= len(sys.argv)-1:
                    iarg += 1
                    template_file = sys.argv[iarg]
                    print('template_file = %r'%(template_file))
            elif arg_str == '-out' or arg_str == '-output':
                if iarg+1 <= len(sys.argv)-1:
                    iarg += 1
                    output_file = sys.argv[iarg]
                    print('output_file = %r'%(output_file))
        else:
            if input_file == '':
                input_file = sys.argv[iarg]
                print('input_file = %r'%(input_file))
            elif template_file == '':
                template_file = sys.argv[iarg]
                print('template_file = %r'%(template_file))
            elif output_file == '':
                output_file = sys.argv[iarg]
                print('output_file = %r'%(output_file))
        # 
        iarg += 1
    # 
    # check user inputs
    if input_file == '' or template_file == '' or output_file == '':
        usage()
        sys.exit(255)
    
    # 
    # read fits
    dataI, headerI = None, None
    with fits.open(input_file) as hdul:
        if 'SCI' in hdul:
            dataI, headerI = hdul['SCI'].data, hdul['SCI'].header
        else:
            hdu0 = hdul[0]
            for hdu in hdul:
                if hdu.header['NAXIS'] >= 2:
                    dataI, headerI = hdu.data, hdu.header
                    break
    #dataI, headerI = fits.getdata(input_file, header=True)
    if dataI is None:
        raise Exception('Error! Could not read an image from {!r}'.format(input_file))
    wcsI = WCS(headerI, naxis=2)
    #print('wcsI.wcs.obsgeo', wcsI.wcs.obsgeo)
    wcsI.wcs.obsgeo[:] = np.nan
    #print('wcsI.wcs.obsgeo', wcsI.wcs.obsgeo)
    # this deals with the error:
    # ".* instance has no attribute 'in_observer_velocity_frame'"
    #print('wcsI.celestial', wcsI.celestial)
    #print('wcsI.has_celestial', wcsI.has_celestial)
    # "ValueError: Input WCS has celestial components but output WCS does not"
    
    # 
    # prepare template
    dataT, headerT = None, None
    with fits.open(template_file) as hdul:
        if 'SCI' in hdul:
            dataT, headerT = hdul['SCI'].data, hdul['SCI'].header
        else:
            hdu0 = hdul[0]
            for hdu in hdul:
                if hdu.header['NAXIS'] >= 2:
                    dataT, headerT = hdu.data, hdu.header
                    break
    #dataT, headerT = fits.getdata(template_file, header=True)
    if dataT is None:
        raise Exception('Error! Could not read an image from {!r}'.format(template_file))
    wcsT = WCS(headerT, naxis=2)
    #print('wcsT.wcs.obsgeo', wcsT.wcs.obsgeo)
    wcsT.wcs.obsgeo[:] = np.nan
    #print('wcsT.wcs.obsgeo', wcsT.wcs.obsgeo)
    
    # 
    # reshape cubeI, as our internal calculation is always 2D plus a channel axis
    shapeI = np.array(list(dataI.shape))
    if len(shapeI) < 2:
        raise Exception('Input data has a wrong dimension %d!'%(len(dataI.shape)))
    elif len(shapeI) == 2:
        nchanI = 1
    else:
        nchanI = np.prod(shapeI[0:-2])
    cubeI = dataI.reshape((nchanI, shapeI[-2], shapeI[-1]))
    
    # 
    # reshape cubeT, we will only use the 2D information in the template data
    shapeT = np.array(list(dataT.shape))
    if len(shapeT) < 2:
        raise Exception('Input data has a wrong dimension %d!'%(len(dataT.shape)))
    elif len(shapeT) == 2:
        nchanT = 1
    else:
        nchanT = np.prod(shapeT[0:-2])
    cubeT = dataT.reshape((nchanT, shapeT[-2], shapeT[-1]))
    
    # 
    # prepare output array
    shapeO = np.concatenate([shapeI[0:-2], [shapeT[-2], shapeT[-1]]])
    cubeO = np.full((nchanI, shapeT[-2], shapeT[-1]), fill_value=0.0)
    
    # 
    # prepare output header
    headerO = copy.deepcopy(headerI)
    wcsheaderT = wcsT.to_header()
    wcsheaderI = wcsI.to_header()
    #print('wcsheaderT', wcsheaderT)
    #print('wcsheaderI', wcsheaderI)
    for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
        if not (key in wcsheaderT):
            if key in headerO:
                del headerO[key]
    if 'CD3_3' in headerO and ('CD3_3' not in wcsheaderT):
        headerO['CDELT3'] = headerO['CD3_3']
        del headerO['CD3_3'] # reproject 
    for key in wcsheaderI:
        if not (key in wcsheaderT):
            if key in headerO:
                del headerO[key]
    for key in wcsheaderT:
        headerO[key] = wcsheaderT[key]
    headerO['HISTORY'] = ''
    headerO['HISTORY'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ' + time.strftime('%Z')
    headerO['HISTORY'] = 'Reprojected the WCS of the input FITS file "%s" to that of the template FITS file "%s".'%(input_file, template_file)
    headerO['HISTORY'] = ''
    
    # 
    # reproject
    for ichan in range(nchanI):
        print('Reprojecting %d/%d plate'%(ichan+1, nchanI))
        cubeO[ichan, :, :] = reproject_interp((cubeI[ichan, :, :], wcsI), wcsT, shape_out=(shapeT[-2], shapeT[-1]), return_footprint=False)
        # 
        #if debug:
        #    if ichan % int(np.ceil(nchanI/15)) == 0:
        #        write_fits_cube_to_file(data=cubeO, header=headerO, output_file=output_file.replace('.fits', '.dump.fits'))
    
    # 
    # reshape
    cubeO = cubeO.reshape(shapeO)
    
    # 
    # cleanup fits header
    if 'WCSAXES' in headerO:
        del headerO['WCSAXES']
    
    # 
    # output
    write_fits_cube_to_file(data=cubeO, header=headerO, output_file=output_file)



