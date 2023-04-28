#!/usr/bin/env python
# 
import os, sys, re, shutil, time, datetime, copy
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS
from astropy.io import fits # fits.ImageHDU
#from astropy.modeling import rotations # rotations.Rotation2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import reproject
from reproject import reproject_interp


def usage():
    print('Usage: ')
    print('    make_a_cutout_image_using_reproject.py input_image.fits center_RA center_Dec FoV_RA FoV_Dec output_image.fits')
    print('Notes: ')
    print('    center_RA and center_Dec must be in units of degrees.')
    print('    FoV_RA and FoV_Dec must be in units of arcsec.')
    print('')



def parse_angle_str(angle_str, default_unit='degree', output_unit='arcsec'):
    angle_arcsec = np.nan
    angle_value = np.nan
    angle_unit = ''
    if re.match(r'^([0-9eE.+-]+)[ \t]*$', angle_str, re.IGNORECASE):
        angle_value = float(re.sub(r'^([0-9eE.+-]+)[ \t]*$', r'\1', angle_str, re.IGNORECASE))
        angle_unit = default_unit
    elif re.match(r'^([0-9eE.+-]+)[ \t]*arcsec$', angle_str, re.IGNORECASE):
        angle_value = float(re.sub(r'^([0-9eE.+-]+)[ \t]*arcsec$', r'\1', angle_str, re.IGNORECASE))
        angle_unit = 'arcsec'
    elif re.match(r'^([0-9eE.+-]+)[ \t]*arcmin$', angle_str, re.IGNORECASE):
        angle_value = float(re.sub(r'^([0-9eE.+-]+)[ \t]*arcmin$', r'\1', angle_str, re.IGNORECASE))
        angle_unit = 'arcmin'
    elif re.match(r'^([0-9eE.+-]+)[ \t]*(deg|degree)$', angle_str, re.IGNORECASE):
        angle_value = float(re.sub(r'^([0-9eE.+-]+)[ \t]*(deg|degree)', r'\1', angle_str, re.IGNORECASE))
        angle_unit = 'degree'
    else:
        raise Exception('Error! Could not parse the input angle_str=%r when calling parse_angle_str(). It should either be a float number or a float number with a unit \'arcsec\', \'arcmin\' or \'degree\'.'%(angle_str))
    # 
    if angle_unit == 'arcsec':
        angle_arcsec = angle_value
    elif angle_unit == 'arcmin':
        angle_arcsec = angle_value * 60.
    elif angle_unit == 'degree':
        angle_arcsec = angle_value * 3600.
    else:
        raise Exception('Error! angle_unit in parse_angle_str() is not one of \'arcsec\', \'arcmin\' or \'degree\'. This should not happen.')
    # 
    if output_unit == 'arcsec':
        return angle_arcsec
    elif output_unit == 'arcsec':
        return angle_arcsec/60.
    elif output_unit == 'deg' or output_unit == 'degree':
        return angle_arcsec/3600.
    else:
        raise Exception('Error! Wrong input output_unit=%r when calling parse_angle_str(). It must be either \'arcsec\', \'arcmin\' or \'degree\'.'%(output_unit))



# read user input
input_image_file = None
center_RA = None
center_Dec = None
FoV_RA = None
FoV_Dec = None
pixel_size = None
output_image_file = None
set_output_bitpix_32 = False
iarg = 1
arg_str = ''
while iarg < len(sys.argv):
    arg_str = sys.argv[iarg].lower()
    arg_str = re.sub(r'^[-]+', '-', arg_str)
    if arg_str == '-set-output-bitpix-32':
        set_output_bitpix_32 = True
        print('Setting set_output_bitpix_32 = %s'%(set_output_bitpix_32))
    elif arg_str == '-input' or arg_str == '-in':
        if iarg+1 < len(sys.argv):
            iarg += 1
            input_image_file = sys.argv[iarg]
            print('Setting input_image_file = %r'%(input_image_file))
    elif arg_str == '-output' or arg_str == '-out':
        if iarg+1 < len(sys.argv):
            iarg += 1
            output_image_file = sys.argv[iarg]
            print('Setting output_image_file = %r'%(output_image_file))
    elif arg_str == '-pixel-size':
        if iarg+1 < len(sys.argv):
            iarg += 1
            pixel_size = float(sys.argv[iarg])
            print('Setting pixel_size = %s [arcsec]'%(pixel_size))
    elif arg_str == '-field-of-view' or arg_str == '-fov':
        if iarg+2 < len(sys.argv):
            iarg += 1
            FoV_RA = float(sys.argv[iarg])
            print('Setting FoV_RA = %s [arcsec]'%(FoV_RA))
            iarg += 1
            FoV_Dec = float(sys.argv[iarg])
            print('Setting FoV_Dec = %s [arcsec]'%(FoV_Dec))
    elif arg_str == '-center' or arg_str == '-cen':
        if iarg+2 < len(sys.argv):
            iarg += 1
            try:
                center_RA = float(sys.argv[iarg])
            except:
                center_RA = str(sys.argv[iarg])
            print('Setting center_RA = %s [degree]'%(center_RA))
            iarg += 1
            try:
                center_Dec = float(sys.argv[iarg])
            except:
                center_Dec = str(sys.argv[iarg])
            print('Setting center_Dec = %s [degree]'%(center_Dec))
    else:
        if input_image_file is None:
            input_image_file = sys.argv[iarg]
            print('Setting input_image_file = %s'%(input_image_file))
        elif center_RA is None:
            try:
                center_RA = float(sys.argv[iarg])
            except:
                center_RA = str(sys.argv[iarg])
            print('Setting center_RA = %s'%(center_RA))
        elif center_Dec is None:
            try:
                center_Dec = float(sys.argv[iarg])
            except:
                center_Dec = str(sys.argv[iarg])
            print('Setting center_RA = %s'%(center_Dec))
        elif FoV_RA is None:
            FoV_RA = parse_angle_str(sys.argv[iarg], default_unit='arcsec', output_unit='arcsec')
            print('Setting FoV_RA = %s'%(FoV_RA))
        elif FoV_Dec is None:
            FoV_Dec = np.abs(parse_angle_str(sys.argv[iarg], default_unit='arcsec', output_unit='arcsec'))
            print('Setting FoV_Dec = %s'%(FoV_Dec))
        elif output_image_file is None:
            output_image_file = sys.argv[iarg]
            print('Setting output_image_file = %r'%(output_image_file))
    iarg += 1


# check user input
if input_image_file is None or center_RA is None or center_Dec is None or FoV_RA is None or FoV_Dec is None or output_image_file is None:
    usage()
    sys.exit()


# read input fits image
regex_match = re.match(r'^(.*.fits)\[(.+)\]$', input_image_file, re.IGNORECASE) # 20230428: input IRAF style fits with extension
if regex_match:
    input_image_file_path = regex_match.group(1)
    input_image_extension = regex_match.group(2)
else:
    input_image_file_path = input_image_file
    input_image_extension = None
with fits.open(input_image_file_path) as hdulist:
    ihdu = 0
    while ihdu < len(hdulist):
        #print(type(hdulist[ihdu]))
        if isinstance(hdulist[ihdu], fits.hdu.image.PrimaryHDU):
            main_header = copy.copy(hdulist[ihdu].header)
        if isinstance(hdulist[ihdu], fits.ImageHDU) or isinstance(hdulist[ihdu], fits.hdu.image.PrimaryHDU):
            if hdulist[ihdu].header['NAXIS'] >= 2:
                print('reading ihdu %d header'%(ihdu))
                header = copy.copy(hdulist[ihdu].header)
                break
        ihdu += 1
    if input_image_extension is None:
        ihdu = input_image_extension
    print('reading ihdu %d data'%(ihdu))
    image = copy.copy(hdulist[ihdu].data)
    for key in hdulist[ihdu].header:
        header[key] = hdulist[ihdu].header[key] # some cases extensions do not have wcs header, this solves that issue.


# determine wcs, pixscale, x_size y_size
wcs = WCS(header, naxis=2)
wcs.sip = None
if pixel_size is None:
    pixscales = proj_plane_pixel_scales(wcs) * 3600.0
    print('pixscales', pixscales)
else:
    pixscales = [-pixel_size, pixel_size]
x_pixsc = np.abs(pixscales[0])
y_pixsc = np.abs(pixscales[1])
x_sizeF = FoV_RA / x_pixsc
y_sizeF = FoV_Dec / y_pixsc
x_size = int(np.ceil(FoV_RA / x_pixsc))
y_size = int(np.ceil(FoV_Dec / y_pixsc))
print('x_size', x_size, 'y_size', y_size)


# convert RA Dec sexagesimal format
if isinstance(center_RA, str) or isinstance(center_Dec, str):
    center_skycoord = SkyCoord(center_RA, center_Dec, unit=(u.hour, u.deg))
    center_RA = center_skycoord.ra.deg
    center_Dec = center_skycoord.dec.deg


# prepare cutout header
#cutout_header = copy.deepcopy(header)
cutout_header = fits.Header()
cutout_header['BITPIX'] = -32
cutout_header['NAXIS'] = 2
cutout_header['NAXIS1'] = x_size
cutout_header['NAXIS2'] = y_size
cutout_header['CTYPE1'] = 'RA---TAN'
cutout_header['CTYPE2'] = 'DEC--TAN'
cutout_header['CUNIT1'] = 'deg'
cutout_header['CUNIT2'] = 'deg'
cutout_header['CDELT1'] = -x_pixsc / 3600.0
cutout_header['CDELT2'] = y_pixsc / 3600.0
cutout_header['CRPIX1'] = (x_size+1)/2. # 1-based number
cutout_header['CRPIX2'] = (y_size+1)/2. # 1-based number
cutout_header['CRVAL1'] = center_RA
cutout_header['CRVAL2'] = center_Dec
cutout_header['RADESYS'] = 'ICRS'
cutout_header['EQUINOX'] = 2000
#'NAXIS','NAXIS1','NAXIS2','CDELT1','CDELT2','CRPIX1','CRPIX2','CRVAL1','CRVAL2',
for key in ['BUNIT','BMAJ','BMIN','BPA','TELESCOP','INSTRUME','FILTER','EXPTIME',
            'DATE-OBS','TIME-OBS','PHOTMODE','PHOTFLAM','PHTFLAM1','PHTFLAM2','PHTRATIO','PHOTFNU','PHOTZPT','PHOTPLAM','PHOTBW']:
    if key in main_header:
        cutout_header[key] = main_header[key]
    if key in header:
        cutout_header[key] = header[key]
cutout_header['HISTORY'] = ''
cutout_header['HISTORY'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ' + time.strftime('%Z')
cutout_header['HISTORY'] = 'Created cutout image at center RA Dec %s %s with FoV size RA Dec %s %s [arcsec] from the input image file "%s"'%(\
                           center_RA, center_Dec, FoV_RA, FoV_Dec, input_image_file)
cutout_header['HISTORY'] = ''
#print(cutout_header)


cutout_image, cutout_footprint = reproject_interp((image, wcs), cutout_header)


# generate fits HDU
output_hdu = fits.PrimaryHDU(data=cutout_image, header=cutout_header)

if set_output_bitpix_32:
    output_hdu.header['BITPIX'] = -32
    output_hdu.data[np.isnan(output_hdu.data)] = 0.0
    output_hdu.data = output_hdu.data.astype(np.float32)


# prepare to write output file
output_file = output_image_file
if output_file.find(os.sep) >= 0:
    if not os.path.isdir(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
if os.path.isfile(output_file):
    shutil.move(output_file, output_file+'.backup')

output_hdu.writeto(output_file, overwrite=True)
print('Written to "%s"'%(output_file))



