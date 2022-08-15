#!/usr/bin/env python
# 
# 2021-03-31 copied from "trim_field_of_view_for_a_cube.py", made it more general for FITS image 2D 3D 4D
# 2021-07-10 added ext_id ext_name
# 
from __future__ import print_function
import os, sys, re, shutil, glob, copy, datetime, time
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
#from astropy.convolution import Gaussian2DKernel, convolve
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
#from spectral_cube import SpectralCube, VaryingResolutionSpectralCube
#from radio_beam import Beam
#from reproject import reproject_interp
#import scipy.stats
import warnings
warnings.filterwarnings('ignore')


def usage():
    print('Usage: ')
    print('  make_a_cutout_image_by_center_RA_Dec_FoV.py \\')
    print('    INPUT_FITS_IMAGE.fits \\')
    print('    OUTPUT_FITS_IMAGE.fits \\')
    print('    -fov "2arcsec" \\')
    print('    [-center RA Dec] \\')
    print('    [-template INPUT_TEMPLATE_IMAGE.fits] \\')
    print('    [-ext extension_number] \\')
    print('    [-extname extension_name] \\')
    print('    [-set-output-bitpix-32] \\')
    print('    [-set-output-naxis-2] ')
    print('Notes: ')
    print('  ')
    print('')



def parse_fov_str(fov_str):
    fov_arcsec = np.nan
    if re.match(r'^([0-9eE.+-]+)[ \t]*$', fov_str, re.IGNORECASE):
        fov_arcsec = float(re.sub(r'^([0-9eE.+-]+)[ \t]*$', r'\1', fov_str, re.IGNORECASE))
    elif re.match(r'^([0-9eE.+-]+)[ \t]*arcsec$', fov_str, re.IGNORECASE):
        fov_arcsec = float(re.sub(r'^([0-9eE.+-]+)[ \t]*arcsec$', r'\1', fov_str, re.IGNORECASE))
    elif re.match(r'^([0-9eE.+-]+)[ \t]*arcmin$', fov_str, re.IGNORECASE):
        fov_arcsec = float(re.sub(r'^([0-9eE.+-]+)[ \t]*arcmin$', r'\1', fov_str, re.IGNORECASE)) * 60.
    elif re.match(r'^([0-9eE.+-]+)[ \t]*(deg|degree)$', fov_str, re.IGNORECASE):
        fov_arcsec = float(re.sub(r'^([0-9eE.+-]+)[ \t]*(deg|degree)', r'\1', fov_str, re.IGNORECASE)) * 3600.
    return fov_arcsec



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
    elif angle_unit == 'deg' or angle_unit == 'degree':
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



def parse_skycoord(input_ra_str, input_dec_str, default_unit='deg', default_frame=FK5, return_skycoord_obj=False):
    # 
    input_ra_str = str(input_ra_str).strip()
    input_dec_str = str(input_dec_str).strip()
    # 
    output_RA_deg = np.nan
    output_Dec_deg = np.nan
    # 
    regex_skycoord_numeric = re.compile(r'^[0-9.+-]+$')
    regex_skycoord_deg = re.compile(r'^[0-9.+-]+(deg|degree)$')
    # 
    unit_ra = u.hour
    unit_dec = u.deg
    # 
    if regex_skycoord_numeric.match(input_ra_str) and regex_skycoord_numeric.match(input_dec_str):
        input_ra_str += default_unit
        input_dec_str += default_unit
        unit_ra = u.deg
    elif regex_skycoord_deg.match(input_ra_str) and regex_skycoord_deg.match(input_dec_str):
        unit_ra = u.deg
    # 
    skycoord_obj = SkyCoord(input_ra_str+' '+input_dec_str, unit=(unit_ra, unit_dec), frame=default_frame)
    # 
    output_RA_deg = skycoord_obj.ra.deg
    output_Dec_deg = skycoord_obj.dec.deg
    # 
    if return_skycoord_obj:
        return skycoord_obj
    # 
    return output_RA_deg, output_Dec_deg



def trim_header_dimension(header, naxis, verbose=False):
    if header['NAXIS'] > naxis:
        old_naxis = header['NAXIS']
        if verbose:
            print('Trimming header from %d to %d dimensions'%(old_naxis, naxis))
        header['NAXIS'] = naxis
        for i in range(naxis+1, old_naxis+1):
            for key in ['NAXIS', 'CDELT', 'CRVAL', 'CRPIX', 'CUNIT', 'CTYPE', 'CROTA']:
                key2 = key+'%d'%(i)
                if key2 in header:
                    if verbose:
                        print('del header[%r]'%(key2))
                    del header[key2]
            for key in ['CD', 'PC']:
                for j in range(1, old_naxis+1):
                    key2 = key+'%d_%d'%(i, j)
                    if key2 in header:
                        if verbose:
                            print('del header[%r]'%(key2))
                        del header[key2]
                    key2 = key+'%d_%d'%(j, i)
                    if key2 in header:
                        if verbose:
                            print('del header[%r]'%(key2))
                        del header[key2]
    return header



def trim_data_dimension(data, header, naxis, verbose=False, warning=True):
    if verbose:
        print('Trimming data from %d to %d dimensions'%(len(data.shape), naxis))
    old_data_shape = data.shape
    for i in range(naxis, len(old_data_shape)+1):
        if warning and old_data_shape[i-1] > 1:
            print('Warning! The %d-%d-th data along the original %d-th dimension will be lost!'%(2, old_data_shape[i-1]+1, i))
        data = data[0]
    header = trim_header_dimension(header, naxis, verbose=verbose)
    return data, header



def read_fits_data(input_file, naxis=None, ext_id=None, ext_name=None):
    """Read a fits data from a fits file.
    If the input fits file has more than three dimension, then we will 
    trim higher dimensions.
    """
    if ext_id is not None:
        print('Reading "%s" ext_id %d'%(input_file, ext_id))
        data, header = fits.getdata(input_file, ext_id, header=True)
    elif ext_name is not None:
        print('Reading "%s" ext_name "%s"'%(input_file, ext_name))
        data, header = fits.getdata(input_file, extname=ext_name, header=True)
    else:
        print('Reading "%s"'%(input_file))
        data, header = fits.getdata(input_file, header=True)
    if naxis is not None:
        if header['NAXIS'] > naxis:
            print('Trimming to %d dimensions'%(naxis))
            old_naxis = header['NAXIS']
            header['NAXIS'] = naxis
            for i in range(naxis+1, old_naxis+1):
                for key in ['NAXIS', 'CDELT', 'CRVAL', 'CRPIX', 'CUNIT', 'CTYPE']:
                    key2 = key+'%d'%(i)
                    if key2 in header:
                        del header[key2]
                for key in ['CD', 'PC', 'PV']:
                    for j in range(1, old_naxis+1):
                        key2 = key+'%d_%d'%(i, j)
                        if key2 in header:
                            del header[key2]
                        key2 = key+'%d_%d'%(j, i)
                        if key2 in header:
                            del header[key2]
                data = data[0]
    return data, header



def write_fits_data_to_file(data, header, output_file, set_output_bitpix_32=False):
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
    
    # Read user inputs
    input_file = ''
    output_file = ''
    template_file = ''
    debug = False
    new_center = None
    new_fov = []
    ext_id = None
    ext_name = None
    set_output_bitpix_32 = False
    set_output_naxis_2 = False
    set_output_naxis_3 = False
    iarg = 1
    while iarg <= len(sys.argv)-1:
        arg_str = sys.argv[iarg].lower()
        arg_str = re.sub(r'^[-]+', '-', arg_str)
        if arg_str == '-debug':
            debug = True
        elif arg_str == '-set-output-bitpix-32':
            set_output_bitpix_32 = True
            print('Setting set_output_bitpix_32 = %s'%(set_output_bitpix_32))
        elif arg_str == '-set-output-naxis-2':
            set_output_naxis_2 = True
            print('Setting set_output_naxis_2 = %s'%(set_output_naxis_2))
        elif arg_str == '-set-output-naxis-3':
            set_output_naxis_3 = True
            print('Setting set_output_naxis_3 = %s'%(set_output_naxis_3))
        elif arg_str == '-center' or arg_str == '-cen':
            if iarg+2 <= len(sys.argv)-1:
                iarg += 1
                input_RA_str = sys.argv[iarg]
                iarg += 1
                input_Dec_str = sys.argv[iarg]
                new_center = parse_skycoord(input_RA_str, input_Dec_str, return_skycoord_obj=True)
                print('Setting new_center = %s'%(new_center))
        elif arg_str == '-fov':
            if iarg+1 <= len(sys.argv)-1 and re.match(r'^([0-9eE.+-]+)[ \t]*(arcsec|arcmin|deg|degree|)$', sys.argv[iarg+1]):
                iarg += 1
                new_fov = [parse_fov_str(sys.argv[iarg])]
                if iarg+1 <= len(sys.argv)-1 and re.match(r'^([0-9eE.+-]+)[ \t]*(arcsec|arcmin|deg|degree|)$', sys.argv[iarg+1]):
                    iarg += 1
                    new_fov.append(parse_fov_str(sys.argv[iarg]))
                else:
                    new_fov.append(new_fov[0])
                print('Setting new_fov = %s [arcsec]'%(new_fov))
        elif arg_str == '-template':
            if iarg+1 <= len(sys.argv)-1 and re.match(r'^.*\.fits(|\.gz)$', sys.argv[iarg+1], re.IGNORECASE):
                iarg += 1
                template_file = sys.argv[iarg]
                print('Setting template_file = "%s"'%(template_file))
        elif arg_str == '-ext':
            if iarg+1 <= len(sys.argv)-1 and re.match(r'^([0-9]+)$', sys.argv[iarg+1]):
                iarg += 1
                ext_id = int(sys.argv[iarg])
                print('Setting ext_id = %d'%(ext_id))
        elif arg_str == '-extname':
            if iarg+1 <= len(sys.argv)-1:
                iarg += 1
                ext_name = str(sys.argv[iarg])
                print('Setting ext_name = %s'%(ext_name))
        else:
            if input_file == '':
                input_file = sys.argv[iarg]
                print('Setting input_file = "%s"'%(input_file))
            elif output_file == '':
                output_file = sys.argv[iarg]
                print('Setting output_file = "%s"'%(output_file))
            else:
                pass
        # 
        iarg += 1


    # Check user input and print usage if necessary
    if input_file == '' or \
       output_file == '' or \
       (len(new_fov) <= 0 and template_file == ''):
        usage()
        sys.exit()
    
    
    # Read template file if given
    if template_file != '':
        template_data, template_header = read_fits_data(template_file)
        template_wcs2D = WCS(template_header, naxis=2)
        template_pixsc = np.abs(proj_plane_pixel_scales(template_wcs2D)*3600.0)
        template_naxis1 = template_header['NAXIS1']
        template_naxis2 = template_header['NAXIS2']
        new_center_ra, new_center_dec = template_wcs2D.wcs_pix2world([(template_naxis1+1.0)/2.0], [(template_naxis2+1.0)/2.0], 1)
        if not np.isscalar(new_center_ra): new_center_ra = new_center_ra[0]
        if not np.isscalar(new_center_dec): new_center_dec = new_center_dec[0]
        new_center = parse_skycoord(str(new_center_ra)+'deg', str(new_center_dec)+'deg', return_skycoord_obj=True)
        new_fov = [template_naxis1*template_pixsc[0], template_naxis2*template_pixsc[1]]
        print('new_center', new_center)
        print('new_fov', new_fov)
    
    
    # Read fits file
    data, header = read_fits_data(input_file, ext_id=ext_id, ext_name=ext_name)
    wcs2D = WCS(header, naxis=2)
    pixsc = np.abs(proj_plane_pixel_scales(wcs2D)*3600.0)
    print('pixsc', pixsc, '[arcsec]')


    # compute new naxis
    new_naxis1 = int(np.ceil(new_fov[0] / pixsc[0]))
    new_naxis2 = int(np.ceil(new_fov[1] / pixsc[1]))
    print('new_naxis1', new_naxis1, 'new_naxis2', new_naxis2)
    paste_rect_x0 = 0
    paste_rect_x1 = new_naxis1-1
    paste_rect_y0 = 0
    paste_rect_y1 = new_naxis2-1
    cut_rect_x0 = 0
    cut_rect_x1 = header['NAXIS1']-1
    cut_rect_y0 = 0
    cut_rect_y1 = header['NAXIS2']-1
    if new_center is not None:
        new_center_x, new_center_y = wcs2D.wcs_world2pix([new_center.ra.deg], [new_center.dec.deg], 0) # origin=0
        if not np.isscalar(new_center_x): new_center_x = new_center_x[0]
        if not np.isscalar(new_center_y): new_center_y = new_center_y[0]
    else:
        new_center_x = (header['NAXIS1']-1) / 2.0
        new_center_y = (header['NAXIS2']-1) / 2.0
    print('new_center_x', new_center_x, 'new_center_y', new_center_y)

    cut_rect_x0 = int(np.ceil(new_center_x - new_naxis1 / 2.0))
    cut_rect_x1 = cut_rect_x0 + new_naxis1 - 1
    cut_rect_y0 = int(np.ceil(new_center_y - new_naxis2 / 2.0))
    cut_rect_y1 = cut_rect_y0 + new_naxis2 - 1
    print('cut_rect_x0, cut_rect_x1, cut_rect_y0, cut_rect_y1 =', cut_rect_x0, cut_rect_x1, cut_rect_y0, cut_rect_y1)
    print('paste_rect_x0, paste_rect_x1, paste_rect_y0, paste_rect_y1 =', paste_rect_x0, paste_rect_x1, paste_rect_y0, paste_rect_y1)
    if cut_rect_x0 < 0:
        paste_rect_x0 = -cut_rect_x0
        cut_rect_x0 = 0
    if cut_rect_x1 > header['NAXIS1']-1:
        paste_rect_x1 = (new_naxis1-1)-(cut_rect_x1-(header['NAXIS1']-1))
        cut_rect_x1 = header['NAXIS1']-1
    if cut_rect_y0 < 0:
        paste_rect_y0 = -cut_rect_y0
        cut_rect_y0 = 0
    if cut_rect_y1 > header['NAXIS2']-1:
        paste_rect_y1 = (new_naxis2-1)-(cut_rect_y1-(header['NAXIS2']-1))
        cut_rect_y1 = header['NAXIS2']-1
    print('cut_rect_x0, cut_rect_x1, cut_rect_y0, cut_rect_y1 =', cut_rect_x0, cut_rect_x1, cut_rect_y0, cut_rect_y1)
    print('paste_rect_x0, paste_rect_x1, paste_rect_y0, paste_rect_y1 =', paste_rect_x0, paste_rect_x1, paste_rect_y0, paste_rect_y1)
    
    
    # prepare output data
    new_data_shape = copy.copy(list(data.shape))
    new_data_shape[-1] = new_naxis1
    new_data_shape[-2] = new_naxis2
    new_data = np.full(new_data_shape, fill_value=np.nan)


    # trim fov pixels
    #if len(data.shape) == 3:
    #    new_data = np.full((data.shape[0], new_naxis2, new_naxis1), fill_value=np.nan)
    #    print('new_data[:, %d:%d, %d:%d] = data[:, %d:%d, %d:%d]'%(paste_rect_y0, paste_rect_y1+1, paste_rect_x0, paste_rect_x1+1, cut_rect_y0, cut_rect_y1+1, cut_rect_x0, cut_rect_x1+1))
    #    new_data[:, paste_rect_y0:paste_rect_y1+1, paste_rect_x0:paste_rect_x1+1] = data[:, cut_rect_y0:cut_rect_y1+1, cut_rect_x0:cut_rect_x1+1]
    #else:
    #    new_data = np.full((new_naxis2, new_naxis1), fill_value=np.nan)
    #    print('new_data[%d:%d, %d:%d] = data[%d:%d, %d:%d]'%(paste_rect_y0, paste_rect_y1+1, paste_rect_x0, paste_rect_x1+1, cut_rect_y0, cut_rect_y1+1, cut_rect_x0, cut_rect_x1+1))
    #    new_data[paste_rect_y0:paste_rect_y1+1, paste_rect_x0:paste_rect_x1+1] = data[cut_rect_y0:cut_rect_y1+1, cut_rect_x0:cut_rect_x1+1]
    for ndidx in np.ndindex(data.shape[0:-2]):
        new_data[ndidx][paste_rect_y0:paste_rect_y1+1, paste_rect_x0:paste_rect_x1+1] = data[ndidx][cut_rect_y0:cut_rect_y1+1, cut_rect_x0:cut_rect_x1+1]



    # Output
    new_header = copy.deepcopy(header)
    new_header['NAXIS1'] = new_naxis1
    new_header['CRPIX1'] = header['CRPIX1']-cut_rect_x0+paste_rect_x0
    new_header['NAXIS2'] = new_naxis2
    new_header['CRPIX2'] = header['CRPIX2']-cut_rect_y0+paste_rect_y0
    new_header['HISTORY'] = ''
    new_header['HISTORY'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ' ' + time.strftime('%Z')
    history_str = 'Cutout with field of view %s x %s [arcsec] '%(new_fov[0], new_fov[1])
    if new_center is not None:
        history_str += 'at the new center RA Dec %s %s '%(new_center.ra.deg, new_center.dec.deg )
    if ext_id is not None:
        history_str += 'from the extension id %d of the input FITS file "%s".'%(ext_id, input_file)
    elif ext_name is not None:
        if 'EXTNAME' in new_header:
            del new_header['EXTNAME']
        history_str += 'from the extension name %s of the input FITS file "%s".'%(ext_name, input_file)
    else:
        history_str += 'from the input FITS file "%s".'%(input_file)
    new_header['HISTORY'] = history_str
    new_header['HISTORY'] = ''

    if set_output_naxis_3:
        data, header = trim_data_dimension(data, header, 3)
    
    if set_output_naxis_2:
        data, header = trim_data_dimension(data, header, 2)

    write_fits_data_to_file(new_data, new_header, output_file, 
        set_output_bitpix_32=set_output_bitpix_32)









