#!/usr/bin/env python 
# 

"""
Measure and renormalize strips in JWST images. 

"""

import os, sys, re, copy, shutil, glob, datetime
import click
import numpy as np
import tqdm
from astropy.io import fits
from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter, gaussian_filter, binary_dilation, generate_binary_structure
from tqdm import tqdm

# jwst
from jwst.datamodels import ImageModel, FlatModel, dqflags
from jwst.datamodels.dqflags import pixel as dqflags_pixel
from jwst.flatfield.flat_field import do_correction, do_flat_field, apply_flat_field
from stdatamodels import util
import crds

# code name and version
CODE_NAME = 'util_remove_stripes_along_four_angles.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221109'
CODE_HOMEPAGE = ''

# logging
import logging
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)

# Rect
class Rect():
    """A Rectangle."""
    def __init__(self, x0, y0, x1, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
    def __str__(self):
        return 'Rect({}, {}, {}, {}, width={}, height={})'.format(
            self.x0, self.y0, self.x1, self.y1, self.width(), self.height())
    def width(self):
        return self.x1-self.x0+1
    def height(self):
        return self.y1-self.y0+1

# Scan
class Scan():
    """Define a Scan of pixels in an image."""
    def __init__(self, image, px, py, angle=None, rect=None, dqmask=None):
        self.px = px
        self.py = py
        self.ipx = np.round(self.px).astype(int)
        self.ipy = np.round(self.py).astype(int)
        self.angle = angle
        self.rect = rect
        if self.rect is None:
            ny, nx = image.shape
            self.rect = Rect(0, 0, nx-1, ny-1)
        ixmin, ixmax = self.rect.x0, self.rect.x1
        iymin, iymax = self.rect.y0, self.rect.y1
        self.imask = np.logical_and.reduce((self.ipy>=iymin, self.ipx>=ixmin, self.ipy<=iymax, self.ipx<=ixmax))
        self.iipy = self.ipy[self.imask]
        self.iipx = self.ipx[self.imask]
        self.data = image[self.iipy, self.iipx].flatten()
        nanmask = np.isnan(self.data)
        if dqmask is not None:
            nanmask = np.logical_or(nanmask, 
                np.invert(dqmask[self.iipy, self.iipx].flatten()))
                # dqmask>0 means good pixel here, invert means bad pixels.
        if len(self.data) > np.count_nonzero(nanmask):
            self.mean, self.median, self.std = sigma_clipped_stats(self.data, 
                mask=nanmask)
                # mask here means bad pixels to exclude.
        else:
            self.mean, self.median, self.std = np.nan, np.nan, np.nan
    def __str__(self):
        if len(self.data) > 0:
            return 'Scan(({}, {}) -> ({}, {}), angle={}, npix={}, median={})'.format(
                self.iipx[0], self.iipy[0], self.iipx[-1], self.iipy[-1], self.angle, len(self.data), self.median)
        else:
            return 'Scan(({}, {}) -> ({}, {}), angle={}, npix={}, median={})'.format(
                self.ipx[0], self.ipy[0], self.ipx[-1], self.ipy[-1], self.angle, len(self.data), self.median)
    def subtract_median_for_image(self, image, median_image=None, masked_image=None):
        if not np.isnan(self.median):
            image[self.iipy, self.iipx] -= self.median
            if masked_image is not None:
                masked_image[self.iipy, self.iipx] -= self.median
        if median_image is not None:
            median_image[self.iipy, self.iipx] += self.median
            return image, median_image
        return image



# Draw scans
def draw_scans(image, angle, rect=None, dqmask=None):
    """Draw scans along an angle for the whole image.
    """
    scans = []
    ny, nx = image.shape
    cos = np.cos(np.deg2rad(angle))
    sin = np.sin(np.deg2rad(angle))
    tan = np.tan(np.deg2rad(angle))
    if np.isclose(cos, 0, rtol=0, atol=1e-6): # angle==90
        for iscan in range(nx):
            px = np.zeros(ny) + iscan
            py = np.arange(ny)
            scans.append(Scan(image, px, py, angle=angle, rect=rect, dqmask=dqmask))
    else:
        nscan = ny + int(np.floor((nx-1)*np.abs(tan)))
        for iscan in range(nscan):
            y0 = -np.floor(nx*tan) + iscan if (tan >= 0.0) else 0 + iscan
            dy = tan
            #if tan < 0.0:
            #    print(y0, y0+(nx-1)*dy, iscan)
            px = np.arange(nx)
            py = np.linspace(y0, y0+(nx-1)*dy, num=nx, endpoint=True)
            scans.append(Scan(image, px, py, angle=angle, rect=rect, dqmask=dqmask))
    return scans



# Apply flat
def apply_flat_to_model(model):
    logger.info('Applying flat')
    import crds
    try:
        crds_context = os.environ['CRDS_CONTEXT']
    except KeyError:
        crds_context = crds.get_default_context()
    # pull flat from CRDS using the current context
    #<DZLIU># enabling both NIRCAM and MIRI, 
    #<DZLIU># following https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/reference_files.html
    instrument_name = model.meta.instrument.name
    if instrument_name.upper() in ['NIRCAM', 'NIRISS']: #<DZLIU>#
        crds_dict = {'INSTRUME':instrument_name, #<DZLIU>#
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'PUPIL':model.meta.instrument.pupil, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    elif instrument_name.upper() == 'MIRI': #<DZLIU>#
        crds_dict = {'INSTRUME':instrument_name, #<DZLIU>#
                     'DETECTOR':model.meta.instrument.detector, 
                     'FILTER':model.meta.instrument.filter, 
                     'BAND':model.meta.instrument.band, 
                     'READPATT':model.meta.exposure.readpatt, 
                     'SUBARRAY':model.meta.subarray.name, 
                     'DATE-OBS':model.meta.observation.date,
                     'TIME-OBS':model.meta.observation.time}
    else:
        raise ValueError("Error! The `instrument_name` (model.meta.instrument.name) can only be NIRCAM and MIRI!")
    logger.info('crds_dict: {}'.format(crds_dict)) #<DZLIU>#
    flats = crds.getreferences(crds_dict, reftypes=['flat'], 
                               context=crds_context)
    # if the CRDS loopup fails, should return a CrdsLookupError, but 
    # just in case:
    try:
        flatfile = flats['flat']
    except KeyError:
        logger.error('Flat was not found in CRDS with the parameters: {}'.format(crds_dict))
        exit()

    logger.info('Using flat: %s'%(os.path.basename(flatfile)))
    with FlatModel(flatfile) as flat:
        # use the JWST Calibration Pipeline flat fielding Step 
        model, applied_flat = do_correction(model, flat)
    
    return model, applied_flat




# main
@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', required=False, type=click.Path(exists=False), default=None)
@click.option('--mask-image', type=click.Path(exists=True), default=None, help='A mask image containing pixels which will not be analyzed, e.g., due to bright emission. If not given we will automatically determine the mask.')
@click.option('--mask-threshold', type=float, default=None, help='Pixels with values above this threshold in the mask image will be not be analyzed.')
@click.option('--apply-flat', is_flag=True, default=None, help='Apply CRDS flat correction before destriping?')
def main(
        input_file, 
        output_file, 
        mask_image, 
        mask_threshold, 
        apply_flat, 
    ):
    """
    Input rate image. 
    """
    with ImageModel(input_file) as model:

        # check if already processed
        for entry in model.history:
            for k,v in entry.items():
                if v.startswith('Removed stripes'):
                    logger.info('{!r} already had stripes removed. Skipping!'.format(input_file))
                    return
        
        # Apply flat (and dq)
        apply_flat_status = None
        if apply_flat:
            model, applied_flat = apply_flat_to_model(model)
            apply_flat_status = model.meta.cal_step.flat_field
        
        # Get image and dqmask
        # dqmask = bitfield_to_boolean_mask(
        #         model.dq,
        #         interpret_bit_flags('~DO_NOT_USE+NON_SCIENCE', # SkyMatchStep().dqbits is '~DO_NOT_USE+NON_SCIENCE' in default
        #             flag_name_map=dqflags_pixel), 
        #         good_mask_value=1,
        #         dtype=np.uint8
        #     ) * np.isfinite(model.data) # following "jwst/skymatch/skymatch_step.py", this is valid data
        dqmask = None # 
        image = model.data.copy()
        #image[~dqmask] = np.nan
        
        # Get image shape
        logger.info('image.shape {}'.format(image.shape))
        ny, nx = image.shape
        
        # Assign rect area(s), each rect will be destriped separately
        rects = []
        instrument_name = model.meta.instrument.name
        if instrument_name.upper() == 'NIRCAM':
            rects.append(Rect(0, 0, nx-1, ny-1)) # x0, y0, x1, y1
        elif instrument_name.upper() == 'MIRI':
            rects.append(Rect(356, 0, nx-1, ny-1)) # MIRI imager area
            rects.append(Rect(7, 746, 280, 1019)) # MIRI Lyot coronagraph
            rects.append(Rect(7, 464, 233, 683)) # MIRI 4QPM coronagraph
            rects.append(Rect(7, 243, 233, 463)) # MIRI 4QPM coronagraph
            rects.append(Rect(7,  21, 233, 242)) # MIRI 4QPM coronagraph
        else:
            raise Exception('The input data can only be NIRCAM or MIRI data!')
        
        # Mask out sources
        image_masked = copy.deepcopy(image)
        if mask_image is not None:
            logger.info('Using the input mask image {!r} to mask out source flux'.format(mask_image))
            mask_image_data = fits.getdata(mask_image)
            if mask_threshold is not None:
                mask_image_bool = (mask_image_data > mask_threshold)
                n_masked = np.count_nonzero(mask_image_bool)
                logger.info('Masked {} pixels with values above the mask threshold {}'.format(n_masked, mask_threshold))
            else:
                mask_image_bool = mask_image_data.astype(bool)
                n_masked = np.count_nonzero(mask_image_bool)
                logger.info('Masked {} pixels'.format(n_masked))
            image_masked[mask_image_bool] = np.nan
        else:
            # overall sigma clip
            logger.info('Auto masking with 3-sigma clipping')
            overall_mean, overall_median, overall_stddev = sigma_clipped_stats(
                image, mask=np.isnan(image))
            mask_image_bool = (np.abs(image_masked - overall_mean) > 3.0 * overall_stddev)
            mask_image_float = mask_image_bool.astype(float)
            mask_image_float = binary_dilation(mask_image_float, generate_binary_structure(2,1), iterations=2)
            mask_image_bool = mask_image_float.astype(bool)
            image_masked[mask_image_bool] = np.nan
            
        # close model
        #model.close()
    
    
    # Grid
    #gy, gx = np.mgrid[0:ny, 0:nx]
    
    # Four angles to scan the image
    angles = [0.0, 90.0, 30.0, 150.0] # order is important!
    
    bkgs = np.zeros((len(rects), len(angles), ny, nx))
    
    # measure striping per rect, per angle
    #for irect in tqdm(range(len(rects))):
    #    for iangle in tqdm(range(len(angles)), leave=False):
    for irect in range(len(rects)):
        for iangle in range(len(angles)):
            rect = rects[irect]
            angle = angles[iangle]
            #bkg = np.zeros(image.shape)
            #msk = np.zeros(image.shape)
            scans = draw_scans(image_masked, angle, rect=rect, dqmask=dqmask)
            for iscan in range(len(scans)):
                scan = scans[iscan]
                scan.subtract_median_for_image(image, median_image=bkgs[irect, iangle], masked_image=image_masked)
                # 20221110 also need to subtract this image_masked
                #if iangle == 2 and iscan > 800 and iscan < 1000 and scan.median > 0.2:
                #    print('iscan', iscan, 'scan', scan, 'scan.data', scan.data.tolist())
                #    raise Exception('')
                # if iscan > 10:
                #     print('iscan', iscan, 'scan', scan)
                #     break
    
    bkg = np.sum(bkgs, axis=(0, 1))
    
    if output_file is None:
        output_file = input_file
    if os.path.abspath(output_file) == os.path.abspath(input_file):
        logger.info('Updating the input file in-place!')
        output_orig = re.sub(r'\.fits$', r'', output_file) + '_before_removing_stripes.fits'
        if os.path.isfile(output_orig):
            logger.info('Found backup file {}, skipping backup.'.format(output_orig))
        else:
            logger.info('Backing up input file as {}'.format(output_orig))
            shutil.copy2(output_file, output_orig) # do not overwrite
    else:
        if os.path.isfile(output_file):
            shutil.move(output_file, output_file+'.backup')
        elif not os.path.isdir(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        shutil.copy2(input_file, output_file) # copy input to output then update the model.data
    
    out_timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with ImageModel(output_file) as out_model:
        validmask = np.isfinite(out_model.data)
        out_model.data[validmask] = image[validmask] # only copy valid pixel values
        #out_model.data = image
        # check flat
        if apply_flat:
            out_model.meta.cal_step.flat_field = apply_flat_status
        # add history entry following CEERS
        stepdescription = f'Removed stripes with {CODE_NAME} ({out_timestamp})'
        software_dict = {'name':CODE_NAME,
                         'author':CODE_AUTHOR,
                         'version':CODE_VERSION,
                         'homepage':CODE_HOMEPAGE}
        substr = util.create_history_entry(stepdescription, software=software_dict)
        out_model.history.append(substr)
        logger.info('Adding model history: {}'.format(substr))
        out_model.save(output_file)
        logger.info('Saved destriped data into {!r}'.format(output_file))
    
    
    output_bkg = re.sub(r'\.fits$', r'', output_file) + '_removing_stripes_bkg.fits'
    shutil.copy2(output_file, output_bkg)
    with ImageModel(output_bkg) as out_bkg_model:
        out_bkg_model.data = bkg
        out_bkg_model.save(output_bkg)
        logger.info('Saved destriping background into {!r}'.format(output_bkg))
    
    
    output_bkgs = re.sub(r'\.fits$', r'', output_file) + '_removing_stripes_bkgs.fits'
    header0 = fits.getheader(output_bkg, 0)
    header1 = fits.getheader(output_bkg, 1)
    hdulist = [fits.PrimaryHDU(data=None, header=header0)]
    for irect in range(len(rects)):
        for iangle in range(len(angles)):
            rect = rects[irect]
            angle = angles[iangle]
            header2 = copy.deepcopy(header1)
            header2['EXTNAME'] = 'RECT_{}_ANGLE_{}'.format(irect+1, iangle+1)
            header2['RECT'] = str(rect)
            header2['ANGLE'] = angle
            hdulist.append(fits.ImageHDU(data=bkgs[irect, iangle], header=header2))
    hdul = fits.HDUList(hdulist)
    hdul.writeto(output_bkgs, overwrite=True)
    hdul.close()
    logger.info('Saved destriping backgrounds into {!r}'.format(output_bkgs))
    
    
    output_mask = re.sub(r'\.fits$', r'', output_file) + '_removing_stripes_mask.fits'
    shutil.copy2(output_file, output_mask)
    with ImageModel(output_mask) as out_mask_model:
        out_mask_model.data = mask_image_bool.astype(int)
        out_mask_model.save(output_mask)
        logger.info('Saved destriping mask into {!r}'.format(output_mask))
    
    
    # output_mask = re.sub(r'\.fits$', r'', output_file) + '_removing_stripes_masked.fits'
    # shutil.copy2(output_file, output_mask)
    # with ImageModel(output_mask) as out_mask_model:
    #     out_mask_model.data = image_masked
    #     out_mask_model.save(output_mask)
    #     logger.info('Saved destriping masked data into {!r}'.format(output_mask))
    



if __name__ == '__main__':
    main()



