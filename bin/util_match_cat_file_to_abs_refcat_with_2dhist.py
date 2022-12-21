#!/usr/bin/env python
# 
import os, sys, re, copy, json, shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_area
from astropy.coordinates import SkyCoord, FK5, match_coordinates_sky
from astropy.nddata.utils import Cutout2D
from astropy.visualization import simple_norm
from astropy.stats import sigma_clipped_stats
from collections import OrderedDict
from matplotlib.patches import Ellipse
import click
import warnings
warnings.simplefilter('ignore', FITSFixedWarning)
import logging
logging.basicConfig(level='INFO')
logger = logging.getLogger()


default_output_suffix = '_new2dhist'
default_outlier_sigma = 5.0
default_output_abs_refcat = None
default_input_index_base = 1
default_output_index_base = 0

def match_cat_file_to_abs_refcat_with_2dhist(
        cat_file, # This should be the tweakreg.catfile format, i.e., two columns, space-sep., 
                  # for a FITS-image file name and a source-finding catalog file name.
                  # Each source-finding catalog file should contain 'x' and 'y' columns, space-sep.
                  # Each FITS-image file must has an 'SCI' extension whose header contains WCS information. 
        abs_refcat, # The tweakreg.abs_refcat catalog file. 
                    # It must include 'RA', 'DEC', and 'ID' (or 'phot_id') columns. 
        output_suffix = default_output_suffix, 
        outlier_sigma = default_outlier_sigma, # filter out outliers of large offsets
        output_abs_refcat = default_output_abs_refcat, # also output the filtered abs_refcat
        input_index_base = default_input_index_base, # the input image catalog listed in catfile should have 1-based x y coordinates
        output_index_base = default_output_index_base, # the output catalog x y coordinate base, default is 0 because jwst pipeline tweakwcs uses that.
    ):

    # check user input
    if output_suffix is None:
        output_suffix = '_new2dhist'
    if outlier_sigma is None:
        outlier_sigma = 5.0

    # read cat_file
    cats = OrderedDict()
    catpaths = OrderedDict()
    with open(cat_file, 'r') as fp:
        for line in fp.readlines():
            imfile, catpath = line.split()
            cats[imfile] = Table.read(catpath)
            catpaths[imfile] = catpath
    refcat = Table.read(abs_refcat)

    # make 2dhist plot
    ncats = len(cats)
    maxsep = 0.7 * u.arcsec
    ipanel = 0
    ncols = 2
    nrows = int(np.ceil(ncats/ncols))
    logger.info('Creating 2Dhist plot with ncols {} nrows {} ncats {}'.format(ncols, nrows, ncats))
    fig = plt.figure(figsize=(2.5+ncols*2.0, 1.0+nrows*2.0))
    axes = []
    for irow in range(nrows):
        for icol in range(ncols):
            axes.append(fig.add_subplot(nrows, ncols, irow*ncols+icol+1))
    axes = np.array(axes)
    #fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(1.0+ncols*2.0, 1.0+nrows*2.0), 
    #                         gridspec_kw=dict(left=0.09, bottom=0.06, right=0.98, top=0.98, wspace=1.0, hspace=0.4)
    #                        )
    #check_numb = 10
    #check_radecs = {}
    catmasks = OrderedDict()
    offsets = OrderedDict()
    matched_abs_refcat_mask = np.full(len(refcat), fill_value=False, dtype=bool)
    matched_abs_refcat_id = OrderedDict()
    matched_abs_refcat_x = OrderedDict()
    matched_abs_refcat_y = OrderedDict()
    for imfile in cats:
        header = fits.getheader(imfile, 'SCI')
        wcs = WCS(header, naxis=2)
        pixsc = np.sqrt(proj_plane_pixel_area(wcs)) * 3600.
        cat = cats[imfile]
        x, y = cat['x'], cat['y']
        if input_index_base == 0:
            x += 1
            y += 1
        ra, dec = wcs.all_pix2world(x, y, 1) # wcs.wcs_pix2world(x, y) is not enough # internally I use DS9 1-based pixcoord
        catcoords = SkyCoord(ra*u.deg, dec*u.deg, frame=FK5)
        refra = refcat['RA'].data
        refdec = refcat['DEC'].data
        refid = refcat['phot_id'].data
        refcatcoords = SkyCoord(refra*u.deg, refdec*u.deg, frame=FK5)
        refx, refy = wcs.all_world2pix(refra, refdec, 1) # internally I use DS9 1-based pixcoord
        idx, d2d, d3d = match_coordinates_sky(catcoords, refcatcoords)
        matches = np.argwhere(np.logical_and(idx>=0, d2d<=maxsep)).ravel()
        d_ra = (ra[matches]-refra[idx[matches]]) * np.cos(np.deg2rad(refdec[idx[matches]])) * 3600.
        d_dec = (dec[matches]-refdec[idx[matches]]) * 3600.
        d_px = (x[matches]-refx[idx[matches]])
        d_py = (y[matches]-refy[idx[matches]])
        lim = np.max(np.abs([d_ra, d_dec])) * 1.2 # axes range +20%
        # 
        for imatch in matches:
            key = refid[imatch]
            #if key not in check_radecs:
            #    check_radecs[key] = OrderedDict()
            #    check_radecs[key]['x'] = np.full([ncats,], fill_value=np.nan)
            #    check_radecs[key]['y'] = np.full([ncats,], fill_value=np.nan)
            #    check_radecs[key]['ra'] = np.full([ncats,], fill_value=np.nan)
            #    check_radecs[key]['dec'] = np.full([ncats,], fill_value=np.nan)
            #    check_radecs[key]['refra'] = np.full([ncats,], fill_value=np.nan)
            #    check_radecs[key]['refdec'] = np.full([ncats,], fill_value=np.nan)
            #check_radecs[key]['x'][ipanel] = x[imatch]
            #check_radecs[key]['y'][ipanel] = y[imatch]
            #check_radecs[key]['ra'][ipanel] = ra[imatch]
            #check_radecs[key]['dec'][ipanel] = dec[imatch]
            #check_radecs[key]['refra'][ipanel] = refra[idx[imatch]]
            #check_radecs[key]['refdec'][ipanel] = refdec[idx[imatch]]
            # 
            matched_abs_refcat_mask[idx[imatch]] = True
        # 
        ax = axes.ravel()[ipanel]
        ax.scatter(d_ra, d_dec, s=3, alpha=0.9)
        ax.set_xlim([lim, -lim])
        ax.set_ylim([-lim, lim])
        # 
        d_ra_mean, d_ra_median, d_ra_sigma = sigma_clipped_stats(d_ra)
        d_dec_mean, d_dec_median, d_dec_sigma = sigma_clipped_stats(d_dec)
        d_px_mean, d_px_median, d_px_sigma = sigma_clipped_stats(d_px)
        d_py_mean, d_py_median, d_py_sigma = sigma_clipped_stats(d_py)
        offsets[imfile] = OrderedDict()
        offsets[imfile]['mean_arcsec'] = [d_ra_mean, d_dec_mean]
        offsets[imfile]['mean_pixel'] = [d_px_mean, d_py_mean] # in detector frame
        offsets[imfile]['median_arcsec'] = [d_ra_median, d_dec_median]
        offsets[imfile]['median_pixel'] = [d_px_median, d_py_median] # in detector frame
        offsets[imfile]['sigma_arcsec'] = [d_ra_sigma, d_dec_sigma]
        offsets[imfile]['sigma_pixel'] = [d_px_sigma, d_py_sigma] # in detector frame
        patch = Ellipse(xy=[d_ra_mean, d_dec_mean], width=2.*outlier_sigma*d_ra_sigma, height=2.*outlier_sigma*d_dec_sigma, 
                        angle=0.0, ec='red', fc='none', alpha=0.8)
        patch = Ellipse(xy=[d_ra_median, d_dec_median], width=2.*outlier_sigma*d_ra_sigma, height=2.*outlier_sigma*d_dec_sigma, 
                        angle=0.0, ec='red', fc='none', alpha=0.8)
        ax.add_patch(patch)
        ax.text(d_ra_mean, d_dec_mean+outlier_sigma*d_dec_sigma, '{:+.4f}, {:+.4f}'.format(d_ra_mean, d_dec_mean),
            ha='center', va='bottom', color='red')
        # 
        # identify outliers
        outliers = ( ((d_ra-d_ra_mean)/(outlier_sigma*d_ra_sigma))**2 + ((d_dec-d_dec_mean)/(outlier_sigma*d_dec_sigma))**2 > 1.0)
        ax.scatter(d_ra[outliers], d_dec[outliers], s=10, marker='*', alpha=0.7, 
                   fc='none', ec='orangered')
        catmask = np.full([len(ra),], fill_value=False, dtype=bool)
        catmask[matches] = np.invert(outliers) # remove outliers, keep good ones
        catmasks[imfile] = catmask
        matched_abs_refcat_id[imfile] = refid[idx[matches][~outliers]]
        matched_abs_refcat_x[imfile] = refx[idx[matches][~outliers]]
        matched_abs_refcat_y[imfile] = refy[idx[matches][~outliers]]
        # 
        ax.plot([-lim, lim], [0., 0.], color='k', ls='dashed', alpha=0.7)
        ax.plot([0., 0.], [-lim, lim], color='k', ls='dashed', alpha=0.7)
        # 
        ax.set_title(imfile, fontsize='small')
        ax.set_aspect(1.0)
        # 
        icol = ipanel%ncols
        irow = ipanel//ncols
        if icol == 0:
            ax.set_ylabel(r'$\Delta$ DEC [arcsec]')
        if irow == nrows-1:
            ax.set_xlabel(r'$\Delta$ RA [arcsec]')
        # 
        ipanel += 1
    
    # remove outliers and write each individual source-finding catalog csv files
    newcatpaths = {}
    for imfile in cats:
        cat = cats[imfile]
        catpath = catpaths[imfile]
        catmask = catmasks[imfile]
        newcat = cat[catmask] # remove outliers, keep good ones, mask means good here
        chkcat = Table({'x': matched_abs_refcat_x[imfile], 'y': matched_abs_refcat_y[imfile]}) # matched refcat for check
        # 
        # output ds9 pix region
        newcatpath = os.path.splitext(catpath)[0] + output_suffix + '.ds9.pix.reg'
        if os.path.isfile(newcatpath):
            shutil.move(newcatpath, newcatpath+'.backup')
        with open(newcatpath, 'w') as fp:
            fp.write('# DS9 region file\n')
            fp.write('global color=green\n')
            fp.write('image\n')
            for k in range(len(newcat)):
                fp.write('circle({}, {}, {}")\n'.format(newcat['x'][k], newcat['y'][k],0.35))
        logger.info('Output to {!r}'.format(newcatpath))
        # 
        # output csv x y pixcoord
        if output_index_base == 0:
            newcat['x'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
            newcat['y'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
        newcatpath = os.path.splitext(catpath)[0] + output_suffix + '.csv'
        if os.path.isfile(newcatpath):
            shutil.move(newcatpath, newcatpath+'.backup')
        newcat.write(newcatpath, format='csv')
        newcatpaths[imfile] = newcatpath
        logger.info('Output to {!r}'.format(newcatpath))
        # 
        # output matched refcat ds9 pix region
        chkcatpath = os.path.splitext(catpath)[0] + output_suffix + '_matched_refcat.ds9.pix.reg'
        if os.path.isfile(chkcatpath):
            shutil.move(chkcatpath, chkcatpath+'.backup')
        with open(chkcatpath, 'w') as fp:
            fp.write('# DS9 region file\n')
            fp.write('global color=yellow\n')
            fp.write('image\n')
            for k in range(len(newcat)):
                fp.write('circle({}, {}, {}")\n'.format(newcat['x'][k], newcat['y'][k],0.35))
        logger.info('Output to {!r}'.format(chkcatpath))
        # 
        # output matched refcat
        chkcatpath = os.path.splitext(catpath)[0] + output_suffix + '_matched_refcat.csv'
        if os.path.isfile(chkcatpath):
            shutil.move(chkcatpath, chkcatpath+'.backup')
        chkcat = Table({'x': matched_abs_refcat_x[imfile], 'y': matched_abs_refcat_y[imfile]})
        if output_index_base == 0:
            chkcat['x'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
            chkcat['y'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
        chkcat.write(chkcatpath, format='csv', overwrite=True)
        logger.info('Output to {!r}'.format(chkcatpath))
    
    # 
    fig.tight_layout(w_pad=3.5)
    
    # write the big two-column cat_file
    cat_filebase, cat_filesuffix = os.path.splitext(cat_file)
    out_file = cat_filebase + output_suffix + cat_filesuffix
    out_json = cat_filebase + output_suffix + '.json'
    out_png = cat_filebase + output_suffix + '.png'
    out_pdf = cat_filebase + output_suffix + '.png'
    fig.savefig(out_png, dpi=300)
    logger.info('Output to {!r}'.format(out_png))
    fig.savefig(out_pdf, dpi=300)
    logger.info('Output to {!r}'.format(out_pdf))
    if os.path.isfile(out_json):
        shutil.move(out_json, out_json+'.backup')
    with open(out_json, 'w') as fp:
        json.dump(offsets, fp, indent=4)
    logger.info('Output to {!r}'.format(out_json))
    if os.path.isfile(out_file):
        shutil.move(out_file, out_file+'.backup')
    with open(out_file, 'w') as fp:
        for imfile in cats:
            newcatpath = newcatpaths[imfile]
            fp.write('{} {}\n'.format(imfile, newcatpath))
    logger.info('Output to {!r}'.format(out_file))
    
    # 
    if output_abs_refcat is None:
        output_abs_refcat = os.path.splitext(abs_refcat)[0] + output_suffix + '.fits' # FITS-format catalog
    if output_abs_refcat != '':
        if 'name' not in refcat.meta: 
            refcat.meta['name'] = os.path.splitext(os.path.basename(abs_refcat))[0]
        refcat = refcat[matched_abs_refcat_mask]
        if os.path.isfile(output_abs_refcat):
            shutil.move(output_abs_refcat, output_abs_refcat+'.backup')
        refcat.write(output_abs_refcat)
        logger.info('Output to {!r}'.format(output_abs_refcat))
    
    # return
    return out_file





@click.command()
@click.argument('cat_file', type=click.Path(exists=True))
@click.argument('abs_refcat', type=click.Path(exists=True))
@click.option('--output-suffix', type=str, default=default_output_suffix)
@click.option('--outlier-sigma', type=float, default=default_outlier_sigma, help='filter outlier by large offsets, in sigma of the offset distribution.')
@click.option('--output-abs-refcat', type=click.Path(exists=False), default=default_output_abs_refcat, help='also output the filtered abs_refcat.')
@click.option('--input-index-base', type=click.Choice([0, 1]), default=default_input_index_base, help='x y pixcoord base for input sextractor catalog listed in catfile. The default is 1 because sextractor uses 1-based pixcoord.')
@click.option('--output-index-base', type=click.Choice([0, 1]), default=default_output_index_base, help='x y pixcoord base for output matched csv. The default is 0 because jwst pipeline tweakwcs correctors.py uses 0-based pixcoord.')
def main(
        cat_file,
        abs_refcat,
        output_suffix,
        outlier_sigma,
        output_abs_refcat, 
        input_index_base, 
        output_index_base, 
    ):

    match_cat_file_to_abs_refcat_with_2dhist(
        cat_file, 
        abs_refcat, 
        output_suffix, 
        outlier_sigma,
        output_abs_refcat, 
        input_index_base, 
        output_index_base, 
    )





if __name__ == '__main__':
    main()









