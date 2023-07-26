#!/usr/bin/env python
# 
# 2023-03-17 fixed ds9 region bug, not affecting csv files.
# 
import os, sys, re, copy, json, shutil
import numpy as np
import matplotlib as mpl
mpl.rcParams['figure.max_open_warning'] = 0
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
import click
import warnings
warnings.simplefilter('ignore', FITSFixedWarning)
import logging
logging.basicConfig(level='INFO')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
logger = logging.getLogger()


default_search_radius = 1.0 # arcsec
default_output_suffix = '_new2dhist'
default_outlier_sigma = 5.0
default_initial_offset = [0.0, 0.0] # arcsec
default_output_dir = None
default_output_abs_refcat = None
default_input_index_base = '1'
default_output_index_base = '0'

def match_cat_file_to_abs_refcat_with_2dhist(
        cat_file, # This should be the tweakreg.catfile format, i.e., two columns, space-sep., 
                  # for a FITS-image file name and a source-finding catalog file name.
                  # Each source-finding catalog file should contain 'x' and 'y' columns, space-sep.
                  # Each FITS-image file must has an 'SCI' extension whose header contains WCS information. 
        abs_refcat, # The tweakreg.abs_refcat catalog file. 
                    # It must include 'RA', 'DEC', and 'ID' (or 'phot_id') columns. 
        search_radius = default_search_radius, 
        output_suffix = default_output_suffix, 
        outlier_sigma = default_outlier_sigma, # filter out outliers of large offsets
        initial_offset = default_initial_offset, 
        output_dir = default_output_dir, 
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
            if line.strip() != '' and not line.startswith('#'):
                imfile, catpath = line.split()
                cats[imfile] = Table.read(catpath)
                catpaths[imfile] = catpath
    if abs_refcat.endswith('.txt') or abs_refcat.endswith('.dat'):
        refcat = Table.read(abs_refcat, format='ascii')
    else:
        refcat = Table.read(abs_refcat)

    # make 2dhist plot
    ncats = len(cats)
    maxsep = search_radius * u.arcsec
    ipanel = 0
    ncols = 2
    nrows = 5
    npages = int(np.ceil(ncats/(ncols*nrows)))
    logger.info('Creating 2Dhist plot with ncols {} * nrows {} * npages {} = ncats {}'.format(ncols, nrows, npages, ncats))
    # 
    figs = []
    axes = []
    for ipage in range(npages):
        fig = plt.figure(figsize=(2.5+ncols*2.0, 1.0+nrows*2.0)) # figsize in inches
        for irow in range(nrows):
            for icol in range(ncols):
                ax = fig.add_subplot(nrows, ncols, irow*ncols+icol+1)
                ax.axis('off')
                axes.append(ax)
        figs.append(fig)
    # 
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
        if str(input_index_base) == '0':
            x += 1
            y += 1
        ra, dec = wcs.all_pix2world(x, y, 1) # wcs.wcs_pix2world(x, y) is not enough # internally I use DS9 1-based pixcoord
        if initial_offset is not None:
            ioffra, ioffdec = initial_offset
            if ioffra != 0 and ioffdec != 0:
                logger.info('Applying initial offset {} {} arcsec for {}'.format(ioffra, ioffdec, imfile))
                dec -= ioffdec / 3600.0
                ra -= ioffra / 3600.0 / np.cos(np.deg2rad(dec))
        catcoords = SkyCoord(ra*u.deg, dec*u.deg, frame=FK5)
        colra = [colname for colname in ['RA', 'ra', 'ALPHA_J2000'] if colname in refcat.colnames][0]
        coldec = [colname for colname in ['DEC', 'Dec', 'dec', 'DELTA_J2000'] if colname in refcat.colnames][0]
        colid = [colname for colname in ['ID', 'id', 'phot_id'] if colname in refcat.colnames][0]
        if colra != 'RA':
            refcat.rename_column(colra, 'RA')
        if coldec != 'DEC':
            refcat.rename_column(coldec, 'DEC')
        if colid != 'ID':
            refcat.rename_column(colid, 'ID')
        refra = refcat['RA'].data
        refdec = refcat['DEC'].data
        refid = refcat['ID'].data
        refcatcoords = SkyCoord(refra*u.deg, refdec*u.deg, frame=FK5)
        refx, refy = wcs.all_world2pix(refra, refdec, 1, quiet=True) # internally I use DS9 1-based pixcoord
        # match coordinate - find the matched source in refcat for each cat source.
        idx, d2d, d3d = match_coordinates_sky(catcoords, refcatcoords)
        matches = np.argwhere(np.logical_and(idx>=0, d2d<=maxsep)).ravel()
        # 20230726: to solve the issue "Number of output coordinates exceeded allocation"
        #           "Multiple sources within specified tolerance matched to a single reference source. "
        #           I need to make sure the cross-matching is one-to-one only.
        for kindex in idx:
            kduplicates = np.argwhere(idx==kindex).ravel() # indices into 'idx' 'd2d' and 'catcoords' arrays
            if len(kduplicates) > 1:
                kkeep = kduplicates[np.argmin(d2d[kduplicates])] # keep min 'd2d' cat source if there are duplicates
                for kk in kduplicates:
                    if kk != kkeep:
                        matches[kk] = False
        # make plot
        ax = axes[ipanel]
        ax.axis('on')
        if len(matches) == 0:
            logger.error('No match is found between the catfile of {!r} and the abs_refcat {!r}!'.format(imfile, abs_refcat))
            raise Exception('No match is found between the catfile of {!r} and the abs_refcat {!r}!'.format(imfile, abs_refcat))
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            d_ra = (ra[matches]-refra[idx[matches]]) * np.cos(np.deg2rad(refdec[idx[matches]])) * 3600.
            d_dec = (dec[matches]-refdec[idx[matches]]) * 3600.
            d_px = (x[matches]-refx[idx[matches]])
            d_py = (y[matches]-refy[idx[matches]])
            #lim = np.max(np.abs([d_ra, d_dec])) * 1.2 # axes range +20%
            # 
            if initial_offset is not None:
                d_ra += ioffra
                d_dec += ioffdec
                d_px += (-ioffra/pixsc)
                d_py += (ioffdec/pixsc)
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
            ax.scatter(d_ra, d_dec, s=3, alpha=0.9)
            #ax.set_xlim([lim, -lim])
            #ax.set_ylim([-lim, lim])
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
            #patch = Ellipse(xy=[d_ra_mean, d_dec_mean], width=2.*outlier_sigma*d_ra_sigma, height=2.*outlier_sigma*d_dec_sigma, 
            #                angle=0.0, ec='red', fc='none', alpha=0.8)
            patch = Ellipse(xy=[d_ra_median, d_dec_median], width=2.*outlier_sigma*d_ra_sigma, height=2.*outlier_sigma*d_dec_sigma, 
                            angle=0.0, ec='red', fc='none', alpha=0.8)
            ax.add_patch(patch)
            ax.text(d_ra_mean, d_dec_mean+outlier_sigma*d_dec_sigma, '{:+.4f}, {:+.4f}'.format(d_ra_mean, d_dec_mean),
                ha='center', va='bottom', color='red')
            # 
            #patch = Ellipse(xy=[d_ra_median, d_dec_median], width=2.*3.*d_ra_sigma, height=2.*3.*d_dec_sigma, 
            #                angle=0.0, ec='red', fc='none', ls='--', alpha=0.8)
            #ax.add_patch(patch)
            # 
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            lim = np.max(np.abs([xlim, ylim])) * 1.2 # axes range +20%
            ax.set_xlim([lim, -lim])
            ax.set_ylim([-lim, lim])
            #ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
            #ax.yaxis.set_major_locator(ticker.MaxNLocator(6))
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
            ax.plot([-lim, lim], [0., 0.], color='k', ls='dashed', lw=0.6, alpha=0.7)
            ax.plot([0., 0.], [-lim, lim], color='k', ls='dashed', lw=0.6, alpha=0.7)
        # 
        imfile_pieces = imfile.split('_')
        imfile_blocks = []
        imfile_block = ''
        while len(imfile_pieces) > 0:
            imfile_piece = imfile_pieces[0]
            if imfile_block == '':
                imfile_block = imfile_piece
                imfile_pieces.pop(0)
            elif len(imfile_block+'_'+imfile_piece) < 35:
                imfile_block = imfile_block+'_'+imfile_piece
                imfile_pieces.pop(0)
            else:
                imfile_blocks.append(imfile_block)
                imfile_block = ''
            if len(imfile_pieces) == 0:
                imfile_blocks.append(imfile_block)
                imfile_block = ''
        imfile_str = '\n_'.join(imfile_blocks)
        ax.set_title(imfile_str, fontsize='small')
        ax.set_aspect(1.0)
        # 
        icol = (ipanel%(ncols*nrows))%ncols
        irow = (ipanel%(ncols*nrows))//ncols
        ipage = ipanel//(ncols*nrows)
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
        if str(output_index_base) == '0':
            newcat['x'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
            newcat['y'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
        newcatpath = os.path.splitext(catpath)[0] + output_suffix + '.csv'
        if os.path.isfile(newcatpath):
            shutil.move(newcatpath, newcatpath+'.backup')
        newcat.write(newcatpath, format='csv')
        newcatpaths[imfile] = newcatpath
        logger.info('Output to {!r}'.format(newcatpath))
        # 
        # get matched refcat -- already got this above
        #chkcat = Table({'x': matched_abs_refcat_x[imfile], 'y': matched_abs_refcat_y[imfile]})
        # 
        # output matched refcat ds9 pix region
        chkcatpath = os.path.splitext(catpath)[0] + output_suffix + '_matched_refcat.ds9.pix.reg'
        if os.path.isfile(chkcatpath):
            shutil.move(chkcatpath, chkcatpath+'.backup')
        with open(chkcatpath, 'w') as fp:
            fp.write('# DS9 region file\n')
            fp.write('global color=yellow\n')
            fp.write('image\n')
            for k in range(len(chkcat)):
                fp.write('circle({}, {}, {}")\n'.format(chkcat['x'][k], chkcat['y'][k],0.35))
        logger.info('Output to {!r}'.format(chkcatpath))
        # 
        # output matched refcat
        if str(output_index_base) == '0':
            chkcat['x'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
            chkcat['y'] -= 1 # jwst tweakreg tweakwcs uses 0-indices, see "tweakwcs/correctors.py" `ra, dec = self._wcs.all_pix2world(x, y, 0)`
        chkcatpath = os.path.splitext(catpath)[0] + output_suffix + '_matched_refcat.csv'
        if os.path.isfile(chkcatpath):
            shutil.move(chkcatpath, chkcatpath+'.backup')
        chkcat.write(chkcatpath, format='csv', overwrite=True)
        logger.info('Output to {!r}'.format(chkcatpath))
    
    # 
    for ipage in range(npages):
        figs[ipage].tight_layout(w_pad=3.5)
    
    # write the big two-column cat_file
    cat_filebase, cat_filesuffix = os.path.splitext(cat_file)
    out_file = cat_filebase + output_suffix + cat_filesuffix
    out_json = cat_filebase + output_suffix + '.json'
    out_png = cat_filebase + output_suffix + '.png'
    out_pdf = cat_filebase + output_suffix + '.pdf'
    
    # save figure in png format if `len(imfiles) <= 10`
    if npages == 1:
        fig.savefig(out_png, dpi=300)
        logger.info('Output to {!r}'.format(out_png))
    
    # save figure in pdf format
    with PdfPages(out_pdf) as pdf:
        for fig in figs:
            plt.figure(fig.number)
            pdf.savefig()
    logger.info('Output to {!r}'.format(out_pdf))
    
    # save json
    if os.path.isfile(out_json):
        shutil.move(out_json, out_json+'.backup')
    with open(out_json, 'w') as fp:
        json.dump(offsets, fp, indent=4)
    logger.info('Output to {!r}'.format(out_json))
    
    # save output catalog
    if os.path.isfile(out_file):
        shutil.move(out_file, out_file+'.backup')
    with open(out_file, 'w') as fp:
        for imfile in cats:
            newcatpath = newcatpaths[imfile]
            fp.write('{} {}\n'.format(imfile, newcatpath))
    logger.info('Output to {!r}'.format(out_file))
    
    # save optimized abs_refcat
    if output_dir is None:
        output_dir = os.path.dirname(abs_refcat)
    if output_abs_refcat is None:
        output_abs_refcat = os.path.join(output_dir, 
            os.path.splitext(os.path.basename(abs_refcat))[0] + output_suffix + '.fits') # FITS-format catalog
    if output_abs_refcat != '':
        if 'name' not in refcat.meta: 
            refcat.meta['name'] = os.path.splitext(os.path.basename(abs_refcat))[0]
        refcat = refcat[matched_abs_refcat_mask]
        if os.path.isfile(output_abs_refcat):
            shutil.move(output_abs_refcat, output_abs_refcat+'.backup')
        refcat.write(output_abs_refcat)
        logger.info('Output to {!r}'.format(output_abs_refcat))
    
    # return
    return out_file, output_abs_refcat





@click.command()
@click.argument('cat_file', type=click.Path(exists=True))
@click.argument('abs_refcat', type=click.Path(exists=True))
@click.option('--search-radius', type=float, default=default_search_radius, help='coordinate matching max search radius in arcsec.')
@click.option('--output-suffix', type=str, default=default_output_suffix, help='output suffix.')
@click.option('--outlier-sigma', type=float, default=default_outlier_sigma, help='filter outlier by large offsets, in sigma of the offset distribution.')
@click.option('--initial-offset', type=float, nargs=2, default=default_initial_offset, help='initial offset (RA_offset, Dec_offset) in arcsec, increasing towards East and North, for all images in the input catfile.')
@click.option('--output-abs-refcat', type=click.Path(exists=False), default=default_output_abs_refcat, help='also output the filtered abs_refcat.')
@click.option('--input-index-base', type=click.Choice(['0', '1']), default=default_input_index_base, help='x y pixcoord base for input sextractor catalog listed in catfile. The default is 1 because sextractor uses 1-based pixcoord.')
@click.option('--output-index-base', type=click.Choice(['0', '1']), default=default_output_index_base, help='x y pixcoord base for output matched csv. The default is 0 because jwst pipeline tweakwcs correctors.py uses 0-based pixcoord.')
def main(
        cat_file,
        abs_refcat,
        search_radius,
        output_suffix,
        outlier_sigma, 
        initial_offset, 
        output_abs_refcat, 
        input_index_base, 
        output_index_base, 
    ):

    match_cat_file_to_abs_refcat_with_2dhist(
        cat_file, 
        abs_refcat, 
        search_radius, 
        output_suffix, 
        outlier_sigma, 
        initial_offset, 
        output_abs_refcat, 
        input_index_base, 
        output_index_base, 
    )





if __name__ == '__main__':
    main()










