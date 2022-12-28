#!/usr/bin/env python
# 
"""
Comparing images from different data reduction by drawing random circles and analyzing the statistics.

Usage: 
    ./util_compare_image_reduction.py image_1.fits image_2.fits image_3.fits output_name

Input: 
    FITS images

Output:
    Analysis files: 
        {output_name}_flux_vs_flux.pdf
        {output_name}_pos_vs_pos.pdf

Example:
    ```
    ./util_compare_image_reduction.py \
        ~/jwst_reduced_data_ceers_by_ceers/hlsp_ceers_jwst_nircam_nircam1_f444w_dr0.5_i2d.fits \
        ~/jwst_reduced_data_ceers_by_daizhong/jw01345_obs001_NIRCAM_F444W_i2d.fits \
        ~/jwst_reduced_data_ceers_by_grizli/ceers-ne-grizli-v4.0-f444w-clear_drc_sci.fits \
        ~/jwst_reduced_data_ceers_by_max/mosaics/CEERS1/mosaic_nircam_f444w_CEERS_without_snowballs_i2d.fits \
        --aperture-size 0.2 \
        --aperture-number 10000 \
        --input-labels 'ceers,daizhong,grizli,max' \
        --plot-orders '0,3,1,2' \
        --aper-corr 0.606 \
        with_aperture_size_0.2_arcsec/output_aperture_statistics_ceers_obs001_f444w
    ```

By Daizhong Liu @MPE. 

Last updates: 
    2022-12-14 
    2022-12-19 

"""
import os, sys, re, shutil
import astropy.units as u
import astropy.constants as const
import click
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, proj_plane_pixel_area
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['axes.labelsize'] = '15' # https://matplotlib.org/users/customizing.html
mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['axes.labelpad'] = 6
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['xtick.labelsize'] = 11
mpl.rcParams['ytick.labelsize'] = 11
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
#mpl.rcParams['grid.color'] = 'b0b0b0'
mpl.rcParams['grid.linestyle'] = '--'
mpl.rcParams['grid.linewidth'] = 0.25
mpl.rcParams['grid.alpha'] = 0.8
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['legend.borderaxespad'] = 0.2 # space between legend border and axis
mpl.rcParams['legend.handletextpad'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.1
#mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rcParams['mathtext.fontset'] = 'dejavuserif'
mpl.rcParams['font.size'] = 12

# code name and version
CODE_NAME = 'util_compare_image_reduction.py'
CODE_AUTHOR = 'Daizhong Liu'
CODE_VERSION = '20221214'
CODE_HOMEPAGE = ''

# logging
import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
logger = logging.getLogger(CODE_NAME)
logger.setLevel(logging.DEBUG)


# defaults
default_aperture_size = 1.5 # arcsec
default_aperture_number = 1000 


# define class
class ImageComparer(object):
    """docstring for ImageComparer"""
    def __init__(self, input_images):
        self.image_files = input_images
        self.hdulists = []
        self.imhdus = []
        self.extnames = []
        self.wcses = []
        self.pixscs = []
        self.rects = []
        for image_file in self.image_files:
            hdul = fits.open(image_file)
            if 'SCI' in hdul:
                extname = 'SCI'
            else:
                extname = ''
                for ihdu in range(len(hdul)):
                    if hdul[ihdu].header['NAXIS'] == 2:
                        extname = hdul[ihdu].name
                        break
                if extname == '':
                    raise Exception('Error! Could not find 2D image in fits files: {!r}'.format(image_file))
            logger.info('Reading {!r} extension {!r}'.format(image_file, extname))
            imhdu = hdul[extname]
            header = hdul[extname].header
            wcs = WCS(header, naxis=2)
            pixsc = np.sqrt(proj_plane_pixel_area(wcs))*3600.
            ny, nx = header['NAXIS2'], header['NAXIS1']
            rect = wcs.wcs_pix2world([1, nx, nx, 1], [1, 1, ny, ny], 1)
            self.hdulists.append(hdul)
            self.imhdus.append(imhdu)
            self.extnames.append(extname)
            self.wcses.append(wcs)
            self.pixscs.append(pixsc)
            self.rects.extend(np.array(rect).T.tolist()) # make rect = [ (x1, y1), (x2, y2), ... ]
        self.rects = np.array(self.rects)
        self.aperture_size = None
        self.apertures_table = None
    # 
    def __enter__(self):
        return self
    # 
    def __exit__(self, type, value, traceback):
        for hdul in self.hdulists:
            hdul.close()
    # 
    def ensure_directory(self,
            output_path,
        ): 
        temp_dir = os.path.dirname(output_path)
        if temp_dir != '':
            if not os.path.isdir(temp_dir):
                os.makedirs(temp_dir)
    # 
    def draw_apertures(self, 
            output_name, 
            aperture_size = default_aperture_size, 
            aperture_number = default_aperture_number, 
        ): 
        self.ensure_directory(output_name)
        apertures_dict = OrderedDict()
        apertures_dict['id'] = []
        apertures_dict['ra'] = []
        apertures_dict['dec'] = []
        for i in range(len(self.imhdus)):
            apertures_dict[f'x_{i}'] = []
            apertures_dict[f'y_{i}'] = []
            apertures_dict[f'rad_{i}'] = []
            apertures_dict[f'pixsc_{i}'] = []
            apertures_dict[f'ix_{i}'] = []
            apertures_dict[f'iy_{i}'] = []
            apertures_dict[f'irad_{i}'] = []
            apertures_dict[f'fluxconv_{i}'] = []
            apertures_dict[f'sum_{i}'] = []
            apertures_dict[f'mean_{i}'] = []
            apertures_dict[f'median_{i}'] = []
            apertures_dict[f'stddev_{i}'] = []
            apertures_dict[f'weighted_x_{i}'] = []
            apertures_dict[f'weighted_y_{i}'] = []
            apertures_dict[f'weighted_ra_{i}'] = []
            apertures_dict[f'weighted_dec_{i}'] = []
        chunk = 100
        count = 0
        # calc boundary
        max_ra = np.min([self.rects[0, 0], self.rects[3, 0]])
        min_ra = np.max([self.rects[1, 0], self.rects[2, 0]])
        min_dec = np.max([self.rects[0, 1], self.rects[1, 1]])
        max_dec = np.min([self.rects[2, 1], self.rects[3, 1]])
        while count < aperture_number:
            logger.info('Processing aperture {}'.format(count))
            rand_ra = np.random.random_sample()
            rand_dec = np.random.random_sample()
            ra = (max_ra - min_ra) * rand_ra + min_ra
            dec = (max_dec - min_dec) * rand_dec + min_dec
            ixyrad_list = []
            fxyrad_list = []
            for i in range(len(self.imhdus)):
                imhdu = self.imhdus[i]
                image = imhdu.data
                header = imhdu.header
                wcs = self.wcses[i]
                pixsc = self.pixscs[i]
                (fx,), (fy,) = wcs.wcs_world2pix([ra], [dec], 0)
                ix, iy = int(np.round(fx)), int(np.round(fy))
                ny, nx = header['NAXIS2'], header['NAXIS1']
                frad = aperture_size/2.0/pixsc
                irad = int(np.ceil(frad))
                if ix > irad and ix < nx-irad-1 and iy > irad and iy < ny-irad-1:
                    if np.all(np.isfinite(image[iy-irad:iy+irad+1, ix-irad:ix+irad+1])): 
                        if np.all(~np.isclose(image[iy-irad:iy+irad+1, ix-irad:ix+irad+1], 0)):
                            ixyrad_list.append([ix, iy, irad])
                            fxyrad_list.append([fx, fy, frad])
            # need all image to have valid pixels
            if len(ixyrad_list) != len(self.imhdus):
                continue
            # draw aperture
            apertures_dict['id'].append(count)
            apertures_dict['ra'].append(ra)
            apertures_dict['dec'].append(dec)
            for i in range(len(self.imhdus)):
                imhdu = self.imhdus[i]
                #image = imhdu.data
                header = imhdu.header
                wcs = self.wcses[i]
                pixsc = self.pixscs[i]
                ny, nx = header['NAXIS2'], header['NAXIS1']
                ix, iy, irad = ixyrad_list[i]
                fx, fy, frad = fxyrad_list[i]
                image = imhdu.data[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                #ggy, ggx = np.mgrid[0:ny, 0:nx]
                #gx = ggx[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                #gy = ggy[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                gy, gx = np.mgrid[iy-irad:iy+irad+1, ix-irad:ix+irad+1]
                grad = np.sqrt((gx-fx)**2 + (gy-fy)**2)
                mask = (grad < frad)
                if header['BUNIT'] == 'MJy/sr':
                    fluxconv = 1
                elif header['BUNIT'] == '10.0*nanoJansky':
                    fluxconv = (10.0*u.nJy / (pixsc*u.arcsec)**2).to(u.MJy/u.sr).value
                else:
                    raise Exception('Error! BUNIT {!r} is not implemented!'.format(header['BUNIT']))
                aperdata = image[mask] * fluxconv
                fraction = 0.5 - (grad - frad)
                fraction[fraction<0.0] = 0.0
                fraction[fraction>1.0] = 1.0
                fluxsum = np.nansum(image * fraction * fluxconv)
                fmask = mask.astype(int).astype(float)
                #print(gx.shape, image.shape, fmask.shape)
                weighted_x = np.nansum(gx*image*fmask) / np.nansum(image*fmask)
                weighted_y = np.nansum(gy*image*fmask) / np.nansum(image*fmask)
                (weighted_ra,), (weighted_dec,) = wcs.wcs_pix2world([weighted_x], [weighted_y], 0)
                apertures_dict[f'x_{i}'].append(fx)
                apertures_dict[f'y_{i}'].append(fy)
                apertures_dict[f'rad_{i}'].append(frad)
                apertures_dict[f'pixsc_{i}'].append(pixsc)
                apertures_dict[f'ix_{i}'].append(ix)
                apertures_dict[f'iy_{i}'].append(iy)
                apertures_dict[f'irad_{i}'].append(irad)
                apertures_dict[f'fluxconv_{i}'].append(fluxconv)
                apertures_dict[f'sum_{i}'].append(fluxsum)
                apertures_dict[f'mean_{i}'].append(np.nanmean(aperdata))
                apertures_dict[f'median_{i}'].append(np.nanmedian(aperdata))
                apertures_dict[f'stddev_{i}'].append(np.nanstd(aperdata))
                apertures_dict[f'weighted_x_{i}'].append(weighted_x)
                apertures_dict[f'weighted_y_{i}'].append(weighted_y)
                apertures_dict[f'weighted_ra_{i}'].append(weighted_ra)
                apertures_dict[f'weighted_dec_{i}'].append(weighted_dec)
            # 
            # dump table
            #if count > 0 and aperture_number > 10 and count % int(aperture_number/10) == 0 :
            #    apertures_table = Table(apertures_dict)
            #    apertures_table.write(f'{output_name}.dump.{count}.csv', overwrite=True)
            #    logger.info('Output to {!r}'.format(f'{output_name}.dump.{count}.csv'))
            # 
            # next count
            count += 1
        # 
        self.aperture_size = aperture_size
        self.apertures_table = Table(apertures_dict)
        self.apertures_table.write(output_name+'.csv', overwrite=True)
        logger.info('Output to {!r}'.format(output_name+'.csv'))
    # 
    def load_apertures(self,
            input_name,
        ):
        # 
        if not os.path.isfile(input_name) and os.path.isfile(input_name+'.csv'):
            input_name += '.csv'
        logger.info('Loading {!r}'.format(input_name))
        self.apertures_table = Table.read(input_name)
        self.aperture_size = self.apertures_table['rad_0'][0] * self.apertures_table['pixsc_0'][0] * 2.0
    # 
    def plot_histogram(self,
            labels, # a list corresponding to _0 _1 _2 _3 columns in the table
            orders, # a list, e.g., [0, 1, 2, 3]
            colors, # a list, e.g., ['red', 'blue', 'orange', 'cyan']
            output_figure, 
            aperture_psf_fraction = 1.0, 
            title = '', 
        ):
        # 
        if self.apertures_table is not None:
            logger.info('Plotting histogram and depth')
            tb = self.apertures_table
            aperture_size = self.aperture_size
            ncols = 1
            nrows = len(labels)
            panel_width = 4.5 # inches
            panel_height = 2.8 # inches
            fig_left = 0.7 # inches
            fig_right = 0.2 # inches
            fig_bottom = 0.6 # inches
            fig_top = 0.4 # inches
            panel_wspace = 0.2 # inches
            panel_hspace = 0.25 # inches
            fig_width = fig_left + panel_width * ncols + panel_wspace * (ncols-1) + fig_right
            fig_height = fig_bottom + panel_height * nrows + panel_hspace * (nrows-1) + fig_top
            logger.debug('fig_width, fig_height: {}, {}'.format(fig_width, fig_height))
            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width, fig_height), 
                gridspec_kw=dict(left=fig_left/fig_width, right=1.-fig_right/fig_width, 
                    bottom=fig_bottom/fig_height, top=1.-fig_top/fig_height, 
                    wspace=panel_wspace/panel_width,
                    hspace=panel_hspace/panel_height)
            )
            # 
            fluxes_list = []
            pixsc_list = []
            min_flux = None
            max_flux = None
            p16_flux = None
            p84_flux = None
            orders = [int(t) for t in orders]
            for i in range(len(labels)):
                # 
                ii = orders.index(i)
                fluxes = tb[f'sum_{ii}']
                pixsc = tb[f'pixsc_{ii}'][0]
                pixsc_list.append(pixsc)
                fluxes = (fluxes*u.MJy/u.sr * (pixsc*u.arcsec)**2).to(u.nJy)
                fluxes_list.append(fluxes)
                if min_flux is None:
                    min_flux = np.min(fluxes)
                    max_flux = np.max(fluxes)
                    p16_flux = np.nanpercentile(fluxes, 16)
                    p84_flux = np.nanpercentile(fluxes, 84)
                else:
                    min_flux = min(min_flux, np.min(fluxes))
                    max_flux = max(max_flux, np.max(fluxes))
                    p16_flux = min(p16_flux, np.nanpercentile(fluxes, 16))
                    p84_flux = max(p84_flux, np.nanpercentile(fluxes, 84))
            # 
            logger.debug('min_flux, max_flux: {}, {}'.format(min_flux, max_flux))
            logger.debug('p16_flux, p84_flux: {}, {}'.format(p16_flux, p84_flux))
            init_sigma = (p84_flux-p16_flux) / 2.0
            bin_min = (p84_flux+p16_flux)/2.0 - init_sigma * 3.0
            bin_max = (p84_flux+p16_flux)/2.0 + init_sigma * 5.0
            logger.debug('bin_min, bin_max: {}, {}'.format(bin_min, bin_max))
            if isinstance(bin_min, u.Quantity):
                bin_min = bin_min.value
            if isinstance(bin_max, u.Quantity):
                bin_max = bin_max.value
            nbins = 151
            # 
            for i in range(len(labels)):
                # 
                ii = orders.index(i)
                ax = axes[i]
                # 
                fluxes = fluxes_list[i]
                pixsc = pixsc_list[i]
                # 
                hist, bin_edges = np.histogram(fluxes, bins=nbins, range=(bin_min, bin_max))
                # 
                x = (bin_edges[0:-1]+bin_edges[1:])/2.0
                y = hist
                #y = hist/np.nanmax(hist)
                ax.step(x, y, where='mid', 
                    label=labels[ii], color=colors[i], alpha=0.8,
                    zorder=50-orders[i],
                    )
                # 
                ipeak = np.argmax(hist)
                ypeak = y[ipeak]
                xpeak = x[ipeak]
                ihalf = np.argwhere(np.logical_and(x<xpeak, y>0.5*ypeak)).ravel()[0]
                xhalf = x[ihalf]
                fitmask = (x < xpeak+(xpeak-xhalf))
                logger.debug('label: {!r}, xpeak: {}, xhalf: {}, fitmask: x < {}, pixsc: {}'.format(
                    labels[ii], xpeak, xhalf, xpeak+(xpeak-xhalf), pixsc))
                # 
                mod = models.Gaussian1D(amplitude=ypeak, mean=0.0, stddev=init_sigma)
                fit = fitting.LevMarLSQFitter()
                fitted_line = fit(mod, x[fitmask], y[fitmask])
                ax.plot(x, fitted_line(x), color=colors[i], ls='dashed', alpha=0.8,
                    zorder=50-orders[i],
                    )
                # 
                mean = fitted_line.mean.value
                sigma = fitted_line.stddev.value
                mean_mag = (mean * u.nJy).to(u.ABmag).value
                sigma_mag = (sigma * u.nJy / aperture_psf_fraction).to(u.ABmag).value
                five_sigma_mag = (5.0*sigma * u.nJy / aperture_psf_fraction).to(u.ABmag).value
                # 
                ax.text(0.99, 0.88, 
                    'histogram: {}'.format(labels[ii]) + '\n' + \
                    'aperture diameter: {:.2f} [arcsec]'.format(aperture_size) + '\n' + \
                    'pixel size: {:.4f} [arcsec]'.format(pixsc) + '\n' + \
                    'mean: {:.5g} [nanoJy]'.format(mean) + '\n' + \
                    'sigma: {:.5g} [nanoJy]'.format(sigma) + '\n' + \
                    'aper corr: {:.3g}'.format(aperture_psf_fraction) + '\n' + \
                    '1-sigma: {:.3f} magAB'.format(sigma_mag) + '\n' + \
                    '5-sigma: {:.3f} magAB'.format(five_sigma_mag),
                    ha='right', va='top', fontsize=9,
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'),
                    transform=ax.transAxes)
                # 
                #ax.set_xlim(ax.get_xlim()[::-1])
                #ax.set_yscale('log')
                if i == len(labels)-1:
                    ax.set_xlabel('flux in aperture [nanoJy]')
                ax.set_ylabel('N')
                ax.legend()
            # 
            # title
            if title != '':
                axes[0].set_title(title)
            # 
            if output_figure.endswith('.pdf') or output_figure.endswith('.png'):
                output_figure = os.path.splitext(output_figure)
            self.ensure_directory(output_figure)
            fig.savefig(output_figure+'.pdf')
            fig.savefig(output_figure+'.png')
            print('Output to {!r}'.format(output_figure+'.*'))

            
        



# 
@click.command()
@click.argument('input_images', type=click.Path(exists=True), nargs=-1)
@click.argument('output_name', type=str, nargs=1)
@click.option('--aperture-size', type=float, default=default_aperture_size, help='aperture diameter in arcsec.')
@click.option('--aperture-number', type=int, default=default_aperture_number, help='aperture number.')
@click.option('--input-labels', type=str, default=None, help='a string for the labels of the input images, using comma to separate.')
@click.option('--plot-orders', type=str, default=None, help='orders for plotting, from 0 to N-1, using comma to separate.')
@click.option('--plot-colors', type=str, default=None, help='colors for plotting, using comma to separate.')
@click.option('--aper-corr', type=float, default=1.0, help='aperture to psf total flux fraction, 0-1.')
@click.option('--title', type=str, default='', help='title of this plot.')
@click.option('--output-dir', type=click.Path(exists=False), default=None, help='not implemented, just use path in `output_name`.')
@click.option('--overwrite/--no-overwrite', is_flag=True, default=False)
@click.option('--verbose/--no-verbose', is_flag=True, default=True)
def main(
        input_images, 
        output_name, 
        aperture_size, 
        aperture_number,
        input_labels,
        plot_orders,
        plot_colors,
        aper_corr,
        title,
        output_dir,
        overwrite,
        verbose, 
    ):
    
    with ImageComparer(input_images) as imcmp:
        if os.path.isfile(output_name+'.csv') and not overwrite:
            imcmp.load_apertures(output_name)
        else:
            imcmp.draw_apertures(output_name, aperture_size=aperture_size, aperture_number=aperture_number)
        # 
        if input_labels is not None:
            if plot_orders is None:
                plot_orders = ','.join(np.arange(len(input_labels)).astype(str))
            if plot_colors is None:
                plot_colors = ','.join(['red', 'blue', 'orange', 'cyan'][0:len(input_labels)])
            imcmp.plot_histogram(
                input_labels.split(','), 
                plot_orders.split(','), 
                plot_colors.split(','), 
                output_name+'_flux_histograms',
                aperture_psf_fraction = aper_corr,
                title = title, 
            )
    



# 
if __name__ == '__main__':
    
    main()



