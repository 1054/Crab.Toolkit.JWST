#!/usr/bin/env python
# 
import os, sys, re, copy, shutil
import astropy.units as u
import click
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits # fits.ImageHDU
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.visualization import simple_norm, ZScaleInterval

import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = '9' # https://matplotlib.org/users/customizing.html
mpl.rcParams['axes.grid'] = False
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['axes.labelpad'] = 10
mpl.rcParams['font.size'] = 10
mpl.rcParams['font.family'] = 'Monaco'
#mpl.rcParams['xtick.direction'] = 'in'
#mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['xtick.minor.visible'] = True
#mpl.rcParams['ytick.minor.visible'] = True
#mpl.rcParams['xtick.labelsize'] = '9'
#mpl.rcParams['ytick.labelsize'] = '9'
#mpl.rcParams['xtick.top'] = True
#mpl.rcParams['ytick.right'] = True
mpl.rcParams['grid.color'] = 'white'
mpl.rcParams['grid.linewidth'] = 0.8
mpl.rcParams['grid.alpha'] = 0.8
mpl.rcParams['legend.fontsize'] = '10'
mpl.rcParams['legend.borderaxespad'] = 0.2 # space between legend border and axis
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['savefig.dpi'] = 300
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
import matplotlib.cm as mplcm
#from matplotlib.cm import register_cmap, cmap_d
import matplotlib.colors as mplcolors
from matplotlib.patches import Ellipse, Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable
#plt.style.use('dark_background')
import warnings
warnings.filterwarnings('ignore')


"""
This script will create a quick-look image in png format.
"""


@click.command()
@click.argument('input_image_file', type=click.Path(exists=True))
@click.argument('output_image_file', type=click.Path(exists=False), required=False, default=None)
@click.option('--figsize', type=float, default=6.0, help='Figure main panel size in inches, not including margins.')
def main(
        input_image_file, 
        output_image_file, 
        figsize, 
    ):
    
    # read input fits image
    with fits.open(input_image_file) as hdulist:
        ihdu = 0
        while ihdu < len(hdulist):
            #print(type(hdulist[ihdu]))
            if isinstance(hdulist[ihdu], fits.hdu.image.PrimaryHDU):
                main_header = copy.copy(hdulist[ihdu].header)
            if isinstance(hdulist[ihdu], fits.ImageHDU) or isinstance(hdulist[ihdu], fits.hdu.image.PrimaryHDU):
                if hdulist[ihdu].header['NAXIS'] > 0:
                    hdu = copy.copy(hdulist[ihdu])
                    image = copy.copy(hdu.data)
                    header = copy.copy(hdu.header)
                    print('Read ihdu {} {}'.format(ihdu, type(hdulist[ihdu])))
                    break
            ihdu += 1

    # determine wcs, pixscale, x_size y_size
    wcs = WCS(header, naxis=2)
    pixscales = proj_plane_pixel_scales(wcs) * 3600.0
    print('pixscales: {}'.format(pixscales))
    x_pixsc = np.abs(pixscales[0])
    y_pixsc = np.abs(pixscales[1])
    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    print('x_size: {}, y_size: {}'.format(x_size, y_size))
    #print('wcs.wcs', wcs.wcs)
    if hasattr(wcs.wcs, 'pc'):
        PA = np.rad2deg(np.arctan2(wcs.wcs.pc[1][0], wcs.wcs.pc[1][1]))
        print('PC2_1: {}, PC2_2: {}, PA: {}'.format(wcs.wcs.pc[1][0], wcs.wcs.pc[1][1], PA))
    elif hasattr(wcs.wcs, 'cd'):
        PA = np.rad2deg(np.arctan2(wcs.wcs.cd[1][0], wcs.wcs.cd[1][1]))
        print('CD2_1: {}, CD2_2: {}, PA: {}'.format(wcs.wcs.cd[1][0], wcs.wcs.cd[1][1], PA))
    elif hasattr(wcs.wcs, 'crota'):
        PA = wcs.wcs.crota[1]
        print('CROTA1: {}, CROTA2: {}, PA: {}'.format(wcs.wcs.crota[0], wcs.wcs.crota[1], PA))
    else:
        raise Exception('The input fits header\'s WCS contain no pc, cd or crota! (checked `wcs.wcs.pc`, `wcs.wcs.cd` and `wcs.wcs.crota`)')
    
    # create image
    fig_margin_left = 1.0 # inch
    fig_margin_right = 0.6 # inch
    fig_margin_bottom = 0.8 # inch
    fig_margin_top = 0.6 # inch
    ncols = 1
    nrows = 1
    panel_aspect_ratio = float(y_size)/float(x_size)
    if panel_aspect_ratio > 1.0:
        panel_width = figsize # inch
        panel_height = panel_width * panel_aspect_ratio # inch
    else:
        panel_height = figsize # inch
        panel_width = panel_height / panel_aspect_ratio # inch
    print('panel_width: {}, panel_height: {}'.format(panel_width, panel_height))
    panel_widths = [panel_width]*ncols
    panel_heights = [panel_height]*nrows
    panel_wspace = 0.1 # inch
    panel_hspace = 0.1 # inch
    fig_width = fig_margin_left + np.sum(panel_widths) + panel_wspace * (ncols-1) + fig_margin_right
    fig_height = fig_margin_bottom + np.sum(panel_heights) * nrows + panel_hspace * (nrows-1) + fig_margin_top
    gridspec_kw = {}
    gridspec_kw['wspace'] = panel_wspace / np.mean(panel_widths) # matplotlib uses the fraction of horizontal spacing to average panel width as the wspace
    gridspec_kw['hspace'] = panel_hspace / np.mean(panel_heights) # matplotlib uses the fraction of vertical spacing to average panel height as the hspace
    gridspec_kw['top'] = 1.0 - (fig_margin_top / fig_height)
    gridspec_kw['bottom'] = fig_margin_bottom / fig_height
    gridspec_kw['right'] = 1.0 - (fig_margin_right / fig_width)
    gridspec_kw['left'] = fig_margin_left / fig_width
    print('gridspec_kw: {}'.format(gridspec_kw))
    fig, ax = plt.subplots(
        ncols=ncols, nrows=nrows, figsize=(fig_width, fig_height), 
        gridspec_kw=gridspec_kw,
        subplot_kw=dict(projection=wcs),
    )
    
    # plot imshow
    zscale_interval = ZScaleInterval()
    min_cut, max_cut = zscale_interval.get_limits(image)
    norm = simple_norm(image, 'linear', min_cut=min_cut, max_cut=max_cut)
    ax.imshow(
        image, origin='lower', interpolation='nearest', 
        norm=norm, cmap='gray', aspect=1,
    )
    
    # plot wcs grid
    ax.coords.grid(True, color='white', ls='dotted')
    ax.coords[0].set_auto_axislabel(False)
    ax.coords[1].set_auto_axislabel(False)
    ax.coords[0].set_axislabel('Right Ascension (J2000)')
    ax.coords[1].set_axislabel('Declination (J2000)')
    ax.coords[0].set_ticks(number=10, color='k', direction='out')
    ax.coords[1].set_ticks(number=10, color='k', direction='out')
    ax.coords[0].set_ticklabel(rotation=-PA-90, bbox=dict(color='w', alpha=0.7, boxstyle='round'), rotation_mode='anchor', pad=0)
    ax.coords[1].set_ticklabel(rotation=-PA, bbox=dict(color='w', alpha=0.7, boxstyle='round'), rotation_mode='anchor', pad=0)
    ax.coords[0].set_major_formatter('hh:mm:ss.s')
    ax.coords[1].set_major_formatter('dd:mm:ss.s')
    ax.coords[0].set_separator(':')
    ax.coords[1].set_separator(':')
    ax.coords[0].ticklabels.simplify_labels = lambda : None
    ax.coords[1].ticklabels.simplify_labels = lambda : None
    #print('ax.coords[0].ticklabels', ax.coords[0].ticklabels)
    
    overlay = ax.get_coords_overlay('fk5')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_auto_axislabel(False)
    overlay[1].set_auto_axislabel(False)
    overlay[0].set_axislabel('')
    overlay[1].set_axislabel('')
    overlay[0].set_ticks(number=10, color='k', direction='out')
    overlay[1].set_ticks(number=10, color='k', direction='out')
    overlay[0].set_ticklabel(rotation=-PA-90, bbox=dict(color='w', alpha=0.7, boxstyle='round'), rotation_mode='anchor', pad=0)
    overlay[1].set_ticklabel(rotation=-PA, bbox=dict(color='w', alpha=0.7, boxstyle='round'), rotation_mode='anchor', pad=0)
    overlay[0].set_major_formatter('hh:mm:ss.s')
    overlay[1].set_major_formatter('dd:mm:ss.s')
    overlay[0].set_separator(':')
    overlay[1].set_separator(':')
    overlay[0].ticklabels.simplify_labels = lambda : None
    overlay[1].ticklabels.simplify_labels = lambda : None
    
    # draw
    fig.canvas.draw()
    
    # # align wcs grid ticklabels
    # #print('ax.coords[0].ticklabels.ha', ax.coords[0].ticklabels.ha) # by looking at astropy.visualization.wcsaxes.ticklabels.py
    # #print(type(overlay[0].ticklabels), overlay[0].ticklabels.__dict__)
    # for iaxis in ax.coords[0].ticklabels.ha:
    #     for ibox in ax.coords[0].ticklabels.ha[iaxis]:
    #         ax.coords[0].ticklabels.ha[iaxis][ibox] = 'right'
    #         ax.coords[0].ticklabels.va[iaxis][ibox] = 'center'
    
    # # draw
    # fig.canvas.flush_events()
    
    # prepare output directory
    if output_image_file is None:
        output_image_file = os.path.splitext(input_image_file)[0] + '.quicklook.png'
    if os.path.isfile(output_image_file):
        shutil.copy2(output_image_file, output_image_file+'.backup')
    else:
        output_dir = os.path.dirname(output_image_file)
        if output_dir != '' and not os.path.isdir(output_dir):
            os.makedirs(output_dir)
    
    # write output file
    fig.savefig(output_image_file)
    print('Created quicklook image "{}"'.format(output_image_file))



if __name__ == '__main__':
    
    main()



