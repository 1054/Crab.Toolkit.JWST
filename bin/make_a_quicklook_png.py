#!/usr/bin/env python
# 
import os, sys, re, shutil, glob, copy, datetime, time
import astropy.units as u
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy.io import fits
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, proj_plane_pixel_area
import warnings
warnings.filterwarnings('ignore')



@click.command()
@click.argument('img_file', type=click.Path(exists=True))
@click.argument('out_figure', type=click.Path(exists=False))
@click.option('--fov', '--fov_arcsec', 'fov_arcsec', type=float, default=5.0)
@click.option('--id', '--id-label', 'id_label', type=str, default=None)
@click.option('--img', '--img-label', 'img_label', type=str, default=None)
def make_cutout_png_figure(
        img_file, 
        out_figure, 
        fov_arcsec, 
        id_label, 
        img_label, 
    ):
    fig, ax = plt.subplots(
        ncols=1, nrows=1, figsize=(5, 5), 
        gridspec_kw=dict(left=0.01, right=0.99, bottom=0.01, top=0.99)
    )
    image, header = fits.getdata(img_file, header=True)
    if fov_arcsec is None:
        wcs = WCS(header, naxis=2)
        wcs.sip = None
        pixsc = proj_plane_pixel_scales(wcs) * 3600.0
        fov_arcsec = header['NAXIS2'] * pixsc[1]
    norm = simple_norm(image, 'asinh', min_percent = 0.1, max_percent = 99.9)
    ax.axis('off')
    ax.imshow(image, origin='lower', interpolation='nearest', norm=norm, cmap='gray', 
        extent=[-fov_arcsec/2.0, fov_arcsec/2.0, -fov_arcsec/2.0, fov_arcsec/2.0])
    ax.set_xlim([-fov_arcsec/2.0, fov_arcsec/2.0])
    ax.set_ylim([-fov_arcsec/2.0, fov_arcsec/2.0])
    if id_label is not None:
        ax.text(0.02, 0.99, id_label, ha='left', va='top', color='w', fontsize='large', transform=ax.transAxes)
    if img_label is not None:
        ax.text(0.02, 0.95, img_label, ha='left', va='top', color='w', fontsize='large', transform=ax.transAxes)
    ax.text(0.98, 0.01, 'FoV: {:.3g}'.format(fov_arcsec), ha='right', va='bottom', color='w', transform=ax.transAxes)
    ax.plot([-fov_arcsec/2.0 + 0.1*fov_arcsec/2.0, -fov_arcsec/2.0 + 0.1*fov_arcsec/2.0], 
            [-fov_arcsec/2.0 + 0.1*fov_arcsec/2.0, -fov_arcsec/2.0 + 0.1*fov_arcsec/2.0 + 1.0], 
            color='w')
    ax.text(-fov_arcsec/2.0 + 0.1*fov_arcsec/2.0 * 0.9, 
            -fov_arcsec/2.0 + 0.1*fov_arcsec/2.0 + 0.5,
            '1 arcsec', ha='center', va='bottom', 
            color='w', rotation=90, rotation_mode='anchor')
    out_figure = re.sub(r'^(.*)\.(pdf|png)$', r'\1', out_figure)
    fig.savefig(out_figure+'.pdf', dpi=300)
    fig.savefig(out_figure+'.png', dpi=300)
    print('Output to {!r}'.format(out_figure+'.png'))







# 
# main
# 
if __name__ == '__main__':
    
    make_cutout_png_figure()









