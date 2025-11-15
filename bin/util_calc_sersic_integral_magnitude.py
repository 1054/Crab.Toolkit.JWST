#!/usr/bin/env python

import astropy.units as u
import numpy as np
import json
from astropy.table import Table
from scipy.special import gamma, gammaincinv
pixsc = 0.030

statmorph_results = Table.read('COSMOS_Web_cutouts_run_statmorph_results_with_radec.csv')
n = statmorph_results["sersic_n"]
r_eff = statmorph_results["sersic_rhalf"]
I_eff = statmorph_results["sersic_amplitude"]

# 
bn = gammaincinv(2.*n, 0.5) # =1.678 for n=1, see https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
alpha = r_eff / (bn**n)
integ = I_eff * 2.0*np.pi * alpha**2 * n * gamma(2*n) * np.exp(bn)

#print('integ: {}'.format(integ))

mag = ((integ * u.MJy/u.sr) * (((pixsc * u.arcsec)**2).to(u.sr))).to(u.ABmag)

#print('mag: {}'.format(mag))

#print('statmorph, nsersic, re, mag277w: {:.4f}, {:.4f} (pixels) = {:.4f} (arcsec), {:.4f}'.format(n, r_eff, r_eff * pixsc, mag))

statmorph_results["sersic_totalflux"] = integ
statmorph_results["sersic_mag"] = mag
statmorph_results.write('COSMOS_Web_cutouts_run_statmorph_results_with_radec_with_mag.csv')
