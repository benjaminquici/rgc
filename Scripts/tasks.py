import numpy as np

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table,Column

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt

from astropy.visualization import simple_norm

from reproject import reproject_interp

import matplotlib as mpl

import numpy.ma as ma

import matplotlib.pyplot as plt

# mpl.use('TkAgg')

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
Popen = subprocess.Popen


def get_glass_region(ra,dec):
    if ra < 351 and ra > 347 and dec > -35 and dec < -32.5 :
        return('A')
    elif ra < 351 and ra > 347 and dec > -32.5 and dec < -30 :
        return('B')
    elif ra < 347 and ra > 343 and dec > -32.5 and dec < -30 :
        return('C')
    elif ra < 347 and ra > 343 and dec > -35 and dec < -32.5 :
        return('D')
    elif ra < 343 and ra > 339 and dec > -32.5 and dec < -30 :
        return('E')
    elif ra < 343 and ra > 339 and dec > -35 and dec < -32.5 :
        return('F')
    else:
        return('No region found')

def get_emu_half(ra):
    if ra <= 345:
        return(132)
    elif ra > 345:
        return(137)

def make_box(ra, dec, a, ralist, declist, col=False):
    decmin,decmax,ramin,ramax = dec-a,dec+a,ra-a,ra+a
    ra_range=ma.masked_outside(ralist,ramin,ramax)
    dec_range=ma.masked_outside(declist,decmin,decmax)
    declist=dec_range[np.where((ra_range.mask==False)&(dec_range.mask==False))]
    ralist=ra_range[np.where((ra_range.mask==False)&(dec_range.mask==False))]
    if col:
        col=col[np.where((ra_range.mask==False)&(dec_range.mask==False))]
        return(ralist, declist, col)
    else:
        return(ralist, declist)

# # Popen("echo Hello")
# plt.scatter([1,2,5],[2,3,3])
# plt.show()