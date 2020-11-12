#!/usr/bin/python3

from tasks import get_glass_region
from tasks import get_emu_half
from tasks import make_box

import numpy as np
import argparse

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table, Column

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt

from astropy.visualization import simple_norm

from reproject import reproject_interp

from create_source_new import generate_cutout

import matplotlib 
from matplotlib.colors import Normalize

def plot_contour(ax, data, levels, lw, color, alpha):
    ax.contour(data,levels=levels,linewidths=lw,colors=color,alpha=alpha)
    return(ax)

def hide_ax(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_facecolor('white')
    ax.patch.set_alpha(0.75)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    return(ax)

def viewer(source_name, imsize, reproject=False):

    path_to_source = '/mnt/e/Thesis/Paper_II/G23_sample_first_pass/{0}/'.format(source_name)
    images = ['MIDAS', 'GMRT', 'EMU', 'GLASS', 'VIKING']
    freq = ['216 MHz', '399 MHz', '887 MHz', '5.5 GHz', '2.12 $\\mu$m']
    image_type = ['r', 'r', 'r', 'r', 'o']
    contour = [True, False, True, True, False]
    n_contour_levs = [1, 3, 3, 3, 0]
    image_rms = [9e-4, 3e-5, 45e-6, 35e-6, 0]
    max_contour_lev = [3, 50, 100, 100, 0]
    contour_colors = ['white', 'green', 'cyan', 'magenta', 'None']
    contour_lw = [1, 1, 1, 1, 1]
    contour_alpha = [0.5, 1, 0.75, 1, 1]

    hdu_all = []
    imdata_all = []
    header_all = []
    w_all = []
    data_name_all = []
    for i in range(0,len(images)):
        hdu = fits.open('{0}/{1}_{2}_{3}.fits'.format(path_to_source, source_name, images[i], imsize))
        hdu_all.append(hdu)
        header_all.append(hdu[0].header)
        imdata_all.append(hdu[0].data)
        w_all.append(wcs.WCS(hdu[0].header, naxis=2))
        data_name_all.append('{} {}'.format(images[i], freq[i]))

    horizontal_pad = [0.08, 0.38, 0.68, 0.08, 0.38]
    vertical_pad = [0.53, 0.53, 0.53, 0.09, 0.09]
    horizontal_pad0 = [0.1, 0.4, 0.7, 0.1, 0.4]
    vertical_pad0 = [0.905, 0.905, 0.905, 0.46, 0.46]
    height, width = 0.43, 0.28

    latleft = [True, False, False, True, False]
    latright = [False, False, False, False, False]
    lonbottom = [False, False, False, True, True]
    lontop = [False, False, False, False, False]

    fig = plt.figure(figsize=(17,10))
    
    for i in range(0, len(images)):
        lal, lar, lob, lot = latleft[i], latright[i], lonbottom[i], lontop[i]
        if reproject == False:
            w = w_all[i]
        else:
            w = w_all[4]
        if image_type[i] == 'r':
            cmap = 'inferno'
        else:
            cmap = 'cubehelix'
        ax = fig.add_axes([horizontal_pad[i], vertical_pad[i], width, height], projection=w)
        norm = simple_norm(imdata_all[i], percent=99.5)
        ax.imshow(imdata_all[i], cmap=cmap, norm=norm)
        ax0 = fig.add_axes([horizontal_pad0[i], vertical_pad0[i], 0.5*width, 0.1*height])
        hide_ax(ax0)
        ax0.text(0.1, 0.3, data_name_all[i], fontsize=18)
        if image_type[i] != 'r':
            for ii in range(0, len(contour)):
                if contour[ii]:
                    imdata, footprint = reproject_interp(hdu_all[ii], header_all[4])
                    levels = np.geomspace(3*image_rms[ii], max_contour_lev[ii]*image_rms[ii], n_contour_levs[ii])
                    plot_contour(ax, imdata, levels, contour_lw[ii], contour_colors[ii], contour_alpha[ii])

        lon = ax.coords['ra']
        lat = ax.coords['dec']
        if lob:
            lon.set_axislabel('Right Ascension',fontsize=20)
            lon.set_major_formatter('hh:mm')
        if lal:
            lat.set_axislabel('Declination',minpad=-1,fontsize=20)
            lat.set_major_formatter('dd:mm')
        lon.set_ticklabel_visible(lob)
        lat.set_ticklabel_visible(lal)
        lon.set_ticks_visible(lot)
        lat.set_ticks_visible(lal)
        overlay = ax.get_coords_overlay('fk5')
        overlay.grid(axes=ax, color='white', ls='dotted',alpha=1)
        overlay['ra'].set_ticklabel_visible(lot)
        
        if lot:
            overlay['ra'].set_major_formatter('hh:mm')
        overlay['dec'].set_ticklabel_visible(lar)
        if lar:
            overlay['dec'].set_major_formatter('dd:mm')
        lon.set_ticklabel(size=20)
        lat.set_ticklabel(size=20)

    plt.show()

viewer('MIDAS_J3425652-332012', 500)