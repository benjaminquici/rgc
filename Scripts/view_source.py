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

atlas = Table.read('/mnt/e/Thesis/Paper_II/Pipeline/Tables/source_cat_first.fits')
name_col, ra_col, dec_col, las_col = atlas['Name'], atlas['ra'], atlas['dec'], atlas['emu_LLS']

levels=np.array([3,4,5,10,20,100])
mean_rms_radio=[0, 0.047, 0.035, 1]
contour_color=['r', 'blue', 'magenta', 'yellow']
lw=[1, 1.5 ,2 , 1]

for i in range(228, 229):#,len(name_col)):
    name, ra, dec, las = name_col[i], ra_col[i], dec_col[i], las_col[i]
    print(name)
    glass_reg = get_glass_region(ra,dec)
    emu_half = get_emu_half(ra)
    mosaic_names = ['G23_VIKING_K-band_reg{0}'.format(glass_reg), 'G23_EMU_887MHz_{0}'.format(emu_half), 'G23_GLASS_5500MHz_reg{0}'.format(glass_reg)]#, 'G23_MIDAS_215MHz']
    hdu_list_0 = []
    for ii in range(0,len(mosaic_names)):
        path_to_mosaic = '/mnt/e/Thesis/mosaics/{0}.fits'.format(mosaic_names[ii])
        hdu0 = generate_cutout('prefix', ra, dec, 1.8*las, path_to_mosaic=path_to_mosaic, return_hdu=True)
        hdu_list_0.append(hdu0)
    vik_imdata, vik_w = hdu_list_0[0][0].data, wcs.WCS(hdu_list_0[0][0].header, naxis=2)
    fig = plt.figure(figsize=(10,10))
    ax0 = fig.add_axes([0.09, 0.09, 0.9, 0.9], projection=vik_w)
    norm = simple_norm(vik_imdata, percent=99.5)
    ax0.imshow(vik_imdata, cmap='cubehelix', norm=norm)
    for ii in range(1,2):#len(mosaic_names)):
        # array, w = hdu_list_0[ii][0].data, wcs.WCS(hdu_list_0[ii][0].header, naxis=2)
        # ax0 = fig.add_axes([0.09, 0.09, 0.9, 0.9], projection=w)
        array, footprint = reproject_interp(hdu_list_0[ii][0], hdu_list_0[0][0].header)
        # ax0.contour(array,levels=levels*mean_rms_radio[ii]*1e-3, colors=contour_color[ii], lw=lw[ii])
        # img_data = get_hdu(images[0],reproject=True)[0]
        alphas = np.ones(array.shape)
        alphas[np.where(array<(5*mean_rms_radio[ii]*1e-3))]=0 # The value used here is the threshold below which pixels are automatically set to 'invisible'
        cmap = matplotlib.cm.inferno
        colors = Normalize(0.0001, 0.004, clip=True)(array) # I think these values work well for the jet and viridis colormap
        colors = cmap(colors)
        colors[...,-1]=alphas
        ax0.imshow(colors,cmap=cmap,alpha=0.8)
        ax0.contour(array, levels=47e-6*np.array([5, 10, 30]),lw=0.5, colors='white', alpha=0.35)
        hdu = fits.PrimaryHDU(array)
        hdu = fits.HDUList([hdu])
        hdu[0].header = hdu_list_0[ii][0].header
        hdu[0].header['BMAJ'] = hdu_list_0[ii][0].header['BMAJ']
        hdu[0].header['BMIN'] = hdu_list_0[ii][0].header['BMIN']
        hdu[0].header['BPA'] = hdu_list_0[ii][0].header['BPA']
        hdu[0].header['BUNIT'] = hdu_list_0[ii][0].header['BUNIT']
        hdu[0].header['CDELT2'] = hdu_list_0[ii][0].header['CDELT2']
        hdu[0].header['CDELT1'] = hdu_list_0[ii][0].header['CDELT1']
        hdu.writeto('/mnt/e/Thesis/Paper_II/Pipeline/Scripts/polygon-flux/image.fits', overwrite=True)

    plt.show()
    # print(ing)
    # plt.savefig('/mnt/e/Thesis/Paper_II/Pipeline/Sources/{0}.png'.format(name), dpi=400)