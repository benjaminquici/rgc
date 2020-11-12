#! /usr/bin/env python

"""
Create and edit the source table.

"""

import numpy as np
import matplotlib.pyplot as mpl


def create_source_table(colname, colunit, coltype, savename, overwrite):
    """
    Creates an empty (len=1) fits table containing the columns parsed through argparse.

    parameters
    ----------
    colname : list
        List containing the column names

    colunit : list
        List containing the column units

    coltype: list
        List containing the column dtypes

    savename: str
        Name of output fits table.

    overwrite: Bool
        Overwriting permissions.
    """
    root = 'E:/Thesis/Paper_II/Pipeline/'
    t = Table()
    for i in range(0,len(colname)):
        if colunit[i] == '-':
            colunit[i] = ''
        t[colname[i]] = Column(np.zeros([1]), unit=colunit[i], dtype=coltype[i])
    savepath = '{0}/Tables/{1}.fits'.format(root, savename)
    t.write(savepath, overwrite=overwrite)
    return

def extract_data_wcs(hdu_list):
    imdata_list, header_list, wcs_list = [], [], []
    for i in range(0,len(hdu_list)):
        hdu = hdu_list[0]
        header = hdu[0].header
        imdata_list.append(hdu[0].data)
        header_list.append(header)
        wcs_list.append(wcs.WCS(header, naxis=2)
    return(imdata_list, header_list, wcs_list)


def plot_images(hdu_list, return_ax=True, scalefig=1):
    """
    Initialises an empty figure with n panels.

    parameters
    ----------
    figwidth :
        width of figure

    figheight :
        height of figure

    hdu_list :
        list containing the hdus we want to plot
    """
    n = len(hdu_list)
    if n == 1:
        imdata_list, header_list, wcs_list = extract_data_wcs(hdu_list)
        fig = plt.figure(figsize=(scalefig*5, scalefig*5))
        ax0 = fig.add_axes([0.1,0.1,0.9,0.9], projection=wcs_list[0])
        ax0.imshow(imdata_list[0], cmap=cmap)

        axes = [ax0]
        latleft = [True]
        latright = [False]
        lonbottom = [True]
        lontop = [False]

    elif n == 2:
        imdata_list, header_list, wcs_list = extract_data_wcs(hdu_list)
        ax0 = fig.add_axes([0.05, 0.05, 0.45, 0.9], projection=wcs_list[0])
        ax1 = fig.add_axes([0.5, 0.05, 0.45, 0.9], projection=wcs_list[1])

        axes = [ax0, ax1]
        latleft = [True, False]
        latright = [False, False]
        lonbottom = [True, True]
        lontop = [False, False]
    # elif n == 3:
    #     # return(fig)
    # elif n == 4:
    #     # return(fig)
    # elif n == 6:
    #     # return(fig)
    # else:
    #     # return(fig)

    axes = [ax0, ax1, ax2, ax3]
    latleft = [True, False, True, False]
    latright = [False, False, False, False]
    lonbottom = [True, True, False, False]
    lontop = [False, False, False, False]

    for ax, lal, lar, lob, lot in zip(axes, latleft, latright, lonbottom, lontop):
        lon = ax.coords['ra']
        lat = ax.coords['dec']
        if lob:
            lon.set_axislabel('Right Ascension', fontsize=10)
            lon.set_major_formatter('hh:mm:ss')
        if lal:
            lat.set_axislabel('Declination', minpad=-1, fontsize=10)
            lat.set_major_formatter('dd:mm:ss')
        lon.set_ticklabel_visible(lob)
        lat.set_ticklabel_visible(lal)
        lon.set_ticks_visible(lot)
        lat.set_ticks_visible(lal)
        overlay = ax.get_coords_overlay('fk5')
        overlay.grid(axes=ax, color='white', ls='dotted', alpha=0.25)
        overlay['ra'].set_ticklabel_visible(lot)

        if lot:
            overlay['ra'].set_major_formatter('hh:mm')
        overlay['dec'].set_ticklabel_visible(lar)
        if lar:
            overlay['dec'].set_major_formatter('dd:mm')
        lon.set_ticklabel(size=10)
        lat.set_ticklabel(size=10)
    return(fig)
