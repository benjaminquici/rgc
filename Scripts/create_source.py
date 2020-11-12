from tasks import get_glass_region
from tasks import get_emu_half
from tasks import make_box

import numpy as np
import argparse

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
mpl.use('TkAgg')

firstsnr = { "marker" : "x" , "linestyle" : "None", "color" : "white" }
restsnr = { "marker" : "None" , "linestyle" : "-", "color" : "white" }
srcmark = { "marker" : "o" , "linestyle" : "None", "color" : "yellow" }
firstexc = { "marker" : "x" , "linestyle" : "None", "color" : "blue" }
restexc = { "marker" : "x" , "linestyle" : "None", "color" : "blue" }

class Coords:
    def __init__(self):
        self.x = []
        self.y = []

class PolyPick:
    def __init__(self, ax=None):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        self.cid = ax.figure.canvas.mpl_connect('button_press_event', self)
# Lines (region to include in remnant)
        self.coords = Coords()
# Points (to find and exclude point sources)
        self.points = Coords()
# Exclude (region to exclude from background fit
        self.exclude = Coords()

    def __call__(self, event):
        if event.inaxes!=self.ax.axes: return
        if event.button == 1:
            #print("Selected coords {0:5.4f}, {1:5.4f}".format(event.xdata,event.ydata))
            self.coords.x.append(event.xdata)
            self.coords.y.append(event.ydata)
        elif event.button == 2:
            #print("Selecting source at {0:5.4f}, {1:5.4f} for fitting".format(event.xdata,event.ydata))
            self.points.x.append(event.xdata)
            self.points.y.append(event.ydata)
        elif event.button == 3:
           # print("Selecting coords at {0:5.4f}, {1:5.4f} to remove from background fit".format(event.xdata,event.ydata))
            self.exclude.x.append(event.xdata)
            self.exclude.y.append(event.ydata)

# Remnant region
# First line: draw a "x" to show it's the first one
        if len(self.coords.x) == 1:
            line, = self.ax.plot(self.coords.x, self.coords.y,**firstsnr)
# Second or further line: draw the whole line
        else:
            line, = self.ax.plot(self.coords.x, self.coords.y,**restsnr)
# Sources
        line, = self.ax.plot(self.points.x, self.points.y,**srcmark)
# Exclusion zone
        if len(self.exclude.x) == 1:
            line, = self.ax.plot(self.exclude.x, self.exclude.y,**firstexc)
# Second or further line: draw the whole line
        else:
            line, = self.ax.plot(self.exclude.x, self.exclude.y,**restexc)

# This makes it plot without needing to change focus back to the terminal
        self.ax.figure.canvas.draw()

def create_IAU_name(ra,dec,prefix=None):
    c = SkyCoord(ra*u.deg,dec*u.deg,frame='fk5')
    if prefix is None:
        prefix = 'GAMA'
    name = '{2}_J{0}{1}'.format(c.ra.to_string(unit=u.hourangle,sep='',precision=0,pad=True),c.dec.to_string(sep='',precision=0,alwayssign=True,pad=True),prefix)
    return(name)

def get_skyview(survey_name,ra,dec,size):
    from astroquery.skyview import SkyView
    pstg_stamp=SkyView.get_images("{0},{1}".format(ra,dec),radius=size*u.arcsecond,survey=survey_name)[0]
    image_header=pstg_stamp[0].header
    # print(image_header)
    if survey_name == 'NVSS':
        image_header['BMAJ']=0.0125
        image_header['BMIN']=0.0125
        image_header['BPA']=0
    # pstg_stamp.writeto(path_to_save_dir+viewer_surveys()[i]+'_{0}.fits'.format(size),overwrite=True)
    return(pstg_stamp)
    # pstg_stamp.close()

def generate_cutout(survey_name, ra, dec, size, path_to_mosaic=False, input_hdu=False, target_name=False, savehdu=False, return_hdu=False):
    if path_to_mosaic:
        hdu = fits.open(path_to_mosaic)
    else:
        hdu = input_hdu
    image_data = np.squeeze(hdu[0].data)
    image_header = hdu[0].header
    w = wcs.WCS(hdu[0].header,naxis=2)
    hdu.close()
    cutout_box = (size*u.arcsecond, size*u.arcsecond)
    position = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='fk5')
    image_data = Cutout2D(image_data, position, cutout_box, wcs=w)
    hdu = fits.PrimaryHDU(data=image_data.data, header=image_data.wcs.to_header())
    newhdu = fits.HDUList([hdu])
    newhdr = newhdu[0].header
    try:
        image_header['BMAJ']
    except KeyError as e:
        print('Beam properties not found.')
    else:
        for fitskey in ['BMAJ', 'BMIN', 'BPA']:
            newhdr[fitskey] = image_header[fitskey]
        newhdr['BUNIT'] = 'JY/BEAM'
    if savehdu:
        if target_name is False:
            target_name = generate_j2000_name(ra, dec, 'MIDAS')
        savename = '{0}_{1}_{2}.fits'.format(target_name, survey_name, size)
        validate_path('{0}/{1}'.format(get_paths('G23_cutouts_first_pass'), target_name), create_path=True)
        hdu.writeto('{0}/{1}/{2}'.format(get_paths('G23_cutouts_first_pass'), target_name, savename), overwrite=True)
    if return_hdu:
        return(newhdu)

def plot_contour(ax, data, rms, num_lev, color, lw, alpha=None):
    faint_levels = rms*np.array([3,4])
    if np.nanmax(data) > 7*rms:
        levels = np.geomspace(5*rms, np.nanmax(data), num_lev)
        ax.contour(data, levels=levels, linewidths=lw, colors=color, alpha=alpha)
    ax.contour(data, levels=faint_levels, linewidths=lw-1, colors=color, alpha=alpha)
    return(ax)

def define_source(ra, dec, options):
    """
    Interactively i) match midas comps to the same source, ii) measure the EMU and/or GLASS largest-linear-size, iii) extract the coordinates of the radio core, iv) extract the coordinates of the host galaxy.

    parameters
    ----------

    """
    root = 'E:\Thesis\mosaics'
    emu_half = get_emu_half(ra)
    glass_reg = get_glass_region(ra, dec)
    surveys = ['G23_VIKING_K-band_reg{0}'.format(glass_reg),'G23_MIDAS_215MHz','G23_EMU_887MHz_{0}'.format(emu_half),'G23_GLASS_5500MHz_reg{0}'.format(glass_reg)]

    hdu_list_0 = []
    imdata_list = []
    header_list = []
    rad = 300
    for i in range(0, len(surveys)):
        path_to_mosaic = '{0}\{1}.fits'.format(root,surveys[i])
        hdu0 = generate_cutout(surveys[i], ra, dec, rad, path_to_mosaic=path_to_mosaic, return_hdu=True)
        # hdu_list_0.append(hdu0)
        if i == 0:
            imdata_list.append(hdu0[0].data)
            header_list.append(hdu0[0].header)
            hdu_list_0.append(hdu0)
        elif i > 0:
            hdu_reproject = generate_cutout(surveys[0], ra, dec, rad, path_to_mosaic='{0}\{1}.fits'.format(root,surveys[0]), return_hdu=True)
            array, footprint = reproject_interp(hdu0[0], hdu_reproject[0].header)
            imdata_list.append(array)
            header_list.append(hdu_reproject[0].header)
            new_hdu = fits.PrimaryHDU(array)
            hdu0 = fits.HDUList([new_hdu])
            hdu0[0].header = hdu_reproject[0].header
            hdu_list_0.append(hdu0)

    hdunvss = get_skyview('NVSS',ra ,dec, rad)
    imdata_nvss, footprint = reproject_interp(hdunvss[0], hdu_reproject[0].header)
    fig = plt.figure(figsize=(20, 20))
    ax0 = fig.add_axes([0.025, 0.035, 0.45, 0.45], projection=wcs.WCS((header_list[0]), naxis=2))
    ax0.imshow(imdata_nvss, cmap='inferno')  # g
    ax1 = fig.add_axes([0.025 + 0.45 + 0.025, 0.035, 0.45, 0.45], projection=wcs.WCS((header_list[0]), naxis=2))
    ax1.imshow(imdata_list[0], cmap='inferno')  # v
    norm = simple_norm(imdata_list[0], percent=99.5)
    ax1.imshow(imdata_list[0], vmin=100, vmax=2000, cmap='cubehelix', norm=norm)
    ax2 = fig.add_axes([0.025, 0.03 + 0.45 + 0.035, 0.45, 0.45], projection=wcs.WCS((header_list[0]), naxis=2))
    ax2.imshow(imdata_list[1], cmap='inferno')  # m
    ax3 = fig.add_axes([0.025 + 0.45 + 0.025, 0.035 + 0.45 + 0.025, 0.45, 0.45],
                       projection=wcs.WCS((header_list[0]), naxis=2))
    plot_contour(ax1, imdata_list[2], 25e-5, 7, 'magenta', 2, alpha=None)
    plot_contour(ax1, imdata_list[3], 25e-5, 7, 'cyan', 2, alpha=None)
    plot_contour(ax1, imdata_list[1], 1e-3, 3, 'yellow', 1, alpha=None)
    ax3.imshow(imdata_list[2], cmap='inferno')  # e
    if options.midas_cat:
        # put this into a function like overlay_coords()
        print('MIDAS catalogue supplied.')
        midas_table_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.csv'.format(options.midas_cat)
        atlas = Table.read(midas_table_path,format='csv')
        midas_ra_col, midas_dec_col = atlas[options.midas_ra], atlas[options.midas_dec]
        midas_ra, midas_dec = make_box(ra, dec,rad*(3./5.)/3600.,midas_ra_col, midas_dec_col)
        ax2.scatter(midas_ra, midas_dec, transform=ax2.get_transform('fk5'), s=100, edgecolor='red', facecolor='None',
                    linestyle='-', linewidth=2, marker='o')
    else:
        print('No MIDAS catalogue given.')

    if options.out_cat:
        print('Output source catalogue supplied.')
        output_table_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.fits'.format(options.out_cat)
        atlas = Table.read(output_table_path)
        output_ra_col, output_dec_col = atlas[options.out_ra], atlas[options.out_dec]
        out_ra, out_dec = make_box(ra, dec, rad*(3./5.) / 3600., output_ra_col, output_dec_col)
        if len(out_ra)>0:
            ax2.scatter(out_ra, out_dec, transform=ax2.get_transform('fk5'), s=100, edgecolor='blue', facecolor='blue',
                    linestyle='-', linewidth=2, marker='x')
    else:
        print('No output source catalogue given.')

    if options.nvss_cat:
        print('NVSS catalogue supplied.')
        nvss_table_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.fits'.format(options.nvss_cat)
        atlas = Table.read(nvss_table_path)
        nvss_ra_col, nvss_dec_col = atlas[options.out_ra], atlas[options.out_dec]
        nvss_ra, nvss_dec = make_box(ra, dec, rad*(3./5.) / 3600., nvss_ra_col, nvss_dec_col)
        if len(nvss_ra)>0:
            ax0.scatter(nvss_ra, nvss_dec, transform=ax0.get_transform('fk5'), s=100, edgecolor='red', facecolor='None',
                    linestyle='-', linewidth=2, marker='o')
    else:
        print('No output source catalogue given.')

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

    polypick = PolyPick(ax2)
    polypick1 = PolyPick(ax3)
    polypick2 = PolyPick(ax0)

    plt.show()
    LAS_points = Coords()
    radio_center_coords = Coords()
    midas_comp_coords = Coords()
    nvss_comp_coords = Coords()
    if len(polypick.coords.x):
        if len(polypick.coords.x):
            hdu_1_list = []
            radio_center_coords.x, radio_center_coords.y = (wcs.WCS((header_list[0]), naxis=2)).wcs_pix2world((polypick1.points.x)[len(polypick1.points.x)-1], (polypick1.points.y)[len(polypick1.points.x)-1], 0)
            LAS_points.x, LAS_points.y = (wcs.WCS((header_list[0]), naxis=2)).wcs_pix2world(polypick.coords.x, polypick.coords.y, 0)
            c1 = SkyCoord((LAS_points.x)[0] * u.deg, (LAS_points.y)[0] * u.deg, frame='fk5')
            c2 = SkyCoord((LAS_points.x)[1] * u.deg, (LAS_points.y)[1] * u.deg, frame='fk5')
            sep = c1.separation(c2)
        if len(polypick.points.x):
            midas_comp_coords.x, midas_comp_coords.y = (wcs.WCS((header_list[0]), naxis=2)).wcs_pix2world((polypick.points.x), (polypick.points.y), 0)
        # else:
        #     pass

        if len(polypick2.points.x):
            nvss_comp_coords.x, nvss_comp_coords.y = (wcs.WCS((header_list[0]), naxis=2)).wcs_pix2world((polypick2.points.x), (polypick2.points.y), 0)
        # else:
        #     pass

        for i in range(0,len(hdu_list_0)):
            hdu1 = generate_cutout('DUMMY', radio_center_coords.x, radio_center_coords.y, sep.arcsecond, input_hdu=hdu_list_0[i], return_hdu=True)
            hdu_1_list.append(hdu1)
        fig = plt.figure(figsize=(20, 10))
        hdu_emu = hdu_1_list[2]
        hdu_glass = hdu_1_list[3]
        header_emu = hdu_emu[0].header
        w_emu = wcs.WCS(header_emu, naxis=2)
        header_glass = hdu_glass[0].header
        w_glass = wcs.WCS(header_glass, naxis=2)
        ax0 = fig.add_axes([0.05, 0.05, 0.45, 0.9], projection=w_emu)
        ax0.imshow(hdu_emu[0].data, cmap='inferno')

        ax1 = fig.add_axes([0.05+0.45, 0.05, 0.45, 0.9], projection=w_glass)
        ax1.imshow(hdu_glass[0].data, cmap='inferno')

        polypick_emu = PolyPick(ax0)
        polypick_glass = PolyPick(ax1)

        plt.show()

        emu_las = Coords()
        glass_las = Coords()
        emu_core = Coords()
        glass_core = Coords()

        if len(polypick_emu.coords.x):
            emu_las.x, emu_las.y = w_emu.wcs_pix2world(polypick_emu.coords.x, polypick_emu.coords.y, 0)
            c1 = SkyCoord((emu_las.x)[0] * u.deg, (emu_las.y)[0] * u.deg, frame='fk5')
            c2 = SkyCoord((emu_las.x)[1] * u.deg, (emu_las.y)[1] * u.deg, frame='fk5')
            sep = c1.separation(c2)
            emu_sep = float('%.2g' % (sep.arcsecond))
            print("{0} arcsecond LLS measured from EMU".format(emu_sep))
        else:
            emu_sep = np.nan
            print("EMU LLS not measured")

        if len(polypick_glass.coords.x):
            glass_las.x, glass_las.y = w_glass.wcs_pix2world(polypick_glass.coords.x, polypick_glass.coords.y, 0)
            c1 = SkyCoord((glass_las.x)[0] * u.deg, (glass_las.y)[0] * u.deg, frame='fk5')
            c2 = SkyCoord((glass_las.x)[1] * u.deg, (glass_las.y)[1] * u.deg, frame='fk5')
            sep = c1.separation(c2)
            glass_sep = float('%.2g' % (sep.arcsecond))
            print("{0} arcsecond LLS measured from GLASS".format(glass_sep))
        else:
            glass_sep = np.nan
            print("GLASS LLS not measured")

        if len(polypick_emu.points.x):
            emu_core.x, emu_core.y = w_emu.wcs_pix2world((polypick_emu.points.x)[len(polypick_emu.points.x) - 1], (polypick_emu.points.y)[len(polypick_emu.points.x) - 1], 0)
            ra_emu_core, dec_emu_core= emu_core.x, emu_core.y
            print("EMU core measured at RA = {0}, Dec = {1}".format(ra_emu_core, dec_emu_core))
        else:
            ra_emu_core, dec_emu_core = np.nan, np.nan
            print("EMU core not measured")

        if len(polypick_glass.points.x):
            glass_core.x, glass_core.y = w_glass.wcs_pix2world((polypick_glass.points.x)[len(polypick_glass.points.x) - 1],(polypick_glass.points.y)[len(polypick_glass.points.x) - 1], 0)
            ra_glass_core, dec_glass_core = glass_core.x, glass_core.y
            print("GLASS core measured at RA = {0}, Dec = {1}".format(ra_glass_core, dec_glass_core))
        else:
            ra_glass_core, dec_glass_core = np.nan, np.nan
            print("GLASS core not measured")

            # now make a zoomed in image around the "host galaxy"
        hdu_2_list = []
        if emu_sep != np.nan:
            host_cutout_rad = emu_sep
        else:
            host_cutout_rad = 60
        for i in range(0,len(hdu_1_list)):
            hdu2 = generate_cutout('DUMMY', radio_center_coords.x, radio_center_coords.y,host_cutout_rad , input_hdu=hdu_1_list[i], return_hdu=True)
            hdu_2_list.append(hdu2)
        fig = plt.figure(figsize=(10, 10))
        hdu_vik, hdu_emu, hdu_glass = hdu_2_list[0], hdu_2_list[2], hdu_2_list[3]
        imdata_vik, header_vik, imdata_emu, imdata_glass = hdu_vik[0].data, hdu_vik[0].header, hdu_emu[0].data, hdu_glass[0].data
        w = wcs.WCS(header_vik, naxis=2)
        norm = simple_norm(imdata_vik, percent=99.5)
        ax0 = fig.add_axes([0.075, 0.075, 0.9, 0.9], projection=w)
        ax0.imshow(imdata_vik, vmin=10, vmax=2000, cmap='cubehelix', norm=norm)
        plot_contour(ax0,imdata_emu,35e-5,7,'magenta',2)
        plot_contour(ax0, imdata_glass, 25e-5, 7, 'cyan', 2)
        if ra_emu_core != np.nan:
            ax0.scatter(ra_emu_core, dec_emu_core, transform=ax0.get_transform('fk5'), s=70, edgecolor='magenta',facecolor="None", linestyle='-', linewidth=1, marker='o')
            # else:
            #     pass
        if ra_glass_core != np.nan:
            ax0.scatter(ra_glass_core, dec_glass_core, transform=ax0.get_transform('fk5'), s=70, edgecolor='cyan',facecolor="None", linestyle='-', linewidth=1, marker='o')
            # else:
            #     pass
        if ra_emu_core != np.nan and ra_glass_core != np.nan:
            ax0.scatter(radio_center_coords.x, radio_center_coords.y, transform=ax0.get_transform('fk5'), s=70, edgecolor='red',facecolor="None", linestyle='-', linewidth=1, marker='o')

        if options.gama_cat:
            print('GAMA23 catalogue supplied.')
            gama_table_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.fits'.format(options.gama_cat)
            atlas = Table.read(gama_table_path)
            gama_ra_col, gama_dec_col = atlas[options.gama_ra], atlas[options.gama_dec]
            gama_ra, gama_dec = make_box(radio_center_coords.x, radio_center_coords.y, 0.5*host_cutout_rad/ 3600., gama_ra_col, gama_dec_col)
            if len(gama_ra) > 0:
                ax0.scatter(gama_ra, gama_dec, transform=ax0.get_transform('fk5'), s=70, edgecolor='blue',facecolor='yellow',linestyle='-', linewidth=1, marker='x')
                # else:
                #     pass
        else:
            print('No GAMA23 catalogue given.')

        polypick_host = PolyPick(ax0)
        plt.show()

        hostcoords_gama = Coords()
        hostcoords_nogama = Coords()
        ra_host, dec_host, host_in_gama_cat = [], [], []
        if len(polypick_host.points.x):
            hostcoords_gama.x, hostcoords_gama.y = w.wcs_pix2world(polypick_host.points.x, polypick_host.points.y, 0)
            for i in range(0,len(hostcoords_gama.x)):
                print('Found GAMA host at RA = {0}, Dec = {1}'.format((hostcoords_gama.x)[i], (hostcoords_gama.y)[i]))
                ra_host.append((hostcoords_gama.x)[i])
                dec_host.append((hostcoords_gama.y)[i])
                host_in_gama_cat.append('y')
        elif len(polypick_host.exclude.x):
            hostcoords_nogama.x, hostcoords_nogama.y = w.wcs_pix2world(polypick_host.exclude.x, polypick_host.exclude.y,0)
            for i in range(0, len(hostcoords_nogama.x)):
                print('Found non-GAMA host at RA = {0}, Dec = {1}'.format((hostcoords_nogama.x)[i], (hostcoords_nogama.y)[i]))
                ra_host.append((hostcoords_nogama.x)[i])
                dec_host.append((hostcoords_nogama.y)[i])
                host_in_gama_cat.append('n')
        else:
            print('No host found.')
            ra_host.append(np.nan)
            dec_host.append(np.nan)
            host_in_gama_cat.append(0)

        print(ra_host, dec_host, host_in_gama_cat)

        # now add properties to the table.
        new_name = create_IAU_name(radio_center_coords.x, radio_center_coords.y, prefix='MIDAS')
        if options.out_cat:
            print('Adding new source to output source table.')
            atlas = Table.read(output_table_path)
            value_added_colnames = [options.out_name,options.out_ra, options.out_dec, options.out_emu_lls, options.out_glass_lls, options.out_emu_core_ra, options.out_emu_core_dec, options.out_glass_core_ra, options.out_glass_core_dec]
            valued_added_cols = [new_name, radio_center_coords.x, radio_center_coords.y, emu_sep, glass_sep, ra_emu_core, dec_emu_core, ra_glass_core, dec_glass_core]
            for i in range(0,len(ra_host)):
                atlas.add_row()
                for j in range(0, len(value_added_colnames)):
                    atlas[value_added_colnames[j]][len(atlas)-1] = valued_added_cols[j]
                atlas[options.out_host_ra][len(atlas)-1] = ra_host[i]
                atlas[options.out_host_dec][len(atlas)-1] = dec_host[i]
                atlas[options.host_in_gama_cat][len(atlas)-1] = host_in_gama_cat[i]
                if options.out_midas_cat:
                    atlas[options.out_nmidas][len(atlas)-1] = len(midas_comp_coords.x)
            # print(atlas[len(atlas)-1])
            atlas.write(output_table_path, overwrite=True)
            # else:
            #     pass
        if options.out_midas_cat:
            print('Adding individual MIDAS component coordinates to verbose table.')
            out_midas_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.fits'.format(options.out_midas_cat)
            atlas = Table.read(out_midas_path)
            if len(midas_comp_coords.x) < 1:
                atlas.add_row()
                atlas[options.out_midas_cat_name][len(atlas) - 1] = new_name
                atlas[options.out_midas_cat_ra][len(atlas) - 1] = -1
                atlas[options.out_midas_cat_dec][len(atlas) - 1] = -1
            else:
                for ra, dec in zip(midas_comp_coords.x, midas_comp_coords.y):
                    atlas.add_row()
                    atlas[options.out_midas_cat_name][len(atlas) - 1] = new_name
                    atlas[options.out_midas_cat_ra][len(atlas) - 1] = ra
                    atlas[options.out_midas_cat_dec][len(atlas) - 1] = dec


            atlas.write(out_midas_path, overwrite=True)
            # else:
            #     pass

        if options.out_nvss_cat:
            print('Adding individual NVSS component coordinates to verbose table.')
            out_nvss_path = 'E:\Thesis\Paper_II\Pipeline\Tables\{0}.fits'.format(options.out_nvss_cat)
            atlas = Table.read(out_nvss_path)
            if len(nvss_comp_coords.x) < 1:
                atlas.add_row()
                atlas[options.out_nvss_cat_name][len(atlas) - 1] = new_name
                atlas[options.out_nvss_cat_ra][len(atlas) - 1] = -1
                atlas[options.out_nvss_cat_dec][len(atlas) - 1] = -1
            else:
                for ra, dec in zip(nvss_comp_coords.x, nvss_comp_coords.y):
                    atlas.add_row()
                    atlas[options.out_nvss_cat_name][len(atlas)-1] = new_name
                    atlas[options.out_nvss_cat_ra][len(atlas)-1] = ra
                    atlas[options.out_nvss_cat_dec][len(atlas) - 1] = dec

            atlas.write(out_nvss_path, overwrite=True)
            # else:
            #     pass
    else:
        print('Skipping source.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prefix_chars='-')
    
    # Target coordinates
    group1 = parser.add_argument_group("RA target")
    group1.add_argument("--ra_target", dest='ra_target', type=float, default=None, help='')

    group2 = parser.add_argument_group("Dec target")
    group2.add_argument("--dec_target", dest='dec_target', type=float, default=None, help='')

    group01 = parser.add_argument_group("Coords list")
    group01.add_argument("--coords_list", dest='coords_list', type=str, default=None, help='')
    #
    # group02 = parser.add_argument_group("Dec target list")
    # group02.add_argument("--dec_target_list", dest='dec_target_list',  nargs='+', type=float, default=None, help='')

    # midas table
    group3 = parser.add_argument_group("MIDAS catalogue")
    group3.add_argument("--midas_cat", dest='midas_cat', type=str, default=None, help='')

    group4 = parser.add_argument_group("MIDAS RA col")
    group4.add_argument("--midas_ra", dest='midas_ra', type=str, default='ra', help='')

    group5 = parser.add_argument_group("MIDAS DEC col")
    group5.add_argument("--midas_dec", dest='midas_dec', type=str, default='dec', help='')

    # GAMA23 galaxy table
    group6 = parser.add_argument_group("GAMA catalogue")
    group6.add_argument("--gama_cat", dest='gama_cat', type=str, default=None, help='')

    group7 = parser.add_argument_group("MIDAS RA col")
    group7.add_argument("--gama_ra", dest='gama_ra', type=str, default='RAcen', help='')

    group8 = parser.add_argument_group("MIDAS DEC col")
    group8.add_argument("--gama_dec", dest='gama_dec', type=str, default='Deccen', help='')

    # output source table
    group9 = parser.add_argument_group("Output source catalogue")
    group9.add_argument("--out_cat", dest='out_cat', type=str, default=None, help='')

    group10 = parser.add_argument_group("Output RA col")
    group10.add_argument("--out_ra", dest='out_ra', type=str, default='ra', help='')

    group11 = parser.add_argument_group("Output Dec col")
    group11.add_argument("--out_dec", dest='out_dec', type=str, default='dec', help='')

    group12 = parser.add_argument_group("Output EMU LLS col")
    group12.add_argument("--out_emu_lls", dest='out_emu_lls', type=str, default='emu_LLS', help='')

    group13 = parser.add_argument_group("Output GLASS LLS col")
    group13.add_argument("--out_glass_lls", dest='out_glass_lls', type=str, default='glass_LLS', help='')

    group14 = parser.add_argument_group("Output EMU core RA col")
    group14.add_argument("--out_emu_core_ra", dest='out_emu_core_ra', type=str, default='ra_core_887', help='')

    group15 = parser.add_argument_group("Output EMU core Dec col")
    group15.add_argument("--out_emu_core_dec", dest='out_emu_core_dec', type=str, default='dec_core_887', help='')

    group16 = parser.add_argument_group("Output GLASS core RA col")
    group16.add_argument("--out_glass_core_ra", dest='out_glass_core_ra', type=str, default='ra_core_5500', help='')

    group17 = parser.add_argument_group("Output GLASS core Dec col")
    group17.add_argument("--out_glass_core_dec", dest='out_glass_core_dec', type=str, default='dec_core_5500', help='')

    group18 = parser.add_argument_group("Output Host RA col")
    group18.add_argument("--out_host_ra", dest='out_host_ra', type=str, default='ra_host', help='')

    group19 = parser.add_argument_group("Output Host Dec col")
    group19.add_argument("--out_host_dec", dest='out_host_dec', type=str, default='dec_host', help='')

    group20 = parser.add_argument_group("Output name col")
    group20.add_argument("--out_name", dest='out_name', type=str, default='Name', help='')

    group21 = parser.add_argument_group("Output gamaid col")
    group21.add_argument("--host_in_gama_cat", dest='host_in_gama_cat', type=str, default='host_in_gama_cat', help='')

    group26 = parser.add_argument_group("Output Nmidas col")
    group26.add_argument("--out_nmidas", dest='out_nmidas', type=str, default='Nmidas_cmps', help='')

    # output midas table
    group22 = parser.add_argument_group("Verbose MIDAS cmps cat")
    group22.add_argument("--out_midas_cat", dest='out_midas_cat', type=str, default='', help='')

    group23 = parser.add_argument_group("MIDAS cmps RA")
    group23.add_argument("--out_midas_cat_ra", dest='out_midas_cat_ra', type=str, default='ra', help='')

    group24 = parser.add_argument_group("MIDAS cmps dec")
    group24.add_argument("--out_midas_cat_dec", dest='out_midas_cat_dec', type=str, default='dec', help='')

    group25 = parser.add_argument_group("Verbose MIDAS cmps cat")
    group25.add_argument("--out_midas_cat_name", dest='out_midas_cat_name', type=str, default='Name', help='')

    # nvss table
    group30 = parser.add_argument_group("NVSS catalogue")
    group30.add_argument("--nvss_cat", dest='nvss_cat', type=str, default=None, help='')

    group31 = parser.add_argument_group("NVSS RA col")
    group31.add_argument("--nvss_ra", dest='nvss_ra', type=str, default='ra', help='')

    group32 = parser.add_argument_group("NVSS DEC col")
    group32.add_argument("--nvss_dec", dest='nvss_dec', type=str, default='dec', help='')

    # output nvss table
    group33 = parser.add_argument_group("Verbose NVSS cmps cat")
    group33.add_argument("--out_nvss_cat", dest='out_nvss_cat', type=str, default='', help='')

    group34 = parser.add_argument_group("MIDAS cmps RA")
    group34.add_argument("--out_nvss_cat_ra", dest='out_nvss_cat_ra', type=str, default='ra', help='')

    group35 = parser.add_argument_group("MIDAS cmps dec")
    group35.add_argument("--out_nvss_cat_dec", dest='out_nvss_cat_dec', type=str, default='dec', help='')

    group36 = parser.add_argument_group("Verbose MIDAS cmps cat")
    group36.add_argument("--out_nvss_cat_name", dest='out_nvss_cat_name', type=str, default='Name', help='')


    # add+
    # run bane option

    options = parser.parse_args()

    if options.ra_target:
        define_source(options.ra_target, options.dec_target, options)
    elif options.coords_list:
        path_to_coords_file = 'E:\Thesis\Paper_II\Pipeline\{0}'.format(options.coords_list)
        with open(path_to_coords_file) as f:
            lines = [line.rstrip() for line in f]
            for i in range(0,len(lines)):
                ra, dec = lines[i].split(',')
                ra, dec = float(ra), float(dec)
                # print(ra, dec)
                define_source(ra, dec, options)