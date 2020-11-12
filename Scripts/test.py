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

import matplotlib as mpl

import os

import warnings

warnings.filterwarnings('ignore')

firstLAS = {"marker": "x", "linestyle": "None", "color": "blue"}
restLAS = {"marker": "x", "linestyle": ":", "color": "blue"}
eradiocore = {"marker": "d", "linestyle": "None", "facecolor": "None", "edgecolor":"red", "s":100}
uradiocore = {"marker": "o", "linestyle": "None", "facecolor": "None", "edgecolor":"red", "s":100}
eradiocore_c = {"marker": ".", "linestyle": "None", "facecolor": "None", "edgecolor":"red", "s":10}
uradiocore_c = {"marker": ".", "linestyle": "None", "facecolor": "None", "edgecolor":"red", "s":10}
radiocenter = {"marker": "o", "linestyle": "None", "color": "red"}
xmatchr = {"marker": "x", "linestyle": "None", "color": "yellow"}
xmatchh = {"marker": "x", "linestyle": "None", "color": "yellow"} # host, with GAMA
xmatchhng = {"marker": "x", "linestyle": "None", "color": "blue"} #host, no GAMA
unrelatedcmps = {"marker": "*", "linestyle": "None", "color": "red"}


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
        self.cid = ax.figure.canvas.mpl_connect('key_press_event', self)
        self.las = Coords()        # Set of points to evaluate the largest angular size
        self.ecore = Coords()      # Extended radio core coordinates
        self.ucore = Coords()      # Unresolved radio core coordinates
        self.geomcent = Coords()   # Geometric radio center
        self.xmatch_r = Coords()   # Radio components to cross match
        self.xmatch_h = Coords()   # Host galaxies in the GAMA23 photometry catalogue to cross match
        self.xmatch_hng = Coords() # Host galaxies not in the GAMA23 photometry catalogue to cross match
        self.blended = Coords()    # Unrelated radio sources blended at lower resolution

    def __call__(self, event):
        if event.inaxes != self.ax.axes: return
        if event.key == 'l':
            self.las.x.append(event.xdata)
            self.las.y.append(event.ydata)
        elif event.key == 'u':
            self.ucore.x.append(event.xdata)
            self.ucore.y.append(event.ydata)
        elif event.key == 'e':
            self.ecore.x.append(event.xdata)
            self.ecore.y.append(event.ydata)
        elif event.key == 'g':
            self.geomcent.x.append(event.xdata)
            self.geomcent.y.append(event.ydata)
        elif event.key == 'x':
            self.xmatch_r.x.append(event.xdata)
            self.xmatch_r.y.append(event.ydata)
        elif event.key == 'h':
            self.xmatch_h.x.append(event.xdata)
            self.xmatch_h.y.append(event.ydata)
        elif event.key == 'n':
            self.xmatch_hng.x.append(event.xdata)
            self.xmatch_hng.y.append(event.ydata)
        elif event.key == 'b':
            self.blended.x.append(event.xdata)
            self.blended.y.append(event.ydata)

        # Some simple logic to draw a line on the canvas (relevant only for the LAS)
        if len(self.las.x) == 1:
            line, = self.ax.plot(self.las.x, self.las.y, **firstLAS)
        # Second or further line: draw the whole line
        else:
            line, = self.ax.plot(self.las.x, self.las.y, **restLAS)

        # This draws points onto the canvas
        self.ax.scatter(self.ecore.x, self.ecore.y, **eradiocore)
        self.ax.scatter(self.ucore.x, self.ucore.y, **uradiocore)
        self.ax.scatter(self.ecore.x, self.ecore.y, **eradiocore_c)
        self.ax.scatter(self.ucore.x, self.ucore.y, **uradiocore_c)
        self.ax.scatter(self.geomcent.x, self.geomcent.y, **radiocenter)
        self.ax.scatter(self.xmatch_r.x, self.xmatch_r.y, **xmatchr)
        self.ax.scatter(self.xmatch_h.x, self.xmatch_h.y, **xmatchh)
        self.ax.scatter(self.xmatch_hng.x, self.xmatch_hng.y, **xmatchhng)
        self.ax.scatter(self.blended.x, self.blended.y, **unrelatedcmps)

        # This makes it plot without needing to change focus back to the terminal
        self.ax.figure.canvas.draw()

def create_IAU_name(ra, dec, prefix=None):
    c = SkyCoord(ra * u.deg, dec * u.deg, frame='fk5')
    if prefix is None:
        prefix = 'GAMA'
    name = '{2}_J{0}{1}'.format(c.ra.to_string(unit=u.hourangle, sep='', precision=0, pad=True),
                                c.dec.to_string(sep='', precision=0, alwayssign=True, pad=True), prefix)
    return (name)

def compute_angular_separation(ra1, dec1, ra2, dec2):
    c1 = SkyCoord(ra1*u.deg, dec1*u.deg, frame='fk5')
    c2 = SkyCoord(ra2*u.deg, dec2*u.deg, frame='fk5')
    sep = c1.separation(c2)
    return(sep)

def get_skyview(survey_name, ra, dec, size):
    from astroquery.skyview import SkyView
    pstg_stamp = SkyView.get_images("{0},{1}".format(ra, dec), radius=size * u.arcsecond, survey=survey_name)[0]
    image_header = pstg_stamp[0].header
    # print(image_header)
    if survey_name == 'NVSS':
        image_header['BMAJ'] = 0.0125
        image_header['BMIN'] = 0.0125
        image_header['BPA'] = 0
    # pstg_stamp.writeto(path_to_save_dir+viewer_surveys()[i]+'_{0}.fits'.format(size),overwrite=True)
    return (pstg_stamp)
    # pstg_stamp.close()

def generate_cutout(survey_name, ra, dec, size, path_to_mosaic=None, input_hdu=None, target_name=False, savehdu=False,return_hdu=False):
    if path_to_mosaic:
        hdu = fits.open(path_to_mosaic)
    else:
        hdu = input_hdu
    image_data = np.squeeze(hdu[0].data)
    image_header = hdu[0].header
    w = wcs.WCS(hdu[0].header, naxis=2)
    # hdu.close()
    cutout_box = (size * u.arcsecond, size * u.arcsecond)
    position = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='fk5')
    image_data = Cutout2D(image_data, position, cutout_box, wcs=w)
    hdu = fits.PrimaryHDU(data=image_data.data, header=image_data.wcs.to_header())
    hdu = fits.HDUList([hdu])
    newhdr = hdu[0].header
    try:
        image_header['BMAJ']
    except KeyError as e:
        print('Beam properties not found, must be VIKING image.')
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
        return(hdu)

def plot_contour(ax, data, rms, num_lev, color, lw, alpha=None):
    faint_levels = rms * np.array([3, 4])
    if np.nanmax(data) > 7 * rms:
        levels = np.geomspace(5 * rms, np.nanmax(data), num_lev)
        ax.contour(data, levels=levels, linewidths=lw, colors=color, alpha=alpha)
    ax.contour(data, levels=faint_levels, linewidths=lw - 1, colors=color, alpha=alpha)
    return (ax)

def extract_radio_properties(hdu, survey_name, rms=None,savedir=False):

    fig = plt.figure(figsize=(18,18))
    imdata, w = hdu[0].data, wcs.WCS(hdu[0].header, naxis=2)
    norm = simple_norm(imdata, percent=99.5)
    ax0 = fig.add_axes([0.1, 0.1, 0.85, 0.85], projection=w)
    ax0.imshow(imdata,cmap='inferno',norm=norm)
    if rms:
        ax0.contour(imdata, levels=[3,5]*rms, colors='white', lw=1)
    polypick = PolyPick(ax0)
    if savedir==True:
        for ext in ['.pdf', '.png']:
            fig.savefig('{0}/{1}_LLS_plot{2}'.format(savedir, survey_name, ext))
    plt.show()
    if len(polypick.las.x):
        LAS_points = Coords()
        LAS_points.x, LAS_points.y = w.wcs_pix2world(polypick.las.x,polypick.las.y, 0)
        las = (compute_angular_separation(LAS_points.x[0], LAS_points.y[0], LAS_points.x[1], LAS_points.y[1])).arcsecond
        las = float('%.4g' % las)
        print('Measured a largest angular size of {0}'' from {1}'.format(las, survey_name))
    else:
        print('No angular size measured.')
        las = np.nan
    if len(polypick.ucore.x):
        core_coords = Coords()
        core_coords.x, core_coords.y = w.wcs_pix2world((polypick.ucore.x)[len(polypick.ucore.x) - 1],
                                                                (polypick.ucore.y)[len(polypick.ucore.x) - 1], 0)
        ra_core, dec_core = core_coords.x, core_coords.y
        core_type = 'unresolved'
        print('Unresolved {0} radio core measured at RA = {1}, Dec = {2}'.format(las, ra_core, dec_core))
    elif len(polypick.ecore.x):
        core_coords = Coords()
        core_coords.x, core_coords.y = w.wcs_pix2world((polypick.ecore.x)[len(polypick.ecore.x) - 1],
                                                       (polypick.ecore.y)[len(polypick.ecore.x) - 1], 0)
        ra_core, dec_core = core_coords.x, core_coords.y
        core_type = 'resolved'
        print('Resolved {0} radio core measured at RA = {1}, Dec = {2}'.format(las, ra_core, dec_core))
    else:
        print('No radio core measured.')
        ra_core, dec_core, core_type = np.nan, np.nan, 0
    return(las, ra_core, dec_core, core_type)

def make_radio_overlay(hdu_list, mosaic_type, survey_prefix, mean_rms, levels, contour_colour, savename=None, reprojected=False, galaxy_cat=None,interactive=True):
    fig = plt.figure(figsize=(16, 16))
    # unbind hotkeys
    if 'l' in plt.rcParams['keymap.yscale']:
        plt.rcParams['keymap.yscale'].remove('l')
        plt.rcParams['keymap.grid'].remove('g')
    # apply some logic to extract the viking hdu from list of hdus
    viking_hdu = (hdu_list[np.where(mosaic_type=='o')][0])
    # now plot the viking image
    viking_imdata, viking_w = viking_hdu[0].data, wcs.WCS(viking_hdu[0].header, naxis=2)
    ax0 = fig.add_axes([0.09, 0.09, 0.9, 0.9], projection=viking_w)
    ax1 = fig.add_axes([0.82, 0.88, 0.1, 0.1])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.patch.set_facecolor('white')
    ax1.patch.set_alpha(0.5)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    norm = simple_norm(viking_imdata, percent=99.5)
    ax0.imshow(viking_imdata, vmax=2000, cmap='gray_r', norm=norm)                                                             # make cmap an argparse option
    radio_hdu = (hdu_list[np.where(mosaic_type == 'r')])
    mean_rms_radio = (mean_rms[np.where(mosaic_type == 'r')])
    contour_color = (contour_colour[np.where(mosaic_type == 'r')])
    survey_prefix_radio = (survey_prefix[np.where(mosaic_type == 'r')])
    for i in range(0, len(radio_hdu)):
        imdata, w = radio_hdu[i][0].data, wcs.WCS(radio_hdu[i][0].header, naxis=2)
        if reprojected:
            ax0.contour(imdata, levels=levels* mean_rms_radio[i] * 1e-3, colors=contour_color[i], lw=0.5)
        else:
            ax0.contour(imdata, transform=ax0.get_transform(w),
                        levels=levels * mean_rms_radio[i] * 1e-3,
                        colors=contour_color[i], lw=0.5)
        ax1.text(0.25, 0.75 - i * 0.3, survey_prefix_radio[i], fontsize=25)
        ax1.scatter(0.1, 0.8 - i * 0.3, marker='o', s=100, facecolor=contour_color[i], edgecolor=contour_color[i])

    # Supply a galaxy catalogue to overlay the positions of host galaxies
    if galaxy_cat:
        #obtain the central position of the image from the fits header
        print('NAXIS1, 2 = {0}, {1}'.format(viking_hdu[0].header['NAXIS1'], viking_hdu[0].header['NAXIS2']))
        ra, dec = viking_w.wcs_pix2world(viking_hdu[0].header['NAXIS1']*0.5,viking_hdu[0].header['NAXIS2']*0.5, 0)

    # save the figure
    if savename:
        plt.savefig('{0}/{1}.png'.format(options.temp_dir, savename),dpi=200)

    # if interactive, then enable interactive plotting via Polypick()
    if interactive:
        return(PolyPick(ax0))

def add_psf_header_keys(newhdu, newhdr, hdu):
    newhdu[0].header = newhdr
    newhdu[0].header['BMAJ'] = hdu[0].header['BMAJ']
    newhdu[0].header['BMIN'] = hdu[0].header['BMIN']
    newhdu[0].header['BPA'] = hdu[0].header['BPA']
    newhdu[0].header['BUNIT'] = hdu[0].header['BUNIT']
    return(newhdu)

def interactive_matching(ra, dec, options):
    """
    Interactively i) match midas comps to the same source, ii) measure the EMU and/or GLASS largest-linear-size, iii) extract the coordinates of the radio core, iv) extract the coordinates of the host galaxy.

    parameters
    ----------

    """
    # Define the postage stamp cutout radius, emu half, and glass region.
    rad = options.cutout_size
    emu_half = get_emu_half(ra)
    glass_reg = get_glass_region(ra, dec)

    # now append the mosaic hdus to a list (hdu_list_0)
    hdu_list_0 = []
    print('Reading in mosaic images.')
    for i in range(0, len(options.mosaic_names)):
        mosaic_name = (options.mosaic_names[i]).replace('glass_reg', glass_reg).replace('emu_half', str(emu_half))
        path_to_mosaic = '{0}/{1}.fits'.format(options.mosaic_dir,mosaic_name)
        hdu0 = generate_cutout(options.survey_prefix[i], ra, dec, rad, path_to_mosaic=path_to_mosaic, return_hdu=True)
        print('    Succesfully created 2D cutout from {0}'.format(mosaic_name))
        hdu_list_0.append(hdu0)

    # need to convert lists into arrays so I can apply masks
    hdu_list_0 = np.array(hdu_list_0)
    mosaic_type = np.array(options.mosaic_type)
    survey_prefix = np.array(options.survey_prefix)
    mean_rms = np.array(options.mean_rms)
    contour_colour = np.array(options.contour_colour)

    # now do interactive things on the canvas.
    viking_hdu = (hdu_list_0[np.where(mosaic_type == 'o')][0])
    viking_imdata, viking_w = viking_hdu[0].data, wcs.WCS(viking_hdu[0].header, naxis=2)
    polypick = make_radio_overlay(hdu_list_0, mosaic_type, survey_prefix, mean_rms, np.array([3,5,10]), contour_colour, savename='test', reprojected=False, interactive=True)
    plt.show()

    if len(polypick.geomcent.x) and len(polypick.las.x):
        radio_center = Coords()
        LAS_points = Coords()
        # Use the latest 'key press event' for the geometric radio center
        radio_center.x, radio_center.y = viking_w.wcs_pix2world((polypick.geomcent.x)[len(polypick.geomcent.x) - 1],
                                                                (polypick.geomcent.y)[len(polypick.geomcent.x) - 1], 0)
        # Convert LAS points from pixel to world coords
        LAS_points.x, LAS_points.y = viking_w.wcs_pix2world(polypick.las.x,polypick.las.y, 0)
        new_cutout_size = (compute_angular_separation(LAS_points.x[0], LAS_points.y[0],
                                                      LAS_points.x[1], LAS_points.y[1])).arcsecond

        # Now go and measure the angular size and extract the radio core position from all surveys in options.lls and options.core
        print("Creating a new {2}'' cutout at RA = {0}, Dec = {1} ".format(float('%.4g' % (radio_center.x)),
                                                                           float('%.4g' % (radio_center.y)),
                                                                           float('%.3g' % (new_cutout_size))))
        measure_radio_properties = list(set(options.lls + options.core))
        lls_list, core_ra_list, core_dec_list, core_type_list = np.zeros([len(options.lls)]), \
                                                                np.zeros([len(options.core)]), np.zeros([len(options.core)]), []
        for i in range(0,len(measure_radio_properties)):
            prefix = options.lls[i]
            if prefix in survey_prefix:
                print('Extracting radio properties from {0}'.format(prefix))
                hdu = hdu_list_0[np.where(survey_prefix==prefix)][0]
                hdu = generate_cutout(prefix, radio_center.x, radio_center.y, new_cutout_size,
                                      input_hdu=hdu, return_hdu=True)
                lls_list[i], core_ra_list[i],  core_dec_list[i], core_type = \
                    extract_radio_properties(hdu, prefix, 1e-3*mean_rms[np.where(survey_prefix==prefix)], savedir=options.temp_dir)
                core_type_list.append(core_type)

        ## TODO: save the cutouts as you go.
        # Now go and interactively match the host galaxy(ies).
        new_cutout_size = 1.2*max(lls_list)
        print("Creating a new {2}'' cutout at RA = {0}, Dec = {1} ".format(float('%.4g' % (radio_center.x)),
                                                                           float('%.4g' % (radio_center.y)),
                                                                           float('%.3g' % (new_cutout_size))))
        # Now reproject radio images in hdu_list_0 and cutout based on measured LAS.
        reprojected_hdu = []
        for i in range(0,len(hdu_list_0)):
            hdu = hdu_list_0[i]
            # If radio image, reproject to viking header.
            if mosaic_type[i] == 'r':
                array, footprint = reproject_interp(hdu[0], viking_hdu[0].header)
                hdu = fits.PrimaryHDU(array)
                hdu = fits.HDUList([hdu])
                hdu = add_psf_header_keys(hdu, viking_hdu[0].header, hdu_list_0[i])
            reprojected_hdu.append(generate_cutout(survey_prefix[i], radio_center.x, radio_center.y, new_cutout_size,
                                  input_hdu=hdu, return_hdu=True))
        reprojected_hdu = np.array(reprojected_hdu)
        polypick = make_radio_overlay(reprojected_hdu, mosaic_type, survey_prefix, mean_rms, np.geomspace(5,100,5), contour_colour, savename='reproject',
                                      reprojected=True,galaxy_cat='Yes',interactive=True)
        plt.show()

    else:
        print('No interesting radio source at this position, skipping.')

    # measure_host()
        # if no 'button press' -> no host
        # if 'left click' -> GAMA23 host exists
        # if 'right click' -> non-GAMA23 host exists

    # match midas_cmps()
        # make cutout r=2.5*LAS. if LAS==False, r = 500''
        # 2-panel figure. Left: overlay, right: midas image w/ midas cat
        # call separate xmatch() function to xmatch user-input cmps with midas cat and write xmatched cmps to out_midas_cat. Include uuids.

    # check_for_blend()
        #if check_for_blend($surv) == True, check for blend

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prefix_chars='-')

    # Target cutout
    group0 = parser.add_argument_group("Input coords")
    group0.add_argument("--ra_target", dest='ra_target', type=float, default=None, help='Target RA in degrees.')
    group0.add_argument("--dec_target", dest='dec_target', type=float, default=None, help='Target Dec in degrees.')
    group0.add_argument("--coords_list", dest='coords_list', type=str, default=None,
                        help='Comma-separated .txt file containing list of RA and Dec in degrees.')
    group0.add_argument("--cutout_size", dest='cutout_size', type=float, default=500,
                        help='Size of initial cutout in arcseconds.')

    # Mosaics
    group1 = parser.add_argument_group('Mosaics')
    group1.add_argument("--mosaic_dir", dest='mosaic_dir', type=str, default='/mnt/e/Thesis/mosaics/',
                        help='Directory containing the mosaic fits images. len(mosaic_path) = N.')
    group1.add_argument("--mosaic_names", dest='mosaic_names', nargs='+', type=str, default=None,
                        help='List the names of each mosaic. len(mosaic_names) = N.')
    group1.add_argument("--survey_prefix", dest='survey_prefix', nargs='+', type=str, default=None,
                        help='Unique prefix for each survey mosaic, e.g. "glass5" for the 5.5 GHz GLASS survey. len(prefix) = N.')
    group1.add_argument("--mosaic_type", dest='mosaic_type', nargs='+', type=str, default=None,
                        help='"r" for radio, "o" for optical/near-IR. len(mosaic_type) = N.')
    group1.add_argument("--mean_rms", dest='mean_rms',  nargs='+', type=float, default=None,
                        help='Average rms, in mJy/beam, for each radio image. Use "-" for any non-radio image.')

    # Tables
    group2 = parser.add_argument_group('Tables')
    group2.add_argument("--table_dir", dest='table_dir', type=str, default=None,
                        help='Directory containing any radio/galaxy catalogs')
    group2.add_argument("--cat_names", dest='cat_names', nargs='+', type=str, default=None,
                        help='Name of catalogues associated with surveys')
    group2.add_argument("--cat_pref", dest='cat_pref', nargs='+', type=str, default=None,
                        help='')
    group2.add_argument("--cat_ra_name", dest='cat_ra_name', nargs='+', type=str, default=None,
                        help='Name of RA column for each table in "cat_names". Use "-" for blank entries.')
    group2.add_argument("--cat_dec_name", dest='cat_dec_name', nargs='+', type=str, default=None,
                        help='Name of Dec column for each table in "cat_names". Use "-" for blank entries.')

    # NVSS. Separate group given because the NVSS images are sourced from Skyview.
    group3 = parser.add_argument_group('NVSS')
    group3.add_argument("--nvss", dest='nvss', action='store_true', default=False,
                        help='(Optional). Whether to include NVSS into the cross matching process.')
    group3.add_argument("--nvss_cat", dest='nvss_cat', type=str, default=None,
                        help='Name of NVSS catalogue. Note, this must exist in the same directory as "cat_names".')
    group3.add_argument("-nvss_prefix", dest='nvss_prefix', type=str, default='nvss',
                        help='Prefix for NVSS image.')

    # Usage
    group8 = parser.add_argument_group('Usage')
    group8.add_argument("--lls", dest='lls', nargs='+', default=None,
                        help='Prefix of survey(s) from which to measure the largest linear size.')
    group8.add_argument("--core", dest='core', nargs='+', default=None,
                        help='Prefix of survey(s) from which to extract properties of the radio core.')

    # Output source table
    group4 = parser.add_argument_group('Output')
    group4.add_argument("--output_dir", dest='output_dir', type=str, default=os.getcwd().replace('Scripts', 'Output'),
                        help='Directory that outputs are written to.')
    group4.add_argument("--create_new", dest='create_new', action='store_true', default=False,
                        help='Whether to create a new blank source table.')
    group4.add_argument("--out_cat", dest='out_cat', type=str, default=None,
                        help='Name of output source catalogue. This must exist within "output_dir"/Tables/')
    group4.add_argument("--out_name", dest='out_name', type=str, default='Name',
                        help='Name of column in output source table containing the radio source name.')
    group4.add_argument("--out_ra", dest='out_ra', type=str, default='ra',
                        help='Name of column in output source table containing the right ascension.')
    group4.add_argument("--out_dec", dest='out_dec', type=str, default='dec',
                        help='Name of column in output source table containing the declination.')
    group4.add_argument("--out_morph", dest='out_morph', type=str, default='Morphology',
                        help='Name of column in output source table containing the radio morphology')
    group4.add_argument("--out_frclass", dest='out_frclass', type=str, default='FR_classification',
                        help='Name of column in output source table containing the FR classification.')
    group4.add_argument("--out_agnphase", dest='out_agnphase', type=str, default='AGN_phase',
                        help='Name of column in output source table containing the AGN phase.')
    group4.add_argument("--out_rgphase", dest='out_rgphase', type=str, default='RG_phase',
                        help='Name of column in output source table containing the radio galaxy phase.')
    group4.add_argument("--out_emu_lls", dest='out_emu_lls', type=str, default='emu_LLS',
                        help='Name of column in output source table containing the EMU largest linear size.')
    group4.add_argument("--out_glass_lls", dest='out_glass_lls', type=str, default='glass_LLS',
                        help='Name of column in output source table containing the GLASS largest linear size.')
    group4.add_argument("--out_emu_core_ra", dest='out_emu_core_ra', type=str, default='ra_core_887',
                        help='Name of column in output source table containing the EMU radio core right ascension. ')
    group4.add_argument("--out_emu_core_dec", dest='out_emu_core_dec', type=str, default='dec_core_887',
                        help='Name of column in output source table containing the EMU radio core declination')
    group4.add_argument("--out_glass_core_ra", dest='out_glass_core_ra', type=str, default='ra_core_5500',
                        help='Name of column in output source table containing the GLASS radio core right ascension. ')
    group4.add_argument("--out_glass_core_dec", dest='out_glass_core_dec', type=str, default='dec_core_5500',
                        help='Name of column in output source table containing the EMU radio core declination. ')
    group4.add_argument("--out_host_ra", dest='out_host_ra', type=str, default='ra_host',
                        help='Name of column in output source table containing the host galaxy right ascension ')
    group4.add_argument("--out_host_dec", dest='out_host_dec', type=str, default='dec_host',
                        help='Name of column in output source table containing the host galaxy declination')
    group4.add_argument("--host_in_gama_cat", dest='host_in_gama_cat', type=str, default='host_in_gama_cat',
                        help='Name of column in output source table which indicates whether the host galaxy is recorded in the GAMA23 photometry catalogue')
    group4.add_argument("--out_nmidas", dest='out_nmidas', type=str, default='Nmidas_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_ngmrt", dest='out_ngmrt', type=str, default='Ngmrt_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nemu", dest='out_nemu', type=str, default='Nemu_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nemu_lowres", dest='out_nemu_lowres', type=str, default='Nemu_lowres_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nglass", dest='out_nglass', type=str, default='Nglass_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nnvss", dest='out_nnvss', type=str, default='Nglass_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nglass_lowres", dest='out_nglass_lowres', type=str, default='Nglass_lowres_cmps',
                        help='Name of in the output source table.')
    group4.add_argument("--out_nhost", dest='out_nhost', type=str, default='Nhosts',
                        help='Name of in the output source table.')

    # Make low_res opts
    group5 = parser.add_argument_group('Make low res')
    group5.add_argument("--convolve_images", dest='convolve_images', nargs='+', type=str, default=None,
                        help='Prefixes of surveys/mosaics to be convolved')
    group5.add_argument("--extract_beam_from", dest='exctract_beam_from', type=str, default=None,
                        help='Low resolution image from which to extract synthesized beam properties')
    group5.add_argument("--beam_shape", dest='beam_shape', nargs='+', type=float, default=None,
                        help='Alternatively supply a synthesized beam for convolving.')
    group5.add_argument("--sourcefind_lowres", dest='sourcefind_lowres', action='store_true', default=True,
                        help='Run BANE and AEGEAN on low resolution images')

    # Cross matching options
    group6 = parser.add_argument_group('Output verbose')
    group6.add_argument("--new_xmatch_cat_verbose", dest='new_xmatch_cat_verbose', action='store_true', default=False,
                        help='Whether to create a new blank verbose table(s) containing cross matched components.')
    group6.add_argument("--xmatch_cat_prefix", dest='xmatch_cat_prefix', nargs='+', type=str, default=None,
                        help='Prefixes of surveys/mosaics to be cross matched. ')
    group6.add_argument("--out_xmatch_ra", dest='out_xmatch_ra', nargs='+', type=str, default=None,
                        help='Name of RA column for each table in "xmatch_cat_prefix". Use "-" for blank entries.')
    group6.add_argument("--out_xmatch_dec", dest='out_xmatch_dec', nargs='+', type=str, default=None,
                        help='Name of Dec column for each table in "xmatch_cat_prefix". Use "-" for blank entries.')

    # Extra
    group7 = parser.add_argument_group('Output verbose')
    group7.add_argument("--temp_dir", dest='temp_dir', default='/tmp/pipeline/', help='Temporary working directory.')
    group7.add_argument("--viking_cmap", dest='viking_cmap', default='gray_r', help='Colormap used to display any VIKING image.')
    group7.add_argument("--contour_colour", dest='contour_colour', nargs='+', type=str, default=None,
                        help='Contour colors used for each survey in "mosaic_names". ')


    options = parser.parse_args()

    if options.ra_target:
        interactive_matching(options.ra_target, options.dec_target, options)

    elif options.coords_list:
        path_to_coords_file = '/mnt/e/Thesis/Paper_II/Pipeline/{0}'.format(options.coords_list)                                 # FIX
        print("Reading in coordinates from {0}".format(path_to_coords_file))
        hdu_mosaics = []
        with open(path_to_coords_file) as f:
            lines = [line.rstrip() for line in f]
            for i in range(0, len(lines)):
                print('+ + + + + + + + + + + + + + + +')
                ra, dec = lines[i].split(',')
                ra, dec = float(ra), float(dec)
                print('Read in new coordinate: RA={0}, Dec={1}'.format(ra, dec))
                interactive_matching(ra, dec, options)