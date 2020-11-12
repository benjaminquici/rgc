import numpy as np
from scipy.optimize import curve_fit 

from astropy.io import fits
from astropy import wcs
from astropy import units as u

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt

from astropy.visualization import simple_norm

from reproject import reproject_interp

from astropy.visualization.wcsaxes import SphericalCircle

import matplotlib.patches as patches

def powlaw(freq,a,alpha): # defining powlaw as S = a*nu^-alpha. Important to have x value first in definition of function.
    return a*(freq**(-alpha))

def powlaw_model(freq,flux,fluxerr,plotting_freq):
	freq,flux,fluxerr=np.array(freq),np.array(flux),np.array(fluxerr)
	popt,pcov = curve_fit(powlaw,freq,flux,sigma=fluxerr)
	perr = np.sqrt(np.diag(pcov))
	fit = powlaw(plotting_freq, *popt)
	print('y-intercept:'+str(popt[0]))
	return(fit,popt[1],perr[1])

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
        # print('Beam properties not found, must be VIKING image.')
        pass
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

def get_skyview(survey_name, ra, dec, size):
    from astroquery.skyview import SkyView
    hdu = SkyView.get_images("{0},{1}".format(ra, dec), radius=size * u.arcsecond, survey=survey_name)[0]
    image_header = hdu[0].header
    return(hdu)

def hide_ax(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_facecolor('white')
    ax.patch.set_alpha(0.75)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    return(ax)
"""
Making plot 1
"""
mosaic_dir = '/mnt/e/Thesis/mosaics/'

ra, dec, size = 8.521, -66.66, 2000
atlbs_lowres_hdu = generate_cutout('atlbs_lr', ra, dec, size, path_to_mosaic='{0}/Low_resolution_region_A.fits'.format(mosaic_dir), return_hdu=True)
atlbs_highres_hdu = generate_cutout('atlbs_hr', ra, dec, size, path_to_mosaic='{0}/High_resolution_region_A.fits'.format(mosaic_dir), return_hdu=True)
gleam_hdu = get_skyview('GLEAM 170-231 MHz', ra, dec, size)
sumss_hdu = get_skyview('SUMSS 843 MHz', ra, dec, size)
wise_hdu = get_skyview('WISE 3.4', ra, dec, size)

fig = plt.figure(figsize=(17,10))

imdata_atlbs_lr, w_atlbs_lr = atlbs_lowres_hdu[0].data, wcs.WCS(atlbs_lowres_hdu[0].header,naxis=2)
imdata_atlbs_hr, w_atlbs_hr = atlbs_highres_hdu[0].data, wcs.WCS(atlbs_highres_hdu[0].header,naxis=2)
imdata_gleam, w_gleam = gleam_hdu[0].data, wcs.WCS(gleam_hdu[0].header,naxis=2)
imdata_sumss, w_sumss = sumss_hdu[0].data, wcs.WCS(sumss_hdu[0].header,naxis=2)
imdata_wise, w_wise = wise_hdu[0].data, wcs.WCS(wise_hdu[0].header,naxis=2)

height, width = 0.43, 0.28
ax0 = fig.add_axes([0.08, 0.53, width, height], projection=w_gleam)
ax00 = fig.add_axes([0.1, 0.905, 0.5*width, 0.1*height])
hide_ax(ax00)
ax00.text(0.1,0.3, 'GLEAM 200 MHz', fontsize=18)

ax1 = fig.add_axes([0.38, 0.53, width, height], projection=w_atlbs_lr)
ax11 = fig.add_axes([0.4, 0.905, 0.8*width, 0.1*height])
hide_ax(ax11)
ax11.text(0.1,0.3, 'ATLBS (low-res) 1.4 GHz', fontsize=18)

ax2 = fig.add_axes([0.68, 0.53, width, height], projection=w_wise)
ax22 = fig.add_axes([0.7, 0.905, 0.2*width, 0.1*height])
hide_ax(ax22)
ax22.text(0.1,0.3, 'WISE', fontsize=18)

ax3 = fig.add_axes([0.08, 0.09, width, height], projection=w_sumss)
ax33 = fig.add_axes([0.1, 0.46, 0.5*width, 0.1*height])
hide_ax(ax33)
ax33.text(0.1,0.3, 'SUMSS 843 MHz', fontsize=18)

ax4 = fig.add_axes([0.38, 0.09, width, height], projection=w_atlbs_lr)
ax44 = fig.add_axes([0.4, 0.46, 0.8*width, 0.1*height])
hide_ax(ax44)
ax44.text(0.1,0.3, 'ATLBS (high-res) 1.4 GHz', fontsize=18)

ax5 = fig.add_axes([0.7, 0.09, width-0.02, height-0.01])
norm = simple_norm(imdata_gleam, percent=99.5)
ax0.imshow(imdata_gleam, cmap='inferno',norm=norm)
norm = simple_norm(imdata_atlbs_lr, percent=99.5)
ax1.imshow(imdata_atlbs_lr, cmap='inferno',norm=norm)
r1 = SphericalCircle((8.4732939*u.deg, -66.5745355*u.deg), 216.381*u.arcsecond, ls='-.', edgecolor='white', facecolor='none', transform=ax1.get_transform('fk5'))
r2 = SphericalCircle((8.5480254*u.deg, -66.7383874*u.deg), 216.381*u.arcsecond, ls='-.', edgecolor='white', facecolor='none', transform=ax1.get_transform('fk5'))
ax1.add_patch(r2)
ax1.add_patch(r1)
ax1.text(8.33602, -66.5361, 'Lobe (N)', transform=ax1.get_transform('fk5'), fontsize=13, color='white')
ax1.text(8.38624, -66.7651, 'Lobe (S)', transform=ax1.get_transform('fk5'), fontsize=13, color='white')
norm = simple_norm(imdata_wise, percent=99.5)
ax2.imshow(imdata_wise, cmap='cubehelix',norm=norm)
# imdata_gleam, footprint = reproject_interp(gleam_hdu[0], wise_hdu[0].header)
# ax2.contour(imdata_gleam, levels=np.geomspace(5,100,3)*7 * 1e-3, colors='magenta', lw=0.5)
imdata_atlbs_lr, footprint = reproject_interp(atlbs_lowres_hdu[0], wise_hdu[0].header)
ax2.contour(imdata_atlbs_lr, levels=np.geomspace(5,100,5)*100 * 1e-6, colors='white', lw=1)
norm = simple_norm(imdata_sumss, percent=99.5)
ax3.imshow(imdata_sumss, cmap='inferno',norm=norm)
norm = simple_norm(imdata_atlbs_hr, percent=99.5)
ax4.imshow(imdata_atlbs_hr, cmap='inferno',norm=norm)
ax5.set_xlabel('Flux density (mJy)', fontsize=15)
ax5.set_ylabel('Frequency (MHz)', fontsize=15)
flux =     [0.67616856,	0.56138563,	0.5285633,	0.5195936,	0.5345396,	0.43678287,	0.43036097,	0.3443639, 	0.37040266, 0.32250196, 0.29184285,	0.31681442,	0.26817632,	0.31625265,	0.2709392, 	0.24196316, 0.23580118, 0.21931256, 0.2594507,	0.20121422, 0.0658*0.4]
err_flux = [0.125143235, 0.087789446, 0.08036295, 0.08411249, 0.046199687, 0.038078204, 0.035572737, 0.035865434, 0.03179306, 0.02948395, 0.027866803, 0.02806203, 0.03391127, 0.031025216, 0.03297701, 0.03457455, 0.041161492, 0.041089114, 0.043635797, 	0.04155100586, 0.5*0.006]
flux_2 = [0.82615054,0.5695683,	0.88496935, 0.6486915, 	0.6991511, 	0.63993436, 0.5994822, 	0.5079539, 0.47598338, 	0.4158937, 	0.488144, 0.48346546, 0.38898316, 0.35547134, 0.36118504, 0.3427996, 0.2574919, 0.27846906, 0.45287484, 0.26566324, 0.0658*0.6]
err_flux_2 = [0.12155395, 0.09108733, 0.08674075, 0.0882207, 0.048740778, 0.03991658, 0.0366194, 0.037787374, 0.033069137, 0.03093353, 0.029237557, 0.029901056, 0.035164967, 0.03227255, 0.034097496, 0.037211772, 0.042341392, 0.043468114, 0.045385703, 0.045460813, 0.5*0.006]
freq = [76, 84, 92, 99, 107, 115, 122, 130, 143, 151, 158, 166, 174, 181, 189, 197, 204, 212, 220, 227, 1400]
flux = 1000*np.array(flux)
err_flux = 1000*np.array(err_flux)
flux_2 = 1000*np.array(flux_2)
err_flux_2 = 1000*np.array(err_flux_2)
plotting_freq = np.geomspace(10, 3000, 100)
fit, alpha, err_alpha = powlaw_model(freq,flux,err_flux,plotting_freq)
fit_2, alpha_2, err_alpha_2 = powlaw_model(freq,flux_2,err_flux_2,plotting_freq)
alpha = float('%.3g' %alpha)
alpha_2 = float('%.3g' %alpha_2)

ax5.errorbar(freq,flux,xerr=0,yerr=err_flux,color='black',capsize=3,linestyle='None',hold=True,fmt='none',alpha=0.75)
ax5.errorbar(freq,flux_2,xerr=0,yerr=err_flux_2,color='black',capsize=3,linestyle='None',hold=True,fmt='none',alpha=0.75)
ax5.scatter(freq,flux, marker='.', s=100,c='C0', label='Lobe (N): data')
ax5.scatter(freq,flux_2, marker='.', s=100,c='C1', label='Lobe (S): data')
ax5.plot(plotting_freq, fit, c='C0', ls=':', lw=1, label='Lobe (N): fit. $\\alpha = -{0}$'.format(alpha))
ax5.plot(plotting_freq, fit_2, c='C1', ls=':', lw=1, label='Lobe (S): fit. $\\alpha = -{0}$'.format(alpha_2))
ax5.set_yscale('log')
ax5.set_xscale('log')
ax5.legend(loc='upper right')
rect = patches.Rectangle((70,2000),210,-2000,linewidth=1,edgecolor='None',facecolor=[0,1,0,0.1])
ax5.add_patch(rect)
ax5.text(100,80, 'GLEAM', fontsize=15)
ax5.tick_params(axis='both',which='both',labelsize=15)
ax5.xaxis.set_tick_params(direction='in',which='minor',length=4,width=1.5)
ax5.xaxis.set_tick_params(direction='in',which='major',length=10,width=2,pad=10)
ax5.yaxis.set_tick_params(direction='in',which='minor',length=4,width=1.5)
ax5.yaxis.set_tick_params(direction='in',which='major',length=10,pad=10,width=2)
ax5.set_xlim(50,2500)
ax5.set_ylim(20, 1500)
# norm = simple_norm(imdata_atlbs_hr, percent=99.5)
# ax1.imshow(imdata_atlbs_hr,cmap='inferno', norm=norm)

axes = [ax0, ax1, ax2, ax3, ax4]#, ax5]
latleft = [True, False, False, True, False]#, False]
latright = [False, False, False, False, False]#, False]
lonbottom = [False, False, False, True, True]#, True]
lontop = [False, False, False, False, False]#, False]

for ax, lal, lar, lob, lot in zip(axes,latleft,latright,lonbottom,lontop):
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

# plt.show()
plt.savefig('/mnt/e/Thesis/choppa.png', dpi=200)