from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.visualization import make_lupton_rgb

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord

from astropy.visualization import simple_norm

from reproject import reproject_interp

import numpy as np
from astropy.nddata import Cutout2D

def generate_cutout(survey_name, ra, dec, size, path_to_mosaic=None, input_hdu=None, target_name=False, savehdu=False,return_hdu=False):
    if path_to_mosaic:
        hdu = fits.open(path_to_mosaic)
    else:
        hdu = input_hdu
    image_data = np.squeeze(hdu[0].data)
    image_header = hdu[0].header
    image_header['CUNIT1'] = 'deg'
    image_header['CUNIT2'] = 'deg'
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
            target_name = create_IAU_name(ra, dec, 'MIDAS')
        savename = '{0}_{1}_{2}.fits'.format(target_name, survey_name, size)
        # validate_path('{0}/{1}'.format(get_paths('G23_cutouts_first_pass'), target_name), create_path=True)
        hdu.writeto('{0}/{1}/{2}'.format('/mnt/e/Thesis/Paper_II/Pipeline/Sources/', target_name, savename), overwrite=True)
    if return_hdu:
        return(hdu)

def get_skyview(survey_name, ra, dec, size):
    from astroquery.skyview import SkyView
    pstg_stamp = SkyView.get_images("{0},{1}".format(ra, dec), radius=size * u.arcsecond, survey=survey_name)[0]
    image_header = pstg_stamp[0].header
    # print(image_header)
    if survey_name == 'NVSS':
        image_header['BMAJ'] = 0.0125
        image_header['BMIN'] = 0.0125
        image_header['BPA'] = 0
    pstg_stamp.writeto('{2}_{0}_{1}_{3}.fits'.format(ra, dec, survey_name, size),overwrite=True)
    return(pstg_stamp)

get_skyview("DSS2 Blue", 204.253958, -29.865417, 20*60)
c = SkyCoord('04:03:54.3 -43:20:56', unit=(u.hourangle, u.deg))
print(c.ra.degree)
print(c.dec.degree)
# ra, dec, coadd_object_ID, xs, ysize
# image_dir = '/mnt/e/temp/swarp/'
# i_hdu = fits.open(image_dir+'i_band.fits')
# j_hdu = fits.open(image_dir+'j_band.fits')
# k_hdu = fits.open(image_dir+'k_band.fits')

# 11.9208,-25.2669,DES0047-2458,12,12
# 11.8326,-25.3302,DES0046-2541,12,12

# j_hdu[0].header['CUNIT1'] ='deg'
# j_hdu[0].header['CUNIT2'] ='deg'


# hdu = generate_cutout('test', 204.25302, -29.866, 15*60, path_to_mosaic=image_dir+'j_band.fits', return_hdu=True)

# # # # i_hdu[0].header['CUNIT1'] ='deg'
# # # # i_hdu[0].header['CUNIT2'] ='deg'

# i_imdata, footprint = reproject_interp(i_hdu, hdu[0].header)

# # print(i_imdata.an

# w = wcs.WCS(hdu[0].header, naxis=2)

# fig = plt.figure(figsize=(10,10))
# ax = fig.add_axes([0.1,0.1,0.8,0.8], projection=w)
# ax.imshow(i_imdata)
# plt.show()