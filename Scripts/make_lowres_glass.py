# create low-res GLASS images and use BANE/AEGEAN to source find and extract fluxes. 

#!/usr/bin/python3

from tasks import get_glass_region
from tasks import get_emu_half

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table, Column

from create_source_new import generate_cutout

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess
Popen = subprocess.Popen


def bane_wrapper(fitsfile, grid, box):
    Popen("BANE {0} --grid {1} {1} --box {2} {2}".format(fitsfile, grid, box), shell=True).wait()

def aegean_wrapper(fitsfile):
    Popen("aegean --autoload --table={0} {0}".format(fitsfile), shell=True).wait()

def convol_wrapper(mirfile, outfile, beam):
    Popen("convol map={} out={} fwhm={},{} pa={} options=final".format(mirfile, outfile, beam[0], beam[1], beam[2]),  shell=True).wait()

def fitsin(fitsfile, mirfile):
    Popen("fits in={} out={} op=xyin".format(fitsfile, mirfile), shell=True).wait()
    
def fitsout(fitsfile, mirfile):
    Popen("fits in={} out={} op=xyout".format(mirfile, fitsfile), shell=True).wait()



ra, dec, size, survey_name = 349.130, -31.601, 600, 'GLASS'

glass_reg = get_glass_region(ra, dec)

# glass_hdu = generate_cutout(survey_name, ra, dec, size, path_to_mosaic='/mnt/e/Thesis/mosaics/G23_GLASS_5500MHz_reg{0}.fits'.format(glass_reg), savehdu=True)
fitsfile = '/mnt/e/temp/MIDAS_J231631-313604_GLASS_600.fits'
mirfile = '/mnt/e/temp/MIDAS_J231631-313604_GLASS_600.im'
outfile = '/mnt/e/temp/MIDAS_J231631-313604_GLASS_600_cv.im'
out_fitsfile = '/mnt/e/temp/MIDAS_J231631-313604_GLASS_600_cv.fits'
fitsin(fitsfile, mirfile)
convol_wrapper(mirfile, outfile, [60,60,0])
fitsout(out_fitsfile, outfile)

hdu = fits.open(out_fitsfile)
header = hdu[0].header
bmaj = header['BMAJ']
cd = header['CDELT2']
grid = int(20.*bmaj/cd)
box = int(5.*grid)
bane_wrapper(out_fitsfile, grid, box)
aegean_wrapper(out_fitsfile)