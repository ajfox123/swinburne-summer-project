import healpy as hp
import numpy as np
from astropy.io import fits
import reproject as rp
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

NSIDE=1024
COORD_SYS='C'
IMAGE_FOLDER = "/data/vast-survey/pilot/EPOCH08/COMBINED/STOKESI_IMAGES"
FOOTPRINT_FOLDER = 'temp_footprints'
NPIX= hp.nside2npix(NSIDE)

def build_footprints():
    for askap_image in tqdm(os.listdir(IMAGE_FOLDER)):
        if not askap_image.endswith('.fits'):
            continue
    
        footprint_file = os.path.join(FOOTPRINT_FOLDER, askap_image)
        if os.path.isfile(footprint_file):
            continue
            
        with fits.open(os.path.join(IMAGE_FOLDER, askap_image)) as fitsfile:
            hdu = fitsfile[0]
            askap_header = hdu.header
            askap_data = hdu.data#[0,0,:,:]
            askap_wcs = WCS(askap_header, naxis=2)
            askap_footprint = np.ones(askap_data.shape)
            askap_footprint[np.isnan(askap_data)] = 0
            out_, footprint_reproject = rp.reproject_to_healpix((askap_footprint,askap_wcs),COORD_SYS,nside=NSIDE)
            hp.write_map(footprint_file, footprint_reproject)

def combine_footprints():
    full_footprint = np.zeros(NPIX)
    
    for footprintfile in tqdm(os.listdir(FOOTPRINT_FOLDER)):
        footprint = hp.read_map(os.path.join(FOOTPRINT_FOLDER, footprintfile))
        full_footprint[np.where(footprint)] = 1.0
        
    hp.write_map('full_VAST_footprint_{}.fits'.format(NSIDE), full_footprint)
    hp.mollview(full_footprint)
    plt.show()
      
if __name__ == '__main__':
    build_footprints()
    combine_footprints()