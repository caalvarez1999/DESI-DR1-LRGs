from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import astropy.coordinates as coords
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits, ascii
from astropy.table import Table, Column
from astroquery.sdss import SDSS

def load_spec(ID):
    def detect_spectrum_window(flux):
        """ Find the spectrum window excluding non positive side regions
        Input: 
            - flux 
        Output:
            - boolean mask with flux > 0 regions
        """
        il   = 0
        ir   = -1
        while True:
            if (flux[il] > 0):
                break
            il += 1
        while True:
            if (flux[ir] > 0):
                break
            ir -= 1        
        flag_window = np.ones_like(flux, dtype=bool)
        flag_window[:il] = 0
        flag_window[ir:] = 0
        return flag_window
    
    """ Loads spectroscopic data from file. """
    hdulist  = fits.open(ID)

    wave = 10.0**hdulist[1].data["loglam"]
    flux = hdulist[1].data["flux"]
    ferr = hdulist[1].data["ivar"]
    wdisp = hdulist[1].data["wdisp"]
    qual = np.zeros_like(ferr)+1

#    mask = detect_spectrum_window(flux)
    
    spectrum = [wave, 
                flux, 
                ferr,
                wdisp,
                qual]
                
    hdulist.close()

    return spectrum


def load_spec(ID):
    try:
        # Open the FITS file and load data
        with fits.open(ID) as hdulist:
            # Assuming "loglam" exists in the file and is in the first HDU
            wave = 10.0**hdulist[1].data["loglam"]
            flux = hdulist[1].data["flux"]
            dflux = hdulist[1].data["ivar"]
            wdisp = hdulist[1].data["wdisp"]
            qual = np.zeros_like(dflux)+1
            qual2 = 1
            return wave, flux, dflux, wdisp, qual, qual2
            
    except (OSError, TypeError) as e:
        print(f"Warning: Could not load file {ID} due to error: {e}")
        wave = np.arange(3700.0, 6700.5, 0.5)
        flux = np.zeros_like(wave)+1.0+np.random.rand()
        dflux = flux*0.1
        wdisp = np.zeros_like(wave)+1.0
        qual = np.zeros_like(dflux)+1
        qual2 = 0
        
        return wave, flux, dflux, wdisp, qual, qual2

