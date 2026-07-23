#!/usr/bin/env python
# coding: utf-8

# We import the needed modules:

# General loads

import numpy as np
import pandas as pd
import astropy.coordinates as coords
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from astropy.table import Table, Column


# Local loads

from pylick.analysis import Galaxy, Catalog
from .corrections import *
import os
import time as T
from pylick.loadspec import load_spec
from .emcee_lick import I_name, I_name29, I_name30, I_nameall, I_nameopen, I_number, I_number29, I_number30, I_numberall, I_numberopen


# Spectro-photometric features and names (ordered by z)

#resdir = './spectra/galxs_ours/res_2003.ecsv'
#res_ = Table.read(resdir)
#namesdir = './spectra/galxs_ours/allnames.npy'
#names_ = np.load(namesdir)


# List of indices that might be useful

list_full = np.arange(0, 30, 1)
list_25 = np.arange(41, 66, 1)
list_27 = np.asarray([i+37 if i<2 else i+37+2 for i in range(27)])


# measure_base calculates lick indices in three different modalities (RAW,
# IDS, MILES) depending on whether a smoothing to degrade the resolution is
# applied or not. The IDS option degrades the spectra to the resolution of
# the Lick-IDS system in which the resolution is wavelength dependent. The MILES
# option degrades them to a constant resolution corresponding to a FWHM dispersion
# of 2.5 \AA.


def SDSS_SN(res, ruta = '../spectra/SDSSFULL/full/', lamb = 3700.0, lamr = 6500.0):
    out = np.zeros(len(res)) #SN, dimension = number of objects
    dout = np.zeros(len(res)) #dSN, " " " "
    for igal in range(len(res)):
        if igal%1000 == 0 or igal == len(res)-1:
            print(str(round(igal/len(res)*100, 2))+'%')
        name = res['name'][igal]
        ID = ruta+str(name) #where to take the spectrum from
        w, f, df, wdisp, qual = load_spec(ID) #function that gives us w, f... from the fits file
        w = C_toair(w) #convert the wavelength to air
        zgal = res['z1'][igal]
        w_RF = w/(1.0+zgal) #restframe wavelengths
        
        #we select the central window (3700A to 6500A)
        cut = np.intersect1d(np.where(w_RF>lamb)[0], np.where(w_RF<lamr)[0])
        cutcut = np.intersect1d(cut, np.where(df>0))
        W = w_RF[cutcut]
        F = f[cutcut]
        dF = df[cutcut]
        out[igal] = np.median(F/dF)
        dout[igal] = np.std(F/dF)
    
    return out, dout



def measure_indices(w, f, df, SN=False, truesampling=False, nsampling=100, index_list=list_full):
    # Calcular SN si no está dado
    if not np.any(SN):
        SN = f / df

    ID = '_'
    n_indices = len(index_list)

    # Caso con NaNs → salida vacía
    if np.any(np.isnan(f)):
        outind = np.full(n_indices, np.nan)
        outdind_perc = np.full(n_indices, np.nan)
        outdind_std = np.full(n_indices, np.nan)
        outSNR = np.full(n_indices, np.nan)
        return outind, outdind_perc, outdind_std, outSNR

    # Medición base (sin muestreo)
    ind = Galaxy(ID, index_list, spec_wave=w, spec_flux=f, spec_err=df,
                 spec_mask=None, meas_method='int', z=0.0)
    
    outind = ind.vals
    outdind_perc = ind.errs
    outdind_std = ind.errs
    outSNR = ind.SNR

    # Muestreo Monte Carlo (si se solicita)
    if truesampling:
        sigma = np.abs(f / SN) + 1e-14  # ruido por pixel (precalculado)
        
        # Generar todos los espectros perturbados de una sola vez
        noise = np.random.normal(0, sigma, size=(nsampling, len(f)))

        local = np.empty((nsampling, n_indices))  # preasignación para eficiencia

        for i in range(nsampling):
            f_ = f + noise[i]
            ind_ = Galaxy(ID, index_list, spec_wave=w, spec_flux=f_, spec_err=df,
                          spec_mask=None, meas_method='int', z=0.0)
            local[i] = ind_.vals

        # Cálculos estadísticos
        outdind_perc = 0.5 * (np.percentile(local, 84, axis=0) - 
                              np.percentile(local, 16, axis=0))
        outdind_std = np.std(local, axis=0)

    return outind, outdind_perc, outdind_std, outSNR




def do(lalista, res, nameout, Nsampling = 100, ruta = '../spectra/SDSSFULL/full/', ognicuanto = 100, resolution = False, sigmaresolution = 1.5, snmedian_ = False):
    ii_ = np.zeros((len(res), len(lalista)))
    ddii_ = np.zeros((len(res), len(lalista), 2))
    iiSNR_ = np.zeros((len(res), len(lalista)))
    names = res['name']
    z1 = res['z1']
    if snmedian_:
        snmedian = res['snmedian']
    
    ladlista = ['d'+i for i in lalista]
    lalistaSNR = ['SNR_'+i for i in lalista]
    iiilalista = I_numberopen(lalista)
    
    
    t0 = T.time()
    tORIGINAL = t0
    
    QUAL2 = np.zeros(len(res))
    
    for igal in range(len(res)):
        if igal%ognicuanto == 0:
            t1 = T.time()
            print(str(round(igal*100.0/len(res), 1))+'% of the galaxies was analysed')
            dTiempo = t1 - t0
            print(str(round(dTiempo/1.0, 0))+' secs have passed from gal '+str(igal-ognicuanto)+' to gal '+str(igal))
            dTotTiempo = t1 - tORIGINAL
            if igal > 0:
                print(str(round(dTotTiempo/igal*len(res)/3600*(1.0-igal/len(res)*1.0), 3))+' hours until the end')
            t0 = t1
        
        ID = ruta+str(names[igal])
        w, f, df, wdisp, qual, qual2 = load_spec(ID)
        
        s = wdisp*np.append(np.diff(w), np.diff(w)[-1])
        
        w = C_toair(w)
        zgal = z1[igal]
        
        W_RF = w/(1.0+zgal)
        
        
        #cut = np.intersect1d(np.where(W_RF >= 3700.0)[0], np.where(W_RF <= 6700.0)[0])
        cut = (W_RF >= 3600.0) & (W_RF <= 6700.0)
        
        laW = W_RF[cut]
        laF = f[cut]
        ladF = df[cut]
        laS = s[cut]
        
        if np.sum(abs(ladF)) < 0.1:
            qual2 = 0
        QUAL2[igal] = qual2
        
        
        if resolution:
            laW, laF, ladF = C_resol(laW, 0.0, laF, ladF, laS, sigmaresolution)
        
        
        if qual2 > 0:
            if snmedian_:
                ii_[igal, :], ddii_[igal, :, 0], ddii_[igal, :, 1], iiSNR_[igal, :] = measure_indices(laW, laF, ladF, SN = snmedian[igal], truesampling = True, nsampling = Nsampling, index_list = iiilalista)
            else:
                ii_[igal, :], ddii_[igal, :, 0], ddii_[igal, :, 1], iiSNR_[igal, :] = measure_indices(laW, laF, ladF, SN = False, truesampling = True, nsampling = Nsampling, index_list = iiilalista)
        
        else:
            ii_[igal, :], ddii_[igal, :, 0], ddii_[igal, :, 1], iiSNR_[igal, :] = np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan
            
    
    OUT = res.copy()
    for i in range(len(lalista)):
        OUT[lalista[i]] = ii_[:, i]
        OUT[ladlista[i]] = ddii_[:, i, 0]
        
        #OUT[lalista[i]] = ii_[:, i]
        #OUT[ladlista[i]+'_perc'] = ddii_[:, i, 0]
        #OUT[ladlista[i]+'_std'] = ddii_[:, i, 1]
        #OUT[lalistaSNR[i]+'Lick'] = iiSNR_[:, i]
        #OUT[lalistaSNR[i]+'_perc'] = abs(ii_[:, i]/ddii_[:, i, 0])
        #OUT[lalistaSNR[i]+'_std'] = abs(ii_[:, i]/ddii_[:, i, 1])
    OUT['QUALITY'] = QUAL2
    
    OUT.write(nameout, overwrite=True)




def doGAMA(lalista, res, nameout, Nsampling = 100, ruta = '../spectra/GAMA/', ognicuanto = 100, resolution = False, sigmaresolution = 2.5):
    
    ii_ = np.zeros((len(res), len(lalista)))
    ddii_ = np.zeros((len(res), len(lalista), 2))
    iiSNR_ = np.zeros((len(res), len(lalista)))
    
    ladlista = ['d'+i for i in lalista]
    lalistaSNR = ['SNR_'+i for i in lalista]
    iiilalista = I_numberopen(lalista)
    
    t0 = T.time()
    tORIGINAL = t0
    
    QUAL = np.zeros(len(res))
    
    for isel in range(len(res)):
        
        if isel%ognicuanto == 0:
            t1 = T.time()
            print(str(round(isel*100.0/len(res), 1))+'% of the galaxies was analysed')
            dTiempo = t1 - t0
            print(str(round(dTiempo/1.0, 0))+' secs have passed from gal '+str(isel-ognicuanto)+' to gal '+str(isel))
            dTotTiempo = t1 - tORIGINAL
            if isel > 0:
                print(str(round(dTotTiempo/isel*len(res)/3600*(1.0-isel/len(res)*1.0), 3))+' hours until the end')
            t0 = t1
        
        
        url = res['URL'].iloc[isel].decode('utf-8')
        nombre = url.rsplit('/', 1)[-1]
        redz = res['z1'].iloc[isel]
        hdu = fits.open(ruta+nombre)[0]
        
        hdr = hdu.header
        data = hdu.data
        gratid = hdr.get("GRATID", "").strip()
        
        # Flux (row 1)
        flux = data[0, :]
        valid = np.isfinite(flux)
        F = flux[valid]
        
        # Error (row 2)
        dF = data[1, valid]
        
        # Build wavelength array
        npix = hdr['NAXIS1']
        crval = hdr['CRVAL1']
        cdelt = hdr['CD1_1']
        crpix = hdr['CRPIX1']
        
        pixels = np.arange(npix)
        wav = crval + (pixels + 1 - crpix) * cdelt
        
        # resoluciones típicas (FWHM en Angstrom)
        if gratid == "385R 580V":
            blue_fwhm = 3.2
            red_fwhm = 5.3
        else:
            blue_fwhm = red_fwhm = np.nan  # desconocido
        
        # convertir a sigma si quieres
        blue_sigma = blue_fwhm / 2.355
        red_sigma = red_fwhm / 2.355
        
        # crear array: lambda<5700 blue, lambda>=5700 red
        sigma = np.where(wav < 5700, blue_sigma, red_sigma)
        S = sigma[valid]
        W = wav[valid]/(1.0+redz)
        
        cut = (W >= 3600.0) & (W <= 6700.0)
        laW = W[cut]
        laF = F[cut]
        ladF = dF[cut]
        laS = S[cut]
        
        qual = 1
        if np.sum(abs(ladF)) < 0.01:
            qual = 0
        
        if resolution:
            laW, laF, ladF = C_resol(laW, 0.0, laF, ladF, laS, sigmaresolution)
        
        
        if qual > 0:
            ii_[isel, :], ddii_[isel, :, 0], ddii_[isel, :, 1], iiSNR_[isel, :] = measure_indices(laW, laF, ladF, SN = False, truesampling = True, nsampling = Nsampling, index_list = iiilalista)
        
        else:
            ii_[isel, :], ddii_[isel, :, 0], ddii_[isel, :, 1], iiSNR_[isel, :] = np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan, np.zeros(len(iiilalista))+np.nan
            print('ojo NaNs everywhere')
        
        QUAL[isel] = qual
    
    OUT = res.copy()
    for i in range(len(lalista)):
        OUT[lalista[i]] = ii_[:, i]
        OUT[ladlista[i]] = ddii_[:, i, 0]
    OUT['qualitySPEC'] = QUAL
    
    Table.from_pandas(OUT).write(nameout, overwrite=True)



def doDESI(lalista, res, nameout, Nsampling = 100, ognicuanto = 100, resolution = False, sigmaresolution = 1.0):
    ii_ = np.zeros((len(res), len(lalista)))
    ddii_ = np.zeros((len(res), len(lalista), 2))
    iiSNR_ = np.zeros((len(res), len(lalista)))
    z1 = res['redshift']
    import time as T
    
    ladlista = ['d'+i for i in lalista]
    lalistaSNR = ['SNR_'+i for i in lalista]
    iiilalista = I_numberopen(lalista)
    snmedian = np.zeros(len(res))
    
    t0 = T.time()
    tORIGINAL = t0
    nmal = 0
    
    QUAL2 = np.zeros(len(res))
    
    for igal in range(len(res)):
        if igal%ognicuanto == 0:
            t1 = T.time()
            print(str(round(igal*100.0/len(res), 1))+'% of the galaxies was analysed')
            dTiempo = t1 - t0
            print(str(round(dTiempo/1.0, 0))+' secs have passed from gal '+str(igal-ognicuanto)+' to gal '+str(igal))
            dTotTiempo = t1 - tORIGINAL
            if igal > 0:
                print(str(round(dTotTiempo/igal*len(res)/3600*(1.0-igal/len(res)*1.0), 3))+' hours until the end')
            t0 = t1
        
        w, f, df, s, m = res['wavelength'][igal], res['flux'][igal], res['ivar'][igal] **(-0.5), res['wave_sigma'][igal], res['mask'][igal]
        valid = (m == 0)
        w, f, df, s = w[valid], f[valid], df[valid], s[valid]
        
        snmedian[igal] = np.median(f/df)

        w = C_toair(w)
        zgal = z1[igal]
        
        W_RF = w/(1.0+zgal)
        
        
        #cut = np.intersect1d(np.where(W_RF >= 3700.0)[0], np.where(W_RF <= 6700.0)[0])
        cut = (W_RF >= 3600.0) & (W_RF <= 6700.0)
        
        laW = W_RF[cut]
        laF = f[cut]
        ladF = df[cut]
        laS = s[cut]
        
        qual2 = 1
        if np.sum(abs(ladF)) < 0.1 or len(laW) < 0.05*len(m) or len(w) < 0.1*len(m):
            qual2 = 0
            nmal += 1
            print('mal espectro numero '+str(nmal))
        QUAL2[igal] = qual2
        
        
        if resolution:
            laW, laF, ladF = C_resol(laW, 0.0, laF, ladF, laS, sigmaresolution)
        
        
        if qual2 > 0:
            try:
                ii_[igal, :], ddii_[igal, :, 0], ddii_[igal, :, 1], iiSNR_[igal, :] = measure_indices(
                                                                                                        laW, laF, ladF,
                                                                                                        SN=snmedian[igal],
                                                                                                        truesampling=True,
                                                                                                        nsampling=Nsampling,
                                                                                                        index_list=iiilalista
                                                                                                      )
            except Exception as e:
                # Si algo falla en measure_indices o los arrays están vacíos
                nmal += 1
                print(f"⚠️ Espectro {igal} fallido ({str(e)[:60]}...) → asignando NaN")
                ii_[igal, :] = np.zeros(len(iiilalista))+np.nan
                ddii_[igal, :, 0] = np.zeros(len(iiilalista))+np.nan
                ddii_[igal, :, 1] = np.zeros(len(iiilalista))+np.nan
                iiSNR_[igal, :] = np.zeros(len(iiilalista))+np.nan
                QUAL2[igal] = 0
        else:
            ii_[igal, :] = np.zeros(len(iiilalista))+np.nan
            ddii_[igal, :, 0] = np.zeros(len(iiilalista))+np.nan
            ddii_[igal, :, 1] = np.zeros(len(iiilalista))+np.nan
            iiSNR_[igal, :] = np.zeros(len(iiilalista))+np.nan
            QUAL2[igal] = 0
            print("Espectro marcado como malo (qual2 = 0)")
    
    OUT = Table()
    OUT['targetid'] = res['targetid']
    OUT['specid'] = res['specid']
    OUT['z1'] = res['redshift']
    OUT['snmedian'] = snmedian
    OUT['qual'] = QUAL2
    for i in range(len(lalista)):
        OUT[lalista[i]] = ii_[:, i]
        OUT[ladlista[i]] = ddii_[:, i, 0]
    
    OUT.write(nameout, overwrite=True)



