# -*- coding: utf-8 -*-
"""
CARLOS ALONSO ÁLVAREZ
Correction Functions
"""

#SOME IMPORTS THAT SHOULD BE GENERALLY OF GOOD USE
from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import astropy.coordinates as coords
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd

from astropy.io import fits, ascii
from astropy.table import Table, Column
from astroquery.sdss import SDSS
import aux.emcee_lick as mclick

from pylick.analysis import Galaxy, Catalog
#------------------------------------------------------------------------------

#INTENTIONALLY LEFT BLANK

#------------------------------------------------------------------------------

#VELOCITY DISPERSION CORRECTION:
#       

#------------------------------------------------------------------------------

#APERTURE CORRECTION:
#       Function that uses the aperture relation in CLEMENS+06 (equation 2)
#   to modify the incides obtained after other corrections have been applied,
#   it accounts for the effect of the aperture of the spectrograph
#   when taking the spectrum.



def C_toair(lambdas):
    #D. C. Morton 2000, The Astrophysical Journal Supplement Series, 130:403-436, 2000 October
    
    s = 1.0e4/lambdas
    n = 1+0.0000834254+0.02406147/(130.0-s**2.0)+0.00015998/(38.9-s**2.0)
    return lambdas/n

def C_tovac(lambdas):
    #D. C. Morton 2000, The Astrophysical Journal Supplement Series, 130:403-436, 2000 October
    
    s = 1.0e4/lambdas
    n = 1+0.00008336624212083+0.02408926869968/(130.1065924522-s**2.0)+0.0001599740894897/(38.92568793293-s**2.0)
    return lambdas*n

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

#RESOLUTION CORRECTION:
#       Function that takes in an array or list of wavelengths,
#   an array of fluxes and an array of the dispersion of the spectrograph (the 
#   latter in units of pixels).
#
#       As an output it provides an array of wavelengths in the resframe (or
#   provided redshift), an array of fluxes and an array of fluxes corrected by
#   the dispersion of the model; in case the resolution of the model is better
#   at a particular wavelength, the test variable prevents the smoothing from 
#   being performed there. By default the variable is active so that this check
#   is performed; it can be deactivated so that the calculations are faster.

def C_resol(lam, z, flux, dflux, sig_spec, sig_M, lam_sig_M=None):
    #lam == array of wavelengths in units of armstrongs.
    #flux == flux (arbitrary units).
    #disp == dispersion in units of pixels.
    #sigM == resolution (standard deviation) of the spectrum (FWHM/2.355) in 
    #        units of armstrongs. This parameter can either be an integer or a 
    #        list or np.narray (in case it is wavelength dependent).
    #lam_sigM == wavelengths along which the model provides a value for a wavelength
    #       dependent standard deviation / resolution.
    
    
    cut = np.asarray([i for i in range(len(lam))])
    
    if type(sig_M)==np.ndarray or type(sig_M)==list:
        if isinstance(lam_sig_M, (np.ndarray, list)) == False or len(sig_M) != len(lam_sig_M):
            print('ERROR: either lam_sigM not provided or dimensions of sigM and lam_sigM do not match')
        else:
            minlam = np.amax([np.amin(lam), np.amin(lam_sig_M*(1.0+z))])
            maxlam = np.amin([np.amax(lam), np.amax(lam_sig_M*(1.0+z))])
            
            cut = np.intersect1d(np.where(lam > minlam)[0], np.where(lam < maxlam)[0])
            
            lam = lam[cut]
            flux = flux[cut]
            sig_spec = sig_spec[cut]
            dflux = dflux[cut]
            
            sig_M = np.interp(lam/(1.0+z), lam_sig_M, sig_M)
            correc_sig = [(-1.0*sig_spec[k]**2.0+sig_M[k]**2.0)**0.5 if sig_M[k]**2.0 > sig_spec[k]**2.0 
                          else 1.0e-6 for k in range(len(sig_spec))]
                
    elif isinstance(sig_M, (float, int)):
        correc_sig = [(-1.0*sig_spec[k]**2.0+sig_M**2.0)**0.5 if sig_M**2.0 > sig_spec[k]**2.0 
                      else 1.0e-6 for k in range(len(sig_spec))]
    
    else:
        print('ERROR: sigM not float/int/array/list') 
    
    correc_flux = np.full_like(flux, np.nan, dtype=float)
    correc_dflux = np.full_like(dflux, np.nan, dtype=float)
    
    mF = np.isfinite(flux)
    mD = np.isfinite(dflux)
    
    for i in range(len(lam)):
        sig = correc_sig[i]
        if (not np.isfinite(sig)) or (sig <= 0):
            continue
    
        K = np.exp(-0.5 * (lam - lam[i])**2 / sig**2)
    
        denF = np.sum(K[mF])
        if denF > 0:
            correc_flux[i] = np.sum(K[mF] * flux[mF]) / denF
    
        denD = np.sum(K[mD])
        if denD > 0:
            correc_dflux[i] = np.sum(K[mD] * dflux[mD]) / denD
    
    return lam, correc_flux, correc_dflux

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#VELOCITY DISPERSION SMOOTHING:
#       Function that produces an smoothed spectra that represents the effect
#   of a broadening due to an inner velocity dispersion of the observed object.
#
#       As an output it provides an array of wavelengths, an array of fluxes 
#   and an array of fluxes corrected by the application of the smoothing.

def SmoothVD(lam, flux, dflux, vDf, vD0 = 0.0, test = True):
    #lam == array of wavelengths in units of armstrongs.
    #flux == flux (arbitrary units).
    #vDf == velocity dispersion at which we aim to smooth the spectrum. It is
    #        expected in units of km s⁻¹
    #vD0 == velocity dispersion at which the spectrum is already smoothed. It is
    #        expected in units of km s⁻¹, if not provided it is assumed that the
    #        provided spectrum is at 0 velocity dispersion.
    #test == boolean variable used to make sure that the velocity dipersion at
    #       which we want to smooth the spectrum is greater than that at which
    #       the provided spectrum is.
    
    correc_flux = np.zeros_like(flux)
    correc_dflux = np.zeros_like(dflux)
    c = 299792.458 #km s-1 !!!!!!!!!!!!!!!
    if vDf > vD0 and test:
        sigmatch = (vDf**2-vD0**2)**(0.5)*lam/c
        for i in range(len(lam)):
            K = np.exp(-1.0*(lam-lam[i])**2.0/(2.0*sigmatch[i]**2.0))
            correc_flux[i] = np.sum(K*flux)/np.sum(K)
            correc_dflux[i] = np.sum(K*dflux)/np.sum(K)
        return lam, correc_flux, correc_dflux
    else:
            print('ERROR: vD0 > vDf')
            print('Initial velocity dispersion = '+str(round(vD0, 0))+
            '. Final velocity dispersion = '+str(round(vDf, 0))+'.')


    
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#VELOCITY DISPERSION CORRECTIONS:


# C_VD is a function aimed to provide a single result (a corrected index plus an
# associated uncertainty) in the form of a (2,) array.
#def C_VD(I_in, vD, arrayVD, arrayC, arraydC, modo = 'AAA'):
#    I_out = np.zeros_like(I_in)
#    elpunto = np.argmin(arrayVD-vD)
#    C = arrayC[elpunto]
#    dC = arraydC[elpunto]
#    if modo == 'AAA':
#        I_out[0] = I_in[0]*C
#        I_out[1] = I_in[1]*C + I_in[0]*dC
#    if modo == 'mag':
#        I_out[0] = I_in[0]+C
#        I_out[1] = I_in[1]+dC
#        
#    return I_out

# C_VD25 is intended to produce simultaneously the correction to an array of dimension
# (NG, 25, 2) where NG is the number of galaxies with 25 measured indices in which the last
# dimension accounts for the value of the index and its uncertainty.
def C_VD25(I_in, resVD, corr_V, corr_C, corr_dC):
    I_out = np.zeros_like(I_in)
    iv = np.zeros(len(resVD))
    for jv in range(len(resVD)):
        iv[jv] = np.argmin(abs(corr_V-resVD[jv]))
    
    indices_25 = [i+37+4 for i in range(25)]
    
    if np.shape(I_in)[1] != len(indices_25):
        print('ERROR: number of indices in I_in not matching functions purpose')
    if np.shape(corr_C)[0] != len(indices_25):
        print('ERROR: number of correction functions not matching functions purpose')
    
    for ii in range(len(indices_25)):
        lasC = [corr_C[ii, int(iv[jjv])] for jjv in range(len(resVD))]
        lasdC = [corr_dC[ii, int(iv[jjv])] for jjv in range(len(resVD))]
        tabindx = open('./pylick/libs/table.dat', 'r')
        for line in tabindx:
            if line[0] != '#':
                l = line.split()
                if l[0] == '0'+str(indices_25[ii]):
                    if l[2] == 'mag':
                        I_out[:, ii, 0] = I_in[:, ii, 0]+lasC
                        I_out[:, ii, 1] = I_in[:, ii, 1]+lasdC
                    else:
                        I_out[:, ii, 0] = I_in[:, ii, 0]*lasC
                        I_out[:, ii, 1] = I_in[:, ii, 1]*lasC+abs(I_in[:, ii, 0])*lasdC
    return I_out



#-------------------------------------------------------------------------------
def C_VD(I_in, resVD, corr_V, corr_C, corr_dC, lista):
    if len(lista) == 25:
        ilista = mclick.I_number(lista)
    elif len(lista) == 30:
        ilista = mclick.I_number30(lista)
    elif len(lista) == 33:
        ilista = mclick.I_numberopen(lista)
    else:
        print('ERROR EN EL NUMERO DE INDICES, consultar corrections.py')
    
    if len(lista) != np.shape(corr_C)[0]:
        print('MAJOR!!! ERROR NUMERO DE INDICES NO MATCHEA NUMERO DE CORRECCIONES!')
    
    itypes = mclick.I_typeopen(lista)
    I_out = np.zeros_like(I_in)
    iv = np.zeros(len(resVD))
    for jv in range(len(resVD)):
        iv[jv] = np.argmin(abs(corr_V-resVD[jv]))
    
    for ii in range(len(lista)):
        lasC = [corr_C[ilista[ii], int(iv[jjv])] for jjv in range(len(resVD))]
        lasdC = [corr_dC[ilista[ii], int(iv[jjv])] for jjv in range(len(resVD))]
        
        if itypes[ii] == "A":
            I_out[:, ii, 0] = I_in[:, ii, 0]*lasC
            I_out[:, ii, 1] = abs(I_in[:, ii, 1]*lasC)+abs(I_in[:, ii, 0]*lasdC)
        else:
            I_out[:, ii, 0] = I_in[:, ii, 0]+lasC
            I_out[:, ii, 1] = I_in[:, ii, 1]+lasdC
    
    return I_out

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#APERTURE CORRECTION:
#       Function that uses the aperture relation in CLEMENS+06 (equation 2)
#   to modify the incides obtained after other corrections have been applied,
#   it accounts for the effect of the aperture of the spectrograph
#   when taking the spectrum.



def C_apert(Iobs, codeI, vD, re):
    #Iobs == Value of the measured index in the galaxy.
    #codeI == Number of the index according to the numbering in table.dat
    #vD == Value of the velocity dispersion of the galaxy.
    #re == Value of the effective radius for that galaxy (see documentation
    #       from SDSS -i.e. Schema Browser- in which this is better explained).

    xyaqui = []
    xytab = open('./libs/xy.dat', 'r')
    for line in xytab:
        if line[0] != '#':
            l = line.split()
            if str(codeI) == l[0]:
                if l[2] == 'YES':
                    Ire10 = Iobs/(1.0+(float(l[3])*np.log10(vD/100.0)+
                    float(l[4]))*(np.log10(15.0/re)/np.log10(5.0)))
                else:
                    Ire10 = Iobs
                    # print('Index with code '+str(codeI)+' has been left '
                    #      + 'as it was: No correction found.')
    
    xytab.close()
    return Ire10

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#FULL CORRECTION:
#       Function that provides a full correction of the spectrum by performing a
#   smoothing of amplitude \sigma_{match}, in which we take the advantage of the
#   natural velocity dispersion as way to get closer to the (downgraded) resolution
#   of the model with respect to which we confronting our indices.

def C_full(lam, z, flux, dflux, sig_spec, sig_IR, vDf, vD0 = 0.0):
    #lam == array of wavelengths in units of armstrongs.
    #z == redshift of the source
    #flux == flux (arbitrary units).
    #dflux == error associated to the flux.
    #sig_spec == dispersion of the spectrum (FWHM/2.355) in units of armstrongs.
    #        since we are using SDSS is to be expected that it is a wavelength dependent
    #        quantity
    #sig_IR == resolution (standard deviation) of the spectrum (FWHM/2.355) in 
    #        units of armstrongs. This parameter can either be an integer or a 
    #        list or np.narray (in case it is wavelength dependent).
    #vD0 == velocity dispersion of the galaxy as measured by SDSS (with the
    #       unit change v = \sigma *c/\lambda, v is equivalent to an amplitude. This
    #       is \sigma = FWHM/2.355 in units of wavelength).
    
    
    
    c = 299792.458 #km s-1 !!!!!!!!!!!!!!!
    cut = np.asarray([i for i in range(len(lam))])
    if vDf >= vD0:
        sigVD2 = lam**2.0 * (vDf**2.0 - vD0**2.0)/c**2.0
    else:
        sigVD2 = np.zeros_like(lam)
        print('CHECK!!! vDf > vD0 ')
    
    if type(sig_spec)==np.ndarray or type(sig_spec)==list:
        if len(sig_spec) != len(lam):
            print('Dimension of sigma of incoming object is different from dimension of its spectrum')
            return lam, flux, dflux
    else:
        sig_spec = np.zeros_like(lam)+sig_spec
    
    if type(sig_IR)==np.ndarray or type(sig_IR)==list:
        if len(sig_IR) != len(lam):
            print('Mi ninio esto no es jauja, pones sig_IR en el mismo tamagno que el espectro entrante')
            return lam, flux, dflux
    else:
        sig_IR = np.zeros_like(lam)+sig_IR
    
    
    sigmatch = [( -1.0*sig_spec[k]**2.0 + sig_IR[k]**2.0 + sigVD2[k])**0.5
                      if -1.0*sig_spec[k]**2.0 + sig_IR[k]**2.0 + sigVD2[k] > 0
                      else 1.0e-6 for k in range(len(lam))]
    
    correc_flux = np.zeros_like(flux)
    correc_dflux = np.zeros_like(flux)
    for i in range(len(lam)):
        K = np.exp(-1.0*(lam-lam[i])**2.0/(2.0*sigmatch[i]**2.0))
        correc_flux[i] = np.sum(K*flux)/np.sum(K)
        correc_dflux[i] = np.sum(K*dflux)/np.sum(K)
    return lam, correc_flux, correc_dflux

#------------------------------------------------------------------------------










'''
#------------------------------------------------------------------------------

#DEPRECATED



#RESOLUTION CORRECTION:
#       Function that takes in an array or list of wavelengths,
#   an array of fluxes and an array of the dispersion of the spectrograph (the 
#   latter in units of pixels).
#
#       As an output it provides an array of wavelengths in the resframe (or
#   provided redshift), an array of fluxes and an array of fluxes corrected by
#   the dispersion of the model; in case the resolution of the model is better
#   at a particular wavelength, the test variable prevents the smoothing from 
#   being performed there. By default the variable is active so that this check
#   is performed; it can be deactivated so that the calculations are faster.

def C_resol(lam, flux, disp, sigM, z = 0.0, lam_sigM=None):
    #lam == array of wavelengths in units of armstrongs.
    #flux == flux (arbitrary units).
    #disp == dispersion in units of pixels.
    #sigM == resolution (standard deviation) of the spectrum (FWHM/2.355) in 
    #        units of armstrongs. This parameter can either be an integer or a 
    #        list or np.narray (in case it is wavelength dependent).
    #z == redshift of the galaxy, 0 by default.
    #lam_sigM == wavelengths along which the model provides a value for a wavelength
    #       dependent standard deviation / resolution.
    
    lambdasRF = lam/(1.0+z)
    dlams = np.append(np.diff(lambdasRF), np.diff(lambdasRF)[-1])
    sigmaspec = disp*dlams
    correc_flux = np.zeros_like(flux)
    if type(sigM)==np.ndarray or type(sigM)==list:
        if isinstance(lam_sigM, (np.ndarray, list)) == False or len(sigM) != len(lam_sigM):
            print('ERROR: either lam_sigM not provided or dimensions of sigM and lam_sigM do not match')
        else:
            newsigM = np.interp(lambdasRF, lam_sigM, sigM)
            correc_sig = [(-1.0*sigmaspec[k]**2.0+newsigM[k]**2.0)**0.5 if newsigM[k]**2.0 > sigmaspec[k]**2.0 
                          else 1.0e-6 for k in range(len(sigmaspec))]
                
    if isinstance(sigM, (float, int)):
        correc_sig = [(-1.0*sigmaspec[k]**2.0+sigM**2.0)**0.5 if sigM**2.0 > sigmaspec[k]**2.0 
                      else 1.0e-6 for k in range(len(sigmaspec))]
    else:
        print('ERROR: sigM badly defined')
        test = False
    
    
    for i in range(len(lambdasRF)):
        K = np.exp(-1.0*(lambdasRF-lambdasRF[i])**2.0/(2.0*correc_sig[i]**2.0))
        correc_flux[i] = np.sum(K*flux)/np.sum(K)
    return lambdasRF, correc_flux


#Function to correct by the resolution an interpolated spectrum.
def Corr_R(waveRF, flux, wdisp, sigM, lam_sigM=None):
    #The array of wavelengths is expected to be writen in the same units of the dispersion array of the model
    #The array of dispersions of the observed spectrum is expected in units of pixels
    
    dlams = np.append(np.diff(waveRF), np.diff(waveRF)[-1])
    sigmaspec = wdisp*dlams
    correc_flux = np.zeros_like(flux)
    if type(sigM)==np.ndarray:
        if isinstance(lam_sigM, np.ndarray) == False or len(sigM) != len(lam_sigM):
            print("ERROR: either lam_sigM not provided or dimensions of sigM and lam_sigM do not match")
        else:
            newsigM = np.interp(waveRF, lam_sigM, sigM)
            correc_sig = (-1.0*sigmaspec**2.0+newsigM**2.0)**0.5
                
    if (type(sigM) == float) or (type(sigM) == int):
        correc_sig = (-1.0*sigmaspec**2.0+sigM**2.0)**0.5
    else:
        print('ERROR: sigM badly defined')
    
    if (type(sigM) == np.ndarray) or (type(sigM) == int) or (type(sigM) == float):
        for i in range(len(waveRF)):
            if np.isscalar(sigM):
                if sigM**2.0<sigmaspec[i]**2.0:
                    print('imaginary sqrt corrected')
                    correc_sig[i] = sigmaspec[i]
            else:
                if newsigM[i]**2.0<sigmaspec[i]**2.0:
                    print('imaginary sqrt corrected')
                    correc_sig[i] = sigmaspec[i]
            K = np.exp(-1.0*(flux-flux[i])**2.0/(2.0*correc_sig[i]**2.0))
            correc_flux[i] = np.sum(K*flux)/np.sum(K)
        return waveRF, flux, correc_flux


#Function to correct by the resolution the spectrum that can be obtained by the
#ID provided as first argument (expected the file path to exist).
def Corr_R_ID(ID, sigM, z = 0.0, lam_sigM=None):
    dat = fits.open(ID)
    lambdas = 10.0**dat[1].data["loglam"]
    sigmapix = dat[1].data["wdisp"]
    lambdasRF = lambdas/(1.0+z)
    dlams = np.append(np.diff(lambdasRF), np.diff(lambdasRF)[-1])
    sigmaspec = sigmapix*dlams
    flux = dat[1].data["flux"]
    correc_flux = np.zeros_like(flux)
    if type(sigM)==np.ndarray:
        if isinstance(lam_sigM, np.ndarray) == False or len(sigM) != len(lam_sigM):
            print("ERROR: either lam_sigM not provided or dimensions of sigM and lam_sigM do not match")
        else:
            newsigM = np.interp(lambdasRF, lam_sigM, sigM)
            correc_sig = (-1.0*sigmaspec**2.0+newsigM**2.0)**0.5
                
    if (type(sigM) == float) or (type(sigM) == int):
        correc_sig = (-1.0*sigmaspec**2.0+sigM**2.0)**0.5
    else:
        print('ERROR: sigM badly defined')
    
    if (type(sigM) == np.ndarray) or (type(sigM) == int) or (type(sigM) == float):
        for i in range(len(lambdasRF)):
            if np.isscalar(sigM):
                if sigM**2.0<sigmaspec[i]**2.0:
                    print('imaginary sqrt corrected')
                    correc_sig[i] = sigmaspec[i]
            else:
                if newsigM[i]**2.0<sigmaspec[i]**2.0:
                    print('imaginary sqrt corrected')
                    correc_sig[i] = sigmaspec[i]
            K = np.exp(-1.0*(flux-flux[i])**2.0/(2.0*correc_sig[i]**2.0))
            correc_flux[i] = np.sum(K*flux)/np.sum(K)
        return lambdasRF, flux, correc_flux
'''
