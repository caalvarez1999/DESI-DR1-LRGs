#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:13:01 2024

@author: carlos
"""

# Vendored from the author's personal library (ARTICLE3/aux13_MCMC.py),
# verbatim other than repointing MYLIBS.TAYLOR_emcee_tz[_widepriors] at
# aux.taylor_tz[_wide] to match this repo's vendoring convention.

from __future__ import division, print_function, absolute_import, unicode_literals

import numpy as np

import emcee as emcee

import aux.taylor_tz_wide as emw
import aux.taylor_tz as emn

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
def MCMC_one(z1, z0, t, dt, grupoVD, priors = '', ORDER = 3, name = '', nwalkers = 40, Niter = 50000, var_nms = ["age160180", "age180200", "age200225", "age225250", "age250280", "age280320", "age320355", "age355400", "Hz0", "qz0", "jz0"], chains_dir = './chains/cosmographic/'):

    if priors == 'wide':
        lib = emw
    else:
        lib = emn


    #DEFINITIONS THAT ARE GOING TO BE USEFUL:
    ndim = len(var_nms)                                                 #NUMBER OF PARAMETERS TO FIT
    print('number of parameters is '+str(ndim))

    carpeta = str(chains_dir)

    #INITIALISE THE MCMC:
    # Initial walkers' position -> shape(pos)=(nwalkers, ndim)
    inipos = [[lib.fpriors(var)[0] + np.random.uniform(0, lib.fpriors(var)[1]-lib.fpriors(var)[0]) for var in var_nms]
              for i in range(nwalkers)]

    from functools import partial
    lnpost=partial(lib.ln_post,var_nms=var_nms,tobs = t, dtobs = dt, z = z1, z0 = z0, grupo = grupoVD, orden = ORDER)

    fullnamefile = carpeta+str(name)                                    #FULL NAME OF THE FILE
    filename_back = fullnamefile+'.h5'                                  #+ .h5
    backend=emcee.backends.HDFBackend(filename_back)                    #Set up the backend
    backend.reset(nwalkers,ndim)
    sampler=emcee.EnsembleSampler(nwalkers, ndim, lnpost, backend=backend, moves=[(emcee.moves.DEMove(),0.8),(emcee.moves.DESnookerMove(),0.2),])

    index=0
    autocorr=np.empty(Niter)
    old_tau=np.inf



    #-------------------------------------------------------------------
    #PERFORM THE MCMC:
    for sample in sampler.sample(inipos, iterations=Niter, progress=True):
        if sampler.iteration%10:
            continue                                                    #CHECK CONVERGENCE EVERY 50 ITERATIONS

        tau=sampler.get_autocorr_time(tol=0)                            #COMPUTE THE MEAN AUTOCORRELATION TIME SO FAR
        autocorr[index]=np.mean(tau)
        index+=1
        if all(tau*50<sampler.iteration) and all(abs(old_tau-tau)/tau<0.01):
            break                                                       #IF THE CHAIN CONVERGED: BREAK
        old_tau=tau
    #-------------------------------------------------------------------



    #COMPUTE THE PROBABILITY THAT A CHI2 WITH dof DEGREES OF FREEDOM IS SMALLER THAN THE OBSERVED CHI2 VALUE
    print('done')
