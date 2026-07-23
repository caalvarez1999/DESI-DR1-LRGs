#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 12:13:01 2024

@author: carlos
"""

# Vendored from the author's personal library (STACKING/MCMC2_MCMC.py),
# trimmed to just MCMC_one (the full t/Z/afe fit -- what
# 9_FULLTMJ.py/pipeline/20_SPSfitTMJ.py use) and its chisquare_index
# dependency. Left out: MCMC_metal/chisquare_index_metal (a Z/afe-only fit
# variant, not used by the TMJ pipeline stage) and two module-level
# `np.load(...)` calls the original file made unconditionally
# (T10_IDS_finer01001001.npy / T10_MILES_finer01001001.npy, ~700 MB each) --
# neither is referenced anywhere inside MCMC_one, they belong to an
# unrelated fit variant in the same file.

from __future__ import division, print_function, absolute_import, unicode_literals

import numpy as np

import emcee as emcee
from aux.emcee_lick import *
from scipy.special import gammainc


#-----------------------------------------------------------------------
#DEFINITION OF THE FUNCTION FOR THE CHI2 COMPUTATION OF EACH SINGLE INDEX
#-----------------------------------------------------------------------
def chisquare_index(index_number, vals, index_value, index_error, modelo, t_array, Z_array, afe_array):
    t=vals[0]
    Z=vals[1]
    alpha=vals[2]
    index_model=I_model(t,Z,alpha,index_number, modelo, t_array, Z_array, afe_array)
    return (index_value-index_model)**2/index_error**2
#-----------------------------------------------------------------------



#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
def MCMC_one(namefile, medidas, elmodelo, T_ARRAY, Z_ARRAY, AFE_ARRAY, lista, nwalkers = 40, Niter = 50000, var_nms = ["t","Z","alpha"], span_pylick = 30):

    #DEFINITIONS THAT ARE GOING TO BE USEFUL:
    ndim = len(var_nms)                                                 #NUMBER OF PARAMETERS TO FIT
                                                                        #NUMBERS ASSOCIATED TO THE INDICES USED IN THE FIT
    if span_pylick == 25:
        indices = I_number(lista)
    if span_pylick == 27:
        indices = I_number27(lista)
    if span_pylick == 29:
        indices = I_number29(lista)
    if span_pylick == 30:
        indices = I_number30(lista)

    measured_indices = np.zeros((span_pylick, 3))
    for ii in range(np.shape(measured_indices)[0]):
        measured_indices[ii, :] = [medidas[ii, 0], medidas[ii, 1], ii]



    #CHECK THAT THE medidas ARRAY HAS AS SHAPE (span_pylick, 2):
    if np.shape(medidas)[0] != span_pylick:
        print('ERROR: input measurements do not contain span_pylick indices')



    #INITIALISE THE MCMC:
    # Initial walkers' position -> shape(pos)=(nwalkers, ndim)
    inipos = [[fpriors(var, T_ARRAY, Z_ARRAY, AFE_ARRAY)[0] + np.random.uniform(0, fpriors(var, T_ARRAY, Z_ARRAY, AFE_ARRAY)[1]-fpriors(var, T_ARRAY, Z_ARRAY, AFE_ARRAY)[0]) for var in var_nms]
              for i in range(nwalkers)]

    from functools import partial
    lnpost=partial(ln_post,var_nms=var_nms,data=get_data(measured_indices, indices),
              modelo = elmodelo, t_array = T_ARRAY, Z_array = Z_ARRAY, afe_array = AFE_ARRAY)


    filename_back = namefile+'.h5'                                      #+ .h5
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
    reader=emcee.backends.HDFBackend(namefile+'.h5',read_only=True)
    chain=reader.get_chain()
    flatchain=reader.get_chain(flat=True)
    logchain=reader.get_log_prob(flat=True)
    ones_2D=np.ones((len(logchain),1))
    logchain_2D=np.reshape(logchain,(len(logchain),1))
    getdist_chain=np.hstack((ones_2D,logchain_2D))
    getdist_chain=np.hstack((getdist_chain,flatchain))
    sample=getdist_chain

    max_logL=max(sample[:,1])
    w=np.where(sample[:,1]==max_logL)[0]

    if len(w)>0:
        w=w[0]

    vals_bf=np.array([sample[w,2],sample[w,3],sample[w,4]]) #t,Z,alpha

    chisquares={}

    for index in lista:
        if span_pylick == 25:
            index_number=I_number([index])
        if span_pylick == 27:
            index_number=I_number27([index])
        if span_pylick == 29:
            index_number=I_number29([index])
        if span_pylick == 30:
            index_number=I_number30([index])

        index_value=measured_indices[index_number,0]
        index_error=measured_indices[index_number,1]

        chisquares['%s'%index]=chisquare_index(index_number,vals_bf,index_value,index_error,
                                               modelo = elmodelo, t_array = T_ARRAY, Z_array = Z_ARRAY, afe_array = AFE_ARRAY)

    chisquare_total=np.sum(list(chisquares.values()))

    dof=len(lista)-3 #Degrees of freedom of chi square

    P=gammainc(dof/2,chisquare_total/2)

    return P
