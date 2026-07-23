# Vendored from the author's personal library (ARTICLE3/aux13_burnin_maker.py),
# trimmed to drop a dead `from getdist import loadMCSamples, plots` (imported,
# never used) and move the module-level `import matplotlib.pyplot` inline to
# the one `seplotea` branch that actually uses it -- same trims as
# aux/burnin.py (the TMJ-stage burn-in maker this mirrors, which fixes
# ndim=3; this one takes ndim as an argument since the cosmographic fit's
# parameter count depends on how many vd groups are included).

from numpy import *
import numpy as np
import h5py
import emcee

#Define the autocorrelation function for a 1D chain
def next_pow_two(n):
    #Given an integer "n", return the smallest powert of two larger than or equal to "n"
    i=1
    while i<n:
        i=i<<1
    return i

def autocorr_func_1d(x,norm=True):
    #Compute the autocorrelation function for a given 1D chain
    x=atleast_1d(x)
    n=next_pow_two(len(x))

    #Compute the Fast Fourier Transform and, from there, the autocorrelation function
    f=fft.fft(x-mean(x),n=2*n)
    acf=fft.ifft(f*conjugate(f))[:len(x)].real
    acf/= 4*n

    #Normalize it (optional)
    if norm:
        acf/=acf[0]
    return acf

def auto_window(taus,c):
    m=arange(len(taus))<c*taus #Create an array of "True" and "False" values to find what is the smallest integer larger than c*tau
    if any(m):
        return argmin(m)
    return len(taus)-1

#Foreman-Mackey
def autocorr_fm(y,c=5): #"y" is a 2D array of shape (nwalkers,Niter)
    f=zeros(y.shape[1])
    for yy in y: #For each walker, compute the autocorrelation and then average them
        f+=autocorr_func_1d(yy)
    f/=len(y)
    taus=2*cumsum(f)-1 #Subtract 1 because we summed over over 0 twice!
    window=auto_window(taus,c)
    return taus[window]

#Select whhat to do:
#1. "read_h5" to read the h5 output of the MCMC and generate the .txt with the loglikelihood for getdist
#2. "generate_burnin" to generate the burned-in .txt file along with the .paramnames and .ranges files
#3. "posterior_plots" to make corner plots or 1D posteriors

def burnin(namefile, name_file_burnin = False, txt = 'no', what1="read_h5", what2="generate_burnin", seplotea=False, extrafiles = False,
           ndim = 11, nwalkers = 40):

    if what1=="read_h5":

        reader=emcee.backends.HDFBackend('%s.h5'%namefile,read_only=True)
        chain=reader.get_chain()

        flatchain=reader.get_chain(flat=True)
        logchain=reader.get_log_prob(flat=True)

        ones_2D=ones((len(logchain),1))
        logchain_2D=reshape(logchain,(len(logchain),1))

        getdist_chain=hstack((ones_2D,logchain_2D))
        getdist_chain=hstack((getdist_chain,flatchain))



        if txt == 'yes':
            savetxt('%s.txt'%namefile,getdist_chain)

        autocorr_param=zeros(ndim)
        autocorr_iter=zeros(len(chain))

        print(len(chain))


        for i in arange(len(chain)): #For every iteration
            if i % 100:
                continue
            cum_chain=chain[:i+1,:,:]
            #print(shape(cum_chain))
            for j in arange(ndim): #For every parameter
                cum_chain_param=cum_chain[:,:,j]
                autocorr_param[j]=autocorr_fm(cum_chain_param.T)
            #autocorr_iter[i]=autocorposterior_plotr_param[9]
            autocorr_iter[i]=mean(autocorr_param)
        burn_in=2*int(max(autocorr_param))
            #otprint(autocorr_iter[i])
            #plt.plot(i+1,autocorr_iter[i],'ro')
            #plt.pause(0.05)

        if seplotea:
            import matplotlib.pyplot as plt
            plt.figure(1)
            plt.ylabel('Mean autocorrelation time')
            plt.xlabel('Number of iterations')
            plt.plot(arange(0,len(chain),1),autocorr_iter,'ro')
            #plt.savefig('autocorrelation_time_20_BS_cross_beta_fixed_30')
            plt.show()


    if what2=="generate_burnin":

        #Select the chain, the number of walkers and the number of burn-in iterations
        name_file = '%s'%namefile  # Name of the .txt file without the '.txt'
        burnin_it = burn_in
        print(burnin_it)


        sample=getdist_chain #Load the unburnt .txt file

        n_param=len(sample[0,:])-2 #Number of parameters of the sample

        burnin_rows=burnin_it*nwalkers #Number of rows to be burnt it
        sample_burnin=sample[burnin_rows:,:]
        if not name_file_burnin:
            name_file_burnin=name_file+'_burnin'
        savetxt(name_file_burnin+'.txt',sample_burnin)
