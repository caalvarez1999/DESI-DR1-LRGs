"""
.. authors::   Michele Moresco <michele.moresco_at_unibo.it>
..             Nicola Borghi   <nicola.borghi6_at_unibo.it>
..             Salvatore Quai  <squai@uvic.ca>
*General purpose*:
***
*Imports*::
	import numpy as np
	import astropy.io.fits as pyfits
	import os
	from firefly_dust import get_dust_radec	
	import astropy.cosmology as cc
	import astropy.units as uu
	import astropy.constants as const
"""

import sys
import warnings
# warnings.filterwarnings("ignore")
from time import perf_counter as clock

import numpy as np
from astropy.io import ascii
from astropy.table import Table, MaskedColumn

from pylick.plot import plot_spectrum
from pylick.indices import IndexLibrary
from pylick.measurements import MeasSpectrum
from pylick._config import __table__, __dirRes__



class Catalog:
	"""

	Function for resampling spectra (and optionally associated
	uncertainties) onto a new wavelength basis.

	Parameters
	----------
	ID : numpy.ndarray
		Galaxy unique identifier(s)

	spec_wave : numpy.ndarray
		1D array containing the current wavelength sampling of the
		spectrum or spectra.

	spec_fluxes : numpy.ndarray
		Array containing spectral fluxes at the wavelengths specified in
		spec_wavs, last dimension must correspond to the shape of
		spec_wavs. Extra dimensions before this may be used to include
		multiple spectra.

	spec_errs : numpy.ndarray (optional)
		Array of the same shape as spec_fluxes containing uncertainties
		associated with each spectral flux value.

	settings : dict

	verbose : bool (optional)
		Setting verbose to False will suppress the default warning about
		new_wavs extending outside spec_wavs and "fill" being used.


	Output
	------

	The output values are stored as attributes of the *** class.

	pylick.val : np.ndarray
		measured indices
	
	pylick.err : np.ndarray 
		measured errors

	"""



	def _print_hello(self):
		self._print_console("################################################################", 'r')
		self._print_console("#### PyLick: A Python program to measure spectral indices.  ####", 'r')
		self._print_console("####                                                        ####", 'r')
		self._print_console("####        authors: M. Moresco, N. Borghi, S. Quai         ####", 'r')
		self._print_console("################################################################\n", 'r')

	def _print_console(self, text, color):
		# if color == 'bad':
		# 	print("\x1b[0;37;41m"+text+"\x1b[0m")
		# if color == 'good':
		# 	return '\x1b[0;37;42m'+text+'\x1b[0m'
		# if color == 'status':
		# 	return '\x1b[1;34;40m'+text+'\x1b[0m'
		if color == 'w':
			print("\x1b[0;49;97m"+text+"\x1b[0m")
		if color == 'c':
			print("\x1b[36m"+text+"\x1b[0m")
		if color == 'r':
			print("\x1b[91m"+text+"\x1b[0m")
		if color == 'y':
			print("\x1b[1;33m"+text+"\x1b[0m")

	def _print_table(self, names, ind, err, BPR):
		res_table = Table(masked=True)
		res_table["Index"] = names
		res_table["Value"] = MaskedColumn(data=ind, mask=np.isnan(ind))
		res_table["Error"] = MaskedColumn(data=err, mask=np.logical_or(np.isnan(err), err==0) )
		res_table["BPR"] = BPR

		# Remove if index is not measured
		res_table.remove_rows(np.where(np.isnan(ind))[-1])
		res_table.pprint_all()

	def _print_ProgressBar(self, iteration, total, prefix = 'Progress:', suffix = '', length = 40, fill = '█'):
		prefix =  ("{:8s}").format(prefix)
		percent = ("{:d}").format(int(100 * (iteration / float(total))))
		filledLength = int(length * iteration // total)
		bar = fill * filledLength + '-' * (length - filledLength)
		# print('\rID {:s}: {:s}\t'.format(prefix, bar), end = "")
		print(f'  ID {prefix}: |{bar}| {percent}% {suffix}', end = "\r")
		# Print New Line on Complete
		if iteration == total: 
			print()		


	def _save_table(self, data_names, data_res, res_prepend, res_prename, res_format, dir_out):

		out_data_nmes = [res_prename+x for x in data_names]

		if res_prepend is None:
			res_table = Table(masked=True)
			res_table["ID"] = self.IDs
		else:
			res_table = Table(res_prepend, masked=True)

		for i, colname in enumerate(out_data_nmes):

			res_table.add_column(data_res[:,i], name=colname)
		
		res_table.write(dir_out+'/out.'+res_format, overwrite=True)


	def __init__(self, IDs, load_spectrum, index_list, meas_method='int', 
				z=None, index_table=None, bpr_thres=None, 
				res_prepend=None, res_prename='', res_format='.ecsv', res_dir=None,
				plot=False, plot_settings={},
				verbose=False):

		# Do extensive checking of possible input errors
		self.IDs			= np.array(IDs, dtype=str)
		self.index_list		= index_list
		self.meas_method	= meas_method
		self.index_table	= index_table 
		self.verbose		= verbose
		self.bpr_thres      = bpr_thres

		self.Nobj			= len(IDs)
		self.Nidx			= len(index_list)
		self.isrest			= True

		if index_table is None:
			index_table = __table__

		if res_dir is None:
			dir_out = '.'

		if self.bpr_thres is None:
			self.bpr_thres = 0.15

		data_res = np.empty((self.Nobj, self.Nidx*3), dtype=float)



		# GO
		t = clock()
		self._print_hello()

		lib 				= IndexLibrary(file_table=index_table, key_list=index_list)
		data_names 			= ['']*self.Nidx*3
		data_names[::3]  	= lib.names
		data_names[1::3] 	= [name+'_err' for name in lib.names]
		data_names[2::3] 	= [name+'_cont' for name in lib.names]

		if not verbose:
			self._print_ProgressBar(0, self.Nobj)

		if hasattr(z, "__len__"):
			self.isrest	= False

	

		for i, ID in enumerate(self.IDs):

			self.spectrum 	= np.array(load_spectrum(ID))
			self.spec_wave 	= self.spectrum[0] if self.isrest else self.spectrum[0]/(1+z[i])
			self.spec_flux 	= self.spectrum[1]
			self.spec_ferr 	= self.spectrum[2]
			self.spec_mask  = self.spectrum[3].astype(bool)

			#self.spec_mask = (self.spec_mask) & (self.spec_flux>0) & (self.spec_ferr>0)
			self.spec_mask = (self.spec_mask) & (self.spec_ferr>0)

			self.spec_ferr[~self.spec_mask] = 9.9*10**99.


			assert self.spec_wave.size == self.spec_flux.size == self.spec_ferr.size, \
					'PyLick: spec_wave, spec_flux, (spec_err) must have the same size'
			
			meas    = MeasSpectrum(self.spec_wave, self.spec_flux, self.spec_ferr, self.spec_mask,
								lib.regions, lib.units, lib.names, self.meas_method,
								BPR_thres=self.bpr_thres, nans=np.nan, verbose=self.verbose)
			
			# Fill results array ind1, ind1_err, ind2, ind2_err, ...
			data_res[i][::3]  = meas.ind
			data_res[i][1::3] = meas.err
			data_res[i][2::3] = meas.ps_cont

			# Plot
			if plot is True:
				plot_spectrum(ID, self.spec_wave, self.spec_flux, self.spec_ferr, self.spec_mask,
							lib.regions, lib.tex_names, lib.units, meas.finite,
							plot_settings)

			# Print status
			if verbose:
				self._print_console("\n < Galaxy ID {:s} >".format(ID),'y')
				self._print_table(lib.names, meas.ind, meas.err, meas.BPR) 	
			else: 
				self._print_ProgressBar(i+1, self.Nobj, ID)

		# Save results
		self._save_table(data_names, data_res, res_prepend, res_prename, res_format, dir_out)


		if True: 		
			self._print_console("\nElapsed time: {:.2f} s\n".format(clock() - t), 'r')


		self.data_names = data_res
		self.data_res   = data_res

		print("") # End






class Galaxy:

	def __init__(self, ID, index_list, spec_wave, spec_flux, spec_err=None, spec_mask=None,
				 z=None, meas_method='int', bpr_thres=None, index_table=None, plot=False, plot_settings={}):

		# Do extensive checking of possible input errors
		self.IDs			= ID
		self.spec_wave		= np.array(spec_wave)
		self.spec_flux		= np.array(spec_flux)
		self.index_list		= index_list
		self.meas_method	= meas_method
		self.index_table	= index_table 
		self.z				= z
		self.bpr_thres      = bpr_thres
		self.Nidx			= len(index_list)

		if self.bpr_thres is None:
			self.bpr_thres = 0.15

		if self.index_table is None:
			self.index_table = __table__


		self.spec_wave = self.spec_wave if self.z is None else self.spec_wave/(1+self.z)
		self.spec_ferr = np.array(spec_err) if spec_err is not None else np.zeros_like(self.spec_wave)
		self.spec_mask = np.array(spec_mask, dtype=bool) if spec_mask is not None else np.zeros_like(self.spec_wave, dtype=bool)
		
				
		assert self.spec_wave.size == self.spec_flux.size == self.spec_ferr.size, \
				'PyLick: spec_wave, spec_flux, (spec_err) must have the same size'

		#self.spec_qual = (~self.spec_mask) & (self.spec_flux>0) & (self.spec_ferr>0)
		self.spec_qual = (~self.spec_mask) & (self.spec_ferr>0)
		self.spec_ferr[~self.spec_qual] = 9.9*10**99.


		lib	= IndexLibrary(file_table=self.index_table, key_list=self.index_list)


		meas    = MeasSpectrum(self.spec_wave, self.spec_flux, self.spec_ferr, self.spec_qual,
								lib.regions, lib.units, lib.names, self.meas_method,
								BPR_thres=self.bpr_thres, nans=np.nan, verbose=False)
			
		self.vals = meas.ind 
		self.errs = meas.err
		self.SNR = meas.snr_id
		

		# Plot
		if plot is True:
			plot_spectrum(ID, self.spec_wave, self.spec_flux, self.spec_ferr, self.spec_mask,
						lib.regions, lib.tex_names, lib.units, meas.finite,
						plot_settings)


		#print("") # End
