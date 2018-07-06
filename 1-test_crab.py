from __future__ import print_function, division

import os
import numpy
from scipy.ndimage.filters import gaussian_filter
import spipy

CWD = os.getcwd()



if __name__ == '__main__':
	#	import the data
	data = spipy.get_data(CWD+'/Inputs/evts_det_spec.fits', 50)
	energy_boundaries = spipy.get_energybins(CWD+'/Inputs/energy_boundaries.fits')

	#	plot the spectrum for pointing 3, detector 4
	spipy.plotspectrum_raw(3, 4, energy_boundaries, data, 'rawspec', 'Outputs')

	
	# e.g find the mean number of counts per detector observed during pointing 3
	pt3_cts = numpy.mean(data[:,3,:], axis = 1)
	d_matrix = numpy.array(pt3_cts)

	# #	plot the map of which detectors saw which photons during pointing 3
	# spipy.plot_map(d_matrix, 3, w=500, dpi=72., title='SPI Hit Map')

	# e.g find the mean number of counts per detector observed the first 50 pointings
	# for i in range(0, 50):
	# 	pt3_cts = numpy.mean(data[:,i,:], axis = 1)
	# 	d_matrix = numpy.array(pt3_cts)

	# 	#	plot the map of which detectors saw which photons
	# 	spipy.plot_map(d_matrix, i, w=500, dpi=72., title='SPI Hit Map')

	# e.g. plot the number of counts in a specific energy channel per detector (this loop generates plots for all 50 channels)
	# for i in range(0,50):
	# 	pt3_cts = data[:,3,i]
	# 	d_matrix = numpy.array(pt3_cts)

		#	plot the map of which detectors saw which photons
	spipy.plot_map(d_matrix, 3, 'Outputs', w=500, dpi=72., title='SPI Hit Map')

	#	import the data that has been fitted to remove the background
	fitted = spipy.get_fitspec(CWD+'/fitted/spectra_Crab.fits')

	#	fit the data with a powerlaw
	fitparams = spipy.fit_noRMF(CWD+'/fitted/spectra_Crab.fits', energy_boundaries)

	#	plot the data and the linear least squares fitted to the data 
	spipy.plot_spectrumfit_noRMF(fitted, energy_boundaries, fitparams, 'Outputs', 'CRAB_specnoRMF', ext = 'png')

	#	import the RMF matrix
	RMF = spipy.get_RMFmatrix(CWD+'/fitted/spectral_response.rmf.fits', log = True)
	
	#	Gaussian smoothing filter to make the plot look a bit better
	sigma = 0.4
	RMF = gaussian_filter(RMF, sigma)


	x = energy_boundaries
	#	get the new energy boundaries for the RMF matrix
	y = spipy.get_RMF_energyboundaries(CWD+'/fitted/spectral_response.rmf.fits')

	#	plot the visualization of the RMF matrix
	spipy.plot_RMFmatrix(RMF, x, y, 'Outputs', 'CRAB_RMFvis')
