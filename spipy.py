from __future__ import print_function, division

import os


from matplotlib import animation
from matplotlib import colors, cm
from matplotlib.collections import RegularPolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

from scipy import optimize

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

import warnings
warnings.filterwarnings("ignore")


def get_data(filename, channels):
	print('Reading in raw data file '+filename+'...')
	detector_counts = Table.read(filename, format = 'fits')
	detectors = 19
	pointings = int(len(detector_counts)/19)
	channels = channels
	#	TODO: automatically read the number of detectors, pointings, channels
	DATA= np.zeros((detectors,pointings,channels))
	print('Building data structure with dimensions:')
	print('Detectors: %.3g' % (detectors))
	print('Pointings: %.3g' % (pointings))
	print('Channels: %.3g' % (channels))

	for i in range(0,detectors):
		for j in range(0,pointings):
			counts = detector_counts[i*j][0]
			for k in range(0,channels):
				DATA[i][j][k] = counts[k]
	return DATA

def get_energybins(filename):
	energy_boundaries = Table.read(filename, format = 'fits')
	#	get the boundaries of the energy bins
	return [i for i in energy_boundaries[1][:]]

def get_fitspec(filename):
	fitted = Table.read(filename, format = 'fits', hdu=2)
	channel = [i for i in fitted[0][:]]
	flux = [i for i in fitted[1][:]]
	flux_err = [i for i in fitted[2][:]]
	dflux = [i for i in fitted[6][:]]
	dflux_err = [i for i in fitted[7][:]]
	return np.column_stack([channel, flux, flux_err, dflux, dflux_err])


def plotspectrum_raw(pointing, detector, energy_boundaries, data, outname, outdir, scalex = 'log'):
	plt.figure(figsize = (7,7))
	plt.hist(energy_boundaries, weights = data[detector][pointing][:], bins = energy_boundaries, color ='k', histtype = 'step')
	plt.xscale(scalex)
	plt.xlim([27, 1400])
	plt.xlabel('E/keV')
	plt.ylabel('total counts')
	plt.savefig(outir+'/'+outname+'.png')
	return print('file saved in current directory')


def get_RMFmatrix(filename, log = True):
	RMF_file = Table.read(filename, format = 'fits', hdu = 2)

	RMFmat_raw = RMF_file[5][:]

	original_bins = max(RMF_file[4][:])
	new_bins = len(RMF_file[0][:])

	RMF = []
	for i in range(0, new_bins):
		row = (RMFmat_raw[i].tolist())
		while len(row)<original_bins:
			row.append(0)
		if log == True:
			row = np.log10(row)
			row[row == -np.inf] = 0
			RMF.append(row)
		elif log == False:
			RMF.append(row)

	RMF = np.asarray(RMF)

	return RMF

def get_RMF_energyboundaries(filename):
	RMF_file = Table.read(filename, format = 'fits', hdu = 2)
	boundaries = RMF_file[0][:] 
	return boundaries

def fit_noRMF(filename, energy_boundaries):

	fitted = get_fitspec(filename)

	#   overly simplistic fitting
	powerlaw = lambda x, amp, index: amp * (x**index)

	#make data linear
	elin = np.log10(energy_boundaries)
	flux_lin = np.log10(fitted[:,3])
	fluxerrerr_lin = np.log10(fitted[:,4])


	fitfunc = lambda p, x: p[0] + p[1] * x
	errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

	pinit = [-2, 0]
	out = optimize.leastsq(errfunc, pinit, args=(elin, flux_lin, fluxerrerr_lin), full_output=1)
	pfinal = out[0]
	covar = out[1]

	index = pfinal[1]
	amp = 10.0**pfinal[0]

	print('Power Law Index: %.3g, Amplitude: %.3g' % (index, amp))

	return [index, amp]

def plot_map(d_matrix, pointing, outdir,
         w=1080,
        dpi=72.,
        title='SPI hit map'):
	"""
	Plot hexagon map where each detector is represented by a hexagon. The hexagon
	color is given by the number of counts incident on the detector (D-Matrix)

	Args:,
	- d_matrix: array contaning the distances between each neuron
	- w: width of the map in inches
	- title: map title

	Returns the Matplotlib SubAxis instance
	"""
	#	Define the SPI detector grid
	grid = {'centers': np.array([[0,0],
	[ 1, 0],
	[0.5,1],
	[-0.5,1],
	[-1,0],
	[-0.5,-1],
	[0.5,-1],
	[1.5,-1],
	[2,0],
	[1.5,1], 
	[1,2], 
	[0,2], 
	[-1,2], 
	[-1.5,1], 
	[-2,0],
	[-1.5,-1],
	[-1,-2], 
	[0,-2], 
	[1,-2]]),
	'x': np.array([ 3.]),
	'y': np.array([ 3.])}
	n_centers = grid['centers']
	x, y = grid['x'], grid['y']
	# Size of figure in inches
	xinch = (x * w / y) / dpi
	yinch = (y * w / x) / dpi
	fig = plt.figure(figsize=(xinch, yinch), dpi=dpi)
	ax = fig.add_subplot(111, aspect='equal')
	# Get pixel size between to data points
	xpoints = n_centers[:, 0]
	ypoints = n_centers[:, 1]
	ax.scatter(xpoints, ypoints, s=0.0, marker='s')
	ax.axis([min(xpoints)-1., max(xpoints)+1.,
	         min(ypoints)-1., max(ypoints)+1.])
	xy_pixels = ax.transData.transform(np.vstack([xpoints, ypoints]).T)
	xpix, ypix = xy_pixels.T

	# In matplotlib, 0,0 is the lower left corner, whereas it's usually the
	# upper right for most image software, so we'll flip the y-coords
	width, height = fig.canvas.get_width_height()
	ypix = height - ypix

	# discover radius and hexagon
	apothem = .9 * (xpix[1] - xpix[0]) / math.sqrt(3)
	area_inner_circle = math.pi * (apothem ** 2)
	collection_bg = RegularPolyCollection(
	    numsides=6,  # a hexagon
	    rotation=0,
	    sizes=(area_inner_circle,),
	    edgecolors = (0, 0, 0, 1),
	    array= d_matrix,
	    cmap = 'inferno',
	    offsets = n_centers,
	    transOffset = ax.transData
	)
	I = ax.add_collection(collection_bg)
	I.set_clim(vmin=1000, vmax=1800)
	ax.axis('off')
	ax.autoscale_view()
	ax.set_title(title)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="10%", pad=0.05)
	plt.colorbar(collection_bg, cax=cax)
	plt.savefig(outdir+'/'+str(pointing)+'_hitmap.png')

	return print(outdir+'/'+str(pointing)+'_hitmap.png')

def plot_RMFmatrix(RMFmatrix, original_energybound, new_energybound, outdir, outname, ext = 'png', colmap = 'cubehelix', energy_scale = 'log'):
	X, Y = np.meshgrid(original_energybound, new_energybound)

	fig = plt.figure(figsize=(7, 7))
	ax = fig.add_subplot(111)
	I= ax.contourf(X, Y, RMFmatrix, cmap=colmap)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	plt.colorbar(I, cax=cax)
	ax.set_xlabel('E/keV')
	ax.set_ylabel('E/keV')
	ax.set_title('RMF')
	ax.set_xscale(energy_scale)
	ax.set_yscale(energy_scale)
	plt.savefig(outdir+'/'+outname+'.'+ext)
	return print('spectrum saved as '+outdir+'/'+outname+'.'+ext)

def plot_spectrum_noRMF(fitted, energy_boundaries, outdir, outname, ext = 'png'):
	plt.figure(figsize = (7,7))
	plt.errorbar(energy_boundaries, fitted[:,3], yerr = fitted[:,4], fmt = 'ko', marker = '+', capsize = 0)
	plt.ylim([1E-6, 2E-2])
	plt.xlim([27, 2E3])
	plt.yscale('log')
	plt.xscale('log')
	plt.savefig(outdir+'/'+outname+'.'+ext)

	return print('spectrum saved as '+outdir+'/'+outname+'.'+ext)

def plot_spectrumfit_noRMF(fitted, energy_boundaries, fitparams, outdir, outname, ext = 'png'):
	powerlaw = lambda x, amp, index: amp * (x**index)
	amp = fitparams[1]
	index = fitparams[0]
	energy_arr = np.linspace(min(energy_boundaries), max(energy_boundaries), 100)
	plt.figure(figsize = (7,7))
	plt.errorbar(energy_boundaries, fitted[:,3], yerr = fitted[:,4], fmt = 'ko', marker = '+', capsize = 0, label = 'data')
	plt.plot(energy_arr, powerlaw(energy_arr, amp, index), 'r', label = 'Powerlaw best fit')
	plt.ylim([9E-7, 2E-2])
	plt.xlim([27, 2E3])
	plt.xlabel('E/eV')
	plt.ylabel('$ph/s/cm^2/keV$')
	plt.title('Non-RMF corrected spectrum')
	plt.yscale('log')
	plt.xscale('log')
	plt.legend(loc = 'best')
	plt.annotate(r'$\alpha$: %.3g  $A_0$: %.3g' % (index, amp), xy=(40, 3E-6), xytext=(40, 3E-6),
            arrowprops=None)
	plt.savefig(outdir+'/'+outname+'.'+ext)

	return print('spectrum saved as '+outdir+'/'+outname+'.'+ext)