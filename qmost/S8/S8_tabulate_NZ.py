import os, sys
import numpy as np
import healpy
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from matplotlib.path import Path 

survey = 'S8'
working_dir = os.path.join(os.environ['HOME'], 'data', '4most', survey)
# name of the sub survey for which the LSM is created
sub_survey_names = np.array([ 'BG', 'LRG', 'ELG', 'QSO', 'LyA'])
# loads the catalog
data = fits.open(os.path.join(working_dir, 'S8_Cosmology_latest_input_cat_4FS_new_format.fits.gz'))[1].data

selections = np.array([(data['SUBSURVEY']==sub_survey_name) for sub_survey_name in sub_survey_names ])

def tabulate_NZ(sub_survey_name, selection):
	out_file = os.path.join(working_dir, survey +'_NZ_'+str(sub_survey_name)+'.fits')
	plot_file = os.path.join(working_dir,survey +'_NZ_'+str(sub_survey_name)+'.png')
	DZ=0.1
	z_bins = np.arange(0,6.1,DZ)
	xx = 0.5*(z_bins[:-1]+z_bins[1:])
	#########
	#########
	# 1. redshift histogram
	#########
	#########
	# create the ARRAY OF PIXEL INDEXES
	ZZ = data['REDSHIFT_ESTIMATE']  [selection]
	NN = np.histogram(ZZ, bins = z_bins)[0]
	#########
	#########
	# 3. Plots the LSM values over RA and DEC
	#########
	#########
	plt.figure(0, (6,6))
	plt.errorbar(xx, NN, xerr=DZ/2., yerr=NN**(-0.5))
	plt.grid()
	plt.ylim((0.9, 1.2*np.max(NN) ))
	plt.xlim((0.0, 1.2*np.max(xx[NN>=1]) ))
	plt.xlabel('redshift')
	plt.ylabel('Counts (dz=0.1)')
	plt.yscale('log')
	plt.title(sub_survey_name+' '+str(np.sum(NN)))
	plt.savefig(plot_file)
	plt.clf()
	#########
	#########
	# 4. Create the fits LSM file
	#########
	#########
	cols = fits.ColDefs([
		fits.Column( "z_min"	 ,unit='' ,format="D", array=z_bins[:-1]), 
		fits.Column( "z_max"     ,unit='' ,format="D", array=z_bins[1:] ), 
		fits.Column( "z_middle"	 ,unit='' ,format="D", array=xx), 
		fits.Column( "counts"    ,unit='' ,format="D", array=NN)
		])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.header['author'] = 'JC'
	tbhdu.header['HIERARCH SUBSURVEY'] = sub_survey_name
	if os.path.isfile(out_file):
		os.remove(out_file)
	tbhdu.writeto(out_file)

for sub_survey_name, selection in zip(sub_survey_names,selections):
	print(sub_survey_name)
	tabulate_NZ(sub_survey_name, selection)