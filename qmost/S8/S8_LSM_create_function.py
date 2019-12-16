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
data = fits.open(os.path.join(working_dir, 'S8_Cosmology_May2nd.fits.gz'))[1].data

selections = np.array([(data['SUBSURVEY']==sub_survey_name) for sub_survey_name in sub_survey_names ])

DESfile = os.path.join(working_dir, 'DES_footprint.txt')
DESra, DESdec = np.loadtxt(DESfile, unpack=True)
DESnedge= len(DESra)
# isinpoly
DESpoly = Path(np.concatenate(
                    (DESra.reshape((DESnedge,1)),
                     DESdec.reshape((DESnedge,1))),
                    axis=1))

def create_LSM(sub_survey_name, selection):
	print(sub_survey_name, selection)
	out_file = os.path.join(working_dir, survey +'_LSM_'+str(sub_survey_name)+'.fits')
	plot_file = os.path.join(working_dir,survey +'_LSM_'+str(sub_survey_name)+'.png')
	#########
	#########
	# 1. First create the grid of pixels over the full sky
	#########
	#########
	# create the ARRAY OF PIXEL INDEXES
	NSIDE=32
	N_pixels_all = healpy.nside2npix(NSIDE)
	pix_ids = np.arange(N_pixels_all)
	# ARRAY OF RIGHT ASCENSIONS OF THE PIXEL CENTERS
	ra_hp = np.array([healpy.pix2ang(NSIDE,pix_id, nest=True)[1]*180./np.pi for pix_id in pix_ids ])
	# ARRAY OF DECLINATION OF THE PIXEL CENTERS
	dec_hp = np.array([ (np.pi/2. - healpy.pix2ang(NSIDE,pix_id, nest=True)[0])*180./np.pi for pix_id in pix_ids ])
	# area subtended by each pixel in deg2
	area_per_pixel =  healpy.nside2pixarea(32)*(180/np.pi)**2  
	# Initialize all LSM values 
	lsm_values = np.zeros(len(ra_hp))
	#########
	#########
	# 2. Retrieve pixels that contain targets
	#########
	#########
	# select coordinates in your catalog file
	RA = data['RA']  [selection]
	DEC = data['DEC'][selection]
	# compute healpix index for each target
	DATA_HEALPIX_32   = healpy.ang2pix(32,   np.pi/2. - DEC*np.pi/180. , RA*np.pi/180. , nest=True)
	# retrieves a unique and sorted list of healpix indexes containing targets
	unique_healpix_32_indexes = np.array(list(set(DATA_HEALPIX_32)))
	unique_healpix_32_indexes.sort()
	#########
	#########
	# 2. Assigns LSM values
	#########
	#########
	# DES footprint
	ra_hp2 = ra_hp[unique_healpix_32_indexes]
	dec_hp2 = dec_hp[unique_healpix_32_indexes]
	tmpradec      = np.zeros((len(ra_hp2),2))

	tmpradec[:,0] = ra_hp2
	tmp = (tmpradec[:,0]>270)
	tmpradec[tmp,0]-= 360.
	tmpradec[:,1] = dec_hp2
	des_survey           = DESpoly.contains_points(tmpradec)

	des_survey    = DESpoly.contains_points(tmpradec)
	unique_healpix_32_indexes_in_DES = unique_healpix_32_indexes[des_survey]
	# 
	N_pixels_with_targets = len(unique_healpix_32_indexes)
	N_pixels_with_targets_in_DES = len(unique_healpix_32_indexes_in_DES)
	print('N pixels outside DES, and inside',N_pixels_with_targets, N_pixels_with_targets_in_DES)
	print('Area outside DES, and inside',N_pixels_with_targets*area_per_pixel, N_pixels_with_targets_in_DES*area_per_pixel)
	norm_per_pixel = 1./(area_per_pixel*N_pixels_with_targets)
	norm_high = norm_per_pixel * 1.2
	norm_low = (1. - norm_high*N_pixels_with_targets_in_DES)/(N_pixels_with_targets-N_pixels_with_targets_in_DES)

	lsm_values = np.zeros(len(ra_hp))
	print('lsm values',norm_high, norm_low)
	lsm_values[unique_healpix_32_indexes] = norm_high
	lsm_values[unique_healpix_32_indexes_in_DES] = norm_low
	print('unique lsm values',np.unique(lsm_values))

	print('sum of LSM values:',np.sum(lsm_values))
	#########
	#########
	# 3. Plots the LSM values over RA and DEC
	#########
	#########
	plt.figure(0, (6,6))
	ok = lsm_values>0
	#plt.scatter(ra_hp, dec_hp, c=np.log10(lsm_values), s=4, marker='s', edgecolors='face')
	plt.scatter(ra_hp[ok], dec_hp[ok], c=np.log10(lsm_values[ok]), s=4, marker='s', edgecolors='face')
	cb = plt.colorbar(shrink=0.8)
	cb.set_label('LSM')
	plt.grid()
	plt.ylim((-90,40))
	plt.xlim((-1,361))
	plt.xlabel('R.A.')
	plt.ylabel('Dec')
	plt.title(sub_survey_name)
	plt.savefig(plot_file)
	plt.clf()
	#########
	#########
	# 4. Create the fits LSM file
	#########
	#########
	cols = fits.ColDefs([
		fits.Column( "healpix_id"	 ,unit='' ,format="K", array=pix_ids), 
		fits.Column( "RA"            ,unit='deg' ,format="D", array=ra_hp ), 
		fits.Column( "DEC"	         ,unit='deg' ,format="D", array=dec_hp), 
		fits.Column( "LSM"           ,unit='' ,format="D", array=lsm_values)
		])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.header['author'] = 'JC'
	tbhdu.header['HIERARCH SUBSURVEY'] = sub_survey_name
	tbhdu.header['NSIDE'] = NSIDE
	tbhdu.header['COORD'] = 'equatorial'
	tbhdu.header['NESTED'] = 'True'
	if os.path.isfile(out_file):
		os.remove(out_file)
	tbhdu.writeto(out_file)

for sub_survey_name, selection in zip(sub_survey_names,selections):
	print(sub_survey_name)
	create_LSM(sub_survey_name, selection)