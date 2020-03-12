"""
What it does
------------

Add all 4MOST columns

Link to template

redden magnitudes...

add flag is_in_footprint

Command to run
--------------

python3 005_3_4most_catalog.py

then concatenate catalogues over the footprint de <10 ?

"""

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from sklearn.neighbors import BallTree

from astropy_healpix import healpy
import sys
import os
import time
from scipy.interpolate import interp1d
from scipy.stats import norm
from astropy.table import Table, Column
from scipy.optimize import curve_fit
import linmix
import astropy.io.fits as fits
import h5py
import numpy as n
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
import extinction

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(
	H0=67.77 * u.km / u.s / u.Mpc,
	Om0=0.307115)  # , Ob0=0.048206)
h = 0.6777
L_box = 1000.0 / h
cosmo = cosmoMD

zs = n.arange(0.00001, 2.1, 0.0001)
dm_itp = interp1d(zs, cosmo.distmod(zs).value)

#import astropy.io.fits as fits
# import all pathes

env = sys.argv[1]  # 'MD04'
print(env)

root_dir = os.path.join(os.environ[env])
dir_2_OUT = os.path.join(root_dir, "cat_SHAM_COSMO")

sub_survey_names = n.array([ 'BG', 'LRG', 'ELG', 'QSO', 'LyA', 'filament_GAL'])

N_subsurvey = {'BG':1, 'filament_GAL':3, 'LRG':2, 'ELG':3, 'QSO':4, 'LyA':5}
priority_values = {'BG':100, 'filament_GAL':80, 'LRG':99, 'ELG':80, 'QSO':97, 'LyA':98}

# reassigns templates correctly
z_all = n.hstack(( 0., n.arange(0.3, 3., 0.2), 3.5, 4.5, 6. ))
zmins = z_all[:-1]
zmaxs = z_all[1:]


N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels)[::-1]:
	print(HEALPIX_id)
	path_2_BG    = os.path.join(dir_2_OUT, 'BG_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_BG_S5 = os.path.join(dir_2_OUT, 'S5GAL_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_LRG   = os.path.join(dir_2_OUT, 'LRG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_ELG   = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_QSO   = os.path.join(dir_2_OUT, 'QSO_' + str(HEALPIX_id).zfill(6) + '.fit')

	t_bg = Table.read(path_2_BG)
	t_bgS5 = Table.read(path_2_BG_S5)
	t_lrg = Table.read(path_2_LRG)
	t_elg = Table.read(path_2_ELG)
	#t_qso_i = Table.read(path_2_QSO)
	#t_qso = t_qso_i[(t_qso_i['Z']<=2.2) & (t_qso_i['Z']>0.8) ]
	#t_lya = t_qso_i[(t_qso_i['Z']> 2.2) & (t_qso_i['Z']<=3.5)]

	def add_4most_columns(subsurvey = 'BG', t_survey = t_bg ):
		#subsurvey = 'BG'
		#t_survey = t_bg 
		#subsurvey = 'LRG'
		#t_survey = t_lrg 
		#subsurvey = 'ELG'
		#t_survey = t_elg 
		#subsurvey = 'QSO'
		#t_survey = t_qso 
		#subsurvey = 'Lya'
		#t_survey = t_lya 
		#subsurvey = 'filament_GAL'
		#t_survey = t_bgS5 
		N_obj = len(t_survey)
		#  limit size of the string columns to the size of the longer string in the corresponding columns. 
		# 'NAME':str, max 256 char
		N1 = n.arange(len(t_survey['EBV']))
		id_list = HEALPIX_id*1e8 + N_subsurvey[subsurvey]*1e6 + N1
		NAME = n.array([ str(int(el)).zfill(11) for el in id_list ])
		t_survey.add_column(Column(name='NAME', data=NAME, unit=''))
		# 'RA':n.float64, 1D
		# 'DEC':n.float64, 1D
		# 'PMRA':n.float32, 1E
		# 'PMDEC':n.float32, 1E
		# 'EPOCH':n.float32, 1E
		PMRA = n.zeros(N_obj)
		t_survey.add_column(Column(name='PMRA', data=PMRA, unit='mas/yr'))
		PMDEC = n.zeros(N_obj)
		t_survey.add_column(Column(name='PMDEC', data=PMDEC, unit='mas/yr'))
		EPOCH = n.ones(N_obj)*2000.
		t_survey.add_column(Column(name='EPOCH', data=EPOCH, unit='yr'))
		# 'RESOLUTION':n.int16, 1I
		RESOLUTION = n.ones(N_obj).astype('int')
		t_survey.add_column(Column(name='RESOLUTION', data=RESOLUTION, unit=''))
		# 'SUBSURVEY':str, max 256 char
		SUBSURVEY = n.ones(N_obj).astype('str')
		SUBSURVEY[:] = subsurvey
		t_survey.add_column(Column(name='SUBSURVEY', data=SUBSURVEY, unit=''))
		# 'PRIORITY':n.int16, 1I
		PRIORITY = n.zeros(N_obj).astype('int') + priority_values[subsurvey]
		t_survey.add_column(Column(name='PRIORITY', data=PRIORITY, unit=''))
		# EBV for templates
		ebv_1000 = (t_survey['EBV']*1000).astype('int')
		#print('EBV', n.min(ebv_1000), n.max(ebv_1000))
		ebv_1_0 = ( ebv_1000 > 1000 ) 
		ebv_0_5 = ( ebv_1000 > 500 ) & ( ebv_1000 <= 1000 ) 
		ebv_0_4 = ( ebv_1000 > 400 ) & ( ebv_1000 <= 500 ) 
		ebv_0_3 = ( ebv_1000 > 300 ) & ( ebv_1000 <= 400 ) 
		ebv_0_2 = ( ebv_1000 > 200 ) & ( ebv_1000 <= 300 ) 
		ebv_0_1 = ( ebv_1000 > 100 ) & ( ebv_1000 <= 200 ) 
		ebv_0_0 = ( ebv_1000 <= 100 ) 
		z_name = lambda z0, z1 : "_zmin_"+str(int(10*z0)).zfill(2)+"_zmax_"+str(int(10*z1)).zfill(2)
		# templates
		template_names = n.zeros(N_obj).astype('U100')
		ruleset_array = n.zeros(N_obj).astype('str')
		# S8 BG or LRG
		if subsurvey == 'BG' or subsurvey == 'LRG':
			ruleset_array[:] = "COSMO_RedGAL"
			for z0,z1 in zip(zmins,zmaxs):
				zsel = (t_survey['Z']>=z0) & (t_survey['Z']<z1)
				if len(zsel.nonzero()[0])>0:
					#ruleset_array[zsel] = "COSMO_RedGAL"
					template_names[(zsel)]               = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_01.fits'  
					template_names[(zsel)&(ebv_0_0)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_01.fits'  
					template_names[(zsel)&(ebv_0_1)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_1.fits'   
					template_names[(zsel)&(ebv_0_2)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_2.fits'   
					template_names[(zsel)&(ebv_0_3)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_3.fits'   
					template_names[(zsel)&(ebv_0_4)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_4.fits'   
					template_names[(zsel)&(ebv_0_5)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_5.fits'   
					template_names[(zsel)&(ebv_1_0)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_1_0.fits'   
		# S5 filament gal
		if subsurvey == 'filament_GAL' :
			ruleset_array[:] = "RedGAL"
			for z0,z1 in zip(zmins,zmaxs):
				zsel = (t_survey['Z']>=z0) & (t_survey['Z']<z1)
				if len(zsel.nonzero()[0])>0:
					#ruleset_array[zsel] = "RedGAL"
					template_names[(zsel)]               = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_01.fits'  
					template_names[(zsel)&(ebv_0_0)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_01.fits'  
					template_names[(zsel)&(ebv_0_1)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_1.fits'   
					template_names[(zsel)&(ebv_0_2)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_2.fits'   
					template_names[(zsel)&(ebv_0_3)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_3.fits'   
					template_names[(zsel)&(ebv_0_4)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_4.fits'   
					template_names[(zsel)&(ebv_0_5)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_0_5.fits'   
					template_names[(zsel)&(ebv_1_0)]     = "4most_"+'LRG'+z_name( z0, z1)+'_EBV_1_0.fits'   
		# S8 ELG
		if subsurvey == 'ELG' :
			ruleset_array[:] = "ELG"
			for z0,z1 in zip(zmins,zmaxs):
				zsel = (t_survey['Z']>=z0)&(t_survey['Z']<z1)
				if len(zsel.nonzero()[0])>0:
					#ruleset_array[zsel] = "ELG"
					template_names[(zsel)]               = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_01.fits' 
					template_names[(zsel)&(ebv_0_0)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_01.fits' 
					template_names[(zsel)&(ebv_0_1)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_1.fits'  
					template_names[(zsel)&(ebv_0_2)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_2.fits'  
					template_names[(zsel)&(ebv_0_3)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_3.fits'  
					template_names[(zsel)&(ebv_0_4)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_4.fits'  
					template_names[(zsel)&(ebv_0_5)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_0_5.fits'  
					template_names[(zsel)&(ebv_1_0)]     = "4most_"+'ELG'+z_name( z0, z1)+'_EBV_1_0.fits'  

		# 'TEMPLATE':str, max 256 char
		t_survey.add_column(Column(name='TEMPLATE', data=template_names, unit=''))
		# 'RULESET':str, max 256 char
		t_survey.add_column(Column(name='RULESET', data=ruleset_array, unit=''))
		# 'REDSHIFT_ESTIMATE':n.float32, 1E
		# 'REDSHIFT_ERROR':n.float32, 1E
		t_survey.add_column(Column(name='REDSHIFT_ESTIMATE', data=t_survey['Z'], unit=''))
		t_survey.add_column(Column(name='REDSHIFT_ERROR', data=n.ones(N_obj), unit=''))
		# 'EXTENT_FLAG': 1I
		# =1
		# 'EXTENT_PARAMETER': 1E
		# =0
		# 'EXTENT_INDEX': 1E
		# =0
		t_survey.add_column(Column(name='EXTENT_FLAG'     , data=n.ones(N_obj).astype('int') , unit=''))
		t_survey.add_column(Column(name='EXTENT_PARAMETER', data=n.zeros(N_obj), unit=''))
		t_survey.add_column(Column(name='EXTENT_INDEX'    , data=n.zeros(N_obj), unit=''))
		# 'MAG':n.float32,
		# 'MAG_ERR':n.float32
		# 'MAG_TYPE': str max 256 char
		r_v=3.1
		a_v = t_survey['EBV'] * r_v
		delta_mag = n.hstack(( n.array([ extinction.fitzpatrick99(n.array([6500.]), el, r_v=3.1, unit='aa') for el in a_v ]) ))
		#rv = av/ebv
		#av = rv x ebv
		extincted_mag = t_survey['rfib'] + delta_mag
		t_survey.add_column(Column(name='MAG', data=extincted_mag, unit='mag'))
		t_survey.add_column(Column(name='MAG_ERR', data=0.01 * n.ones(N_obj), unit='mag'))
		MAG_TYPE = n.ones(N_obj).astype('str')
		MAG_TYPE[:] = 'DECam_r_AB'
		t_survey.add_column(Column(name='MAG_TYPE', data=MAG_TYPE, unit=''))
		# 'REDDENING':n.float32, 1E
		t_survey.add_column(Column(name='REDDENING',data=t_survey['EBV'], unit='mag'))
		# 'DATE_EARLIEST':n.float64, JulianDate decimal days # 01-Nov-2022
		# 'DATE_LATEST':n.float64, JulianDate decimal days # 02-Feb-2033
		t_survey.add_column(Column(name='DATE_EARLIEST',data=22305 * n.ones(N_obj), unit='d'))
		t_survey.add_column(Column(name='DATE_LATEST'  ,data=33033 * n.ones(N_obj), unit='d'))
		return t_survey

	#path_2_BG    = os.path.join(dir_2_OUT, 'BG_'  + str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_BG_S5 = os.path.join(dir_2_OUT, 'S5GAL_'  + str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_LRG   = os.path.join(dir_2_OUT, 'LRG_' + str(HEALPIX_id).zfill(6) + '.fit')
	#path_2_ELG   = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')
	#
	subsurvey = 'BG'
	t_survey = t_bg 
	t_out = add_4most_columns(subsurvey = subsurvey, t_survey = t_survey)
	t_out.write (path_2_BG  , overwrite=True)
	
	subsurvey = 'LRG'
	t_survey = t_lrg 
	t_out = add_4most_columns(subsurvey = subsurvey, t_survey = t_survey)
	t_out.write (path_2_LRG  , overwrite=True)
	
	subsurvey = 'ELG'
	t_survey = t_elg 
	t_out = add_4most_columns(subsurvey = subsurvey, t_survey = t_survey)
	t_out.write (path_2_ELG  , overwrite=True)
	
	subsurvey = 'filament_GAL'
	t_survey = t_bgS5 
	t_out = add_4most_columns(subsurvey = subsurvey, t_survey = t_survey)
	t_out.write (path_2_BG_S5  , overwrite=True)
