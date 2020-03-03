"""
What it does
------------

Add realistic magnitudes to BG, LRG and ELG

Command to run
--------------

python3 005_2_all_magnitudes.py environmentVAR 

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

env = 'MD10' # sys.argv[1]  # 'MD04'
print(env)

stilts_cmd = 'stilts'

root_dir = os.path.join(os.environ[env])

path_2_kids = os.path.join(root_dir, '4most_s58_kidsdr4.fits') 
# '/home/comparat/data/4most/S8/4most_s58_kidsdr4.fits'
# parent data :
t_kids = Table.read(path_2_kids)
#good_phot = ( t_kids['rtot']>0 ) & (t_kids['gr'] > -10 ) & (t_kids['ri'] > -10 )  & (t_kids['iz'] > -10 )  & (t_kids['zy'] > -10 )  & (t_kids['yj'] > -10 )  & (t_kids['jh'] > -10 )  & (t_kids['hks'] > -10 ) & (t_kids['zphot']>0)

plotDir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'logMs-MK-fits')
dir_2_OUT = os.path.join(root_dir, "cat_SHAM_COSMO")

def add_all_MAG_ELG_only(HEALPIX_id):
	path_2_ELG   = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')
	t_elg = Table.read(path_2_ELG)
	def add_magnitudes(t_survey = t_elg, flag = 'iselg', t_kids_i=t_kids):
		sim_redshift = t_survey['Z']
		sim_sfr = (t_survey['SFR'] - t_survey['SFR'].min() ) / ( t_survey['SFR'].max() - t_survey['SFR'].min() )
		NN, bb = n.histogram(sim_sfr, bins=1000)
		cdf_itp = interp1d(bb, n.hstack((0,n.cumsum(NN)))/n.sum(NN) )
		sim_sfr_0_1 = cdf_itp(sim_sfr)
		#
		sim_k_mag_i = t_survey['K_mag_abs']+dm_itp(t_survey['Z'])
		sim_k_mag = (sim_k_mag_i - sim_k_mag_i.min()) / (sim_k_mag_i.max() - sim_k_mag_i.min())
		NN, bb = n.histogram(sim_k_mag, bins=1000)
		cdf_itpK = interp1d(bb, n.hstack((0,n.cumsum(NN)))/n.sum(NN) )
		sim_k_mag_0_1 = cdf_itpK(sim_k_mag)
		# normalize data arrays to the min max of the simulation
		s_kids = (t_kids_i[flag]) & (t_kids_i['zphot']<=n.max(sim_redshift)) & (t_kids_i['zphot']>=n.min(sim_redshift)) 
		t_kids2 = t_kids_i[s_kids]
		# K magnitude
		kids_k_mag_i = t_kids2['rtot']-(t_kids2['ri']+t_kids2['iz']+t_kids2['zy']+t_kids2['yj']+t_kids2['jh']+t_kids2['hks'])
		kids_k_mag = (kids_k_mag_i - kids_k_mag_i.min()) / (kids_k_mag_i.max() - kids_k_mag_i.min())
		NN, bb = n.histogram(kids_k_mag, bins=1000)
		cdf_itpINV_kids_KMAG = interp1d(n.hstack((0,n.cumsum(NN)))/n.sum(NN), bb)
		# REDSHIFT
		kids_photoz = t_kids2['zphot']
		# SFR / COLOR
		kids_sfr =  (t_kids2['gr'] - t_kids2['gr'].min() ) / ( t_kids2['gr'].max() - t_kids2['gr'].min() )
		NN, bb = n.histogram(kids_sfr, bins=1000)
		cdf_itpINV_kids_SFR = interp1d(n.hstack((0,n.cumsum(NN)))/n.sum(NN), bb)
		# input kids into a tree
		Tree_kids = BallTree(n.transpose([kids_photoz, kids_k_mag, kids_sfr]))
		sim_k_mag_NORM = cdf_itpINV_kids_KMAG(sim_k_mag_0_1)
		sim_sfr_NORM = cdf_itpINV_kids_SFR(sim_sfr_0_1)
		DATA = n.transpose([sim_redshift, sim_k_mag_NORM, sim_sfr_NORM])
		kids_ID = n.arange(len(t_kids2['rtot']))
		#
		test = Tree_kids.query(DATA, k=1, return_distance = True)
		ids_out = Tree_kids.query(DATA, k=1, return_distance = False)
		ids = n.hstack((ids_out))
		id_to_map = kids_ID[ids]
		if 'rtot' in t_survey.columns :
			t_survey['rtot'] = t_kids2['rtot'][id_to_map]
			t_survey['rfib'] = t_kids2['rfib'][id_to_map]
			t_survey['ug'  ] = t_kids2['ug'][id_to_map]
			t_survey['gr'  ] = t_kids2['gr'][id_to_map]
			t_survey['ri'  ] = t_kids2['ri'][id_to_map]
			t_survey['iz'  ] = t_kids2['iz'][id_to_map]
			t_survey['zy'  ] = t_kids2['zy'][id_to_map]
			t_survey['yj'  ] = t_kids2['yj'][id_to_map]
			t_survey['jh'  ] = t_kids2['jh'][id_to_map]
			t_survey['hks' ] = t_kids2['hks'][id_to_map]
		else:
			t_survey.add_column(Column(name='rtot', data=t_kids2['rtot'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='rfib', data=t_kids2['rfib'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='ug',data=t_kids2['ug'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='gr',data=t_kids2['gr'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='ri', data=t_kids2['ri'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='iz', data=t_kids2['iz'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='zy', data=t_kids2['zy'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='yj', data=t_kids2['yj'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='jh', data=t_kids2['jh'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='hks', data=t_kids2['hks'][id_to_map], unit='mag'))
		return t_survey #.write(path_2_out)

	t_out = add_magnitudes(t_survey = t_elg, flag = 'iss8elg', t_kids_i = t_kids)
	t_out.write (path_2_ELG  , overwrite=True)

def add_all_MAG(HEALPIX_id):
	#HEALPIX_id = 0 
	path_2_BG    = os.path.join(dir_2_OUT, 'BG_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_BG_S5 = os.path.join(dir_2_OUT, 'S5GAL_'  + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_LRG   = os.path.join(dir_2_OUT, 'LRG_' + str(HEALPIX_id).zfill(6) + '.fit')
	path_2_ELG   = os.path.join(dir_2_OUT, 'ELG_' + str(HEALPIX_id).zfill(6) + '.fit')

	t_bg = Table.read(path_2_BG)
	t_bgS5 = Table.read(path_2_BG_S5)
	t_lrg = Table.read(path_2_LRG)
	t_elg = Table.read(path_2_ELG)

	def add_magnitudes(t_survey = t_elg, flag = 'iselg', t_kids_i=t_kids):
		#t_survey = t_elg
		#flag = 'iss8elg'
		#t_kids_i=t_kids
		# normalize simulated arrays to 0-1
		sim_redshift = t_survey['Z']
		sim_sfr = (t_survey['SFR'] - t_survey['SFR'].min() ) / ( t_survey['SFR'].max() - t_survey['SFR'].min() )
		NN, bb = n.histogram(sim_sfr, bins=1000)
		cdf_itp = interp1d(bb, n.hstack((0,n.cumsum(NN)))/n.sum(NN) )
		sim_sfr_0_1 = cdf_itp(sim_sfr)

		sim_k_mag_i = t_survey['K_mag_abs']+dm_itp(t_survey['Z'])
		sim_k_mag = (sim_k_mag_i - sim_k_mag_i.min()) / (sim_k_mag_i.max() - sim_k_mag_i.min())
		NN, bb = n.histogram(sim_k_mag, bins=1000)
		cdf_itpK = interp1d(bb, n.hstack((0,n.cumsum(NN)))/n.sum(NN) )
		sim_k_mag_0_1 = cdf_itpK(sim_k_mag)
		# normalize data arrays to the min max of the simulation
		s_kids = (t_kids_i[flag]) & (t_kids_i['zphot']<=n.max(sim_redshift)) & (t_kids_i['zphot']>=n.min(sim_redshift)) 
		t_kids2 = t_kids_i[s_kids]
		# K magnitude
		kids_k_mag_i = t_kids2['rtot']-(t_kids2['ri']+t_kids2['iz']+t_kids2['zy']+t_kids2['yj']+t_kids2['jh']+t_kids2['hks'])
		kids_k_mag = (kids_k_mag_i - kids_k_mag_i.min()) / (kids_k_mag_i.max() - kids_k_mag_i.min())
		NN, bb = n.histogram(kids_k_mag, bins=1000)
		cdf_itpINV_kids_KMAG = interp1d(n.hstack((0,n.cumsum(NN)))/n.sum(NN), bb)
		# REDSHIFT
		kids_photoz = t_kids2['zphot']
		# SFR / COLOR
		kids_sfr =  (t_kids2['gr'] - t_kids2['gr'].min() ) / ( t_kids2['gr'].max() - t_kids2['gr'].min() )
		NN, bb = n.histogram(kids_sfr, bins=1000)
		cdf_itpINV_kids_SFR = interp1d(n.hstack((0,n.cumsum(NN)))/n.sum(NN), bb)
		# input kids into a tree
		Tree_kids = BallTree(n.transpose([kids_photoz, kids_k_mag, kids_sfr]))
		sim_k_mag_NORM = cdf_itpINV_kids_KMAG(sim_k_mag_0_1)
		sim_sfr_NORM = cdf_itpINV_kids_SFR(sim_sfr_0_1)
		DATA = n.transpose([sim_redshift, sim_k_mag_NORM, sim_sfr_NORM])
		kids_ID = n.arange(len(t_kids2['rtot']))
		"""
		p.figure(2, (6,6))
		p.hist(sim_redshift , bins=100, normed=True, cumulative=True, histtype='step', label='sim' )
		p.hist(kids_photoz , bins=100, normed=True, cumulative=True, histtype='step', label='data' )
		p.ylabel('cdf')
		p.xlabel(r'redshift')
		p.grid()
		p.legend(frameon=False, loc=0)
		p.savefig(os.path.join( plotDir, "3d_hist_REDSHIFT.png"))
		p.clf()

		p.figure(2, (6,6))
		p.hist(sim_sfr_NORM , bins=100, normed=True, cumulative=True, histtype='step', label='sim' )
		p.hist(kids_sfr , bins=100, normed=True, cumulative=True, histtype='step', label='data' )
		p.ylabel('cdf')
		p.xlabel(r'SFR proxy')
		p.grid()
		p.legend(frameon=False, loc=0)
		p.savefig(os.path.join( plotDir, "3d_hist_SFR.png"))
		p.clf()

		p.figure(2, (6,6))
		p.hist(sim_k_mag_NORM , bins=100, normed=True, cumulative=True, histtype='step', label='sim' )
		p.hist(kids_k_mag , bins=100, normed=True, cumulative=True, histtype='step', label='data' )
		p.ylabel('cdf')
		p.xlabel(r'kmag')
		p.grid()
		p.legend(frameon=False, loc=0)
		p.savefig(os.path.join( plotDir, "3d_hist_Kmag.png"))
		p.clf()
		"""
		#Tree_sim = BallTree(DATA)
		test = Tree_kids.query(DATA, k=1, return_distance = True)
		ids_out = Tree_kids.query(DATA, k=1, return_distance = False)
		ids = n.hstack((ids_out))
		id_to_map = kids_ID[ids]
		if 'rtot' in t_survey.columns :
			t_survey['rtot'] = t_kids2['rtot'][id_to_map]
			t_survey['rfib'] = t_kids2['rfib'][id_to_map]
			t_survey['ug'  ] = t_kids2['ug'][id_to_map]
			t_survey['gr'  ] = t_kids2['gr'][id_to_map]
			t_survey['ri'  ] = t_kids2['ri'][id_to_map]
			t_survey['iz'  ] = t_kids2['iz'][id_to_map]
			t_survey['zy'  ] = t_kids2['zy'][id_to_map]
			t_survey['yj'  ] = t_kids2['yj'][id_to_map]
			t_survey['jh'  ] = t_kids2['jh'][id_to_map]
			t_survey['hks' ] = t_kids2['hks'][id_to_map]
		else:
			t_survey.add_column(Column(name='rtot', data=t_kids2['rtot'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='rfib', data=t_kids2['rfib'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='ug',data=t_kids2['ug'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='gr',data=t_kids2['gr'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='ri', data=t_kids2['ri'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='iz', data=t_kids2['iz'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='zy', data=t_kids2['zy'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='yj', data=t_kids2['yj'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='jh', data=t_kids2['jh'][id_to_map], unit='mag'))
			t_survey.add_column(Column(name='hks', data=t_kids2['hks'][id_to_map], unit='mag'))
		return t_survey #.write(path_2_out)

	t_out = add_magnitudes(t_survey = t_elg, flag = 'iss8elg', t_kids_i = t_kids)
	t_out.write (path_2_ELG  , overwrite=True)

	t_out = add_magnitudes(t_survey = t_lrg, flag = 'iss8lrg', t_kids_i = t_kids)
	t_out.write (path_2_LRG  , overwrite=True)

	t_out = add_magnitudes(t_survey = t_bg, flag = 'iss8bg', t_kids_i = t_kids)
	t_out.write (path_2_BG  , overwrite=True)

	t_out = add_magnitudes(t_survey = t_bgS5, flag = 'isR195', t_kids_i = t_kids)
	t_out.write (path_2_BG_S5  , overwrite=True)


N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels):
	print(HEALPIX_id)
	add_all_MAG(HEALPIX_id)
	#add_all_MAG_ELG_only(HEALPIX_id)


