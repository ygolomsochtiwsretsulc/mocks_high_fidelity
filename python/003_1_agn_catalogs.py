"""
What it does
------------

Creates a AGN fits catalog for each healpix 768 pixel of 13 deg2 (NSIDE=8)

Command to run
--------------

python3 003_1_agn_catalogs.py environmentVAR fileBasename

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib


"""
from astropy_healpix import healpy
from astropy.table import Table, Column
import sys, os, time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.io.fits as fits
#import h5py
import numpy as n
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits

# import all pathes

env = sys.argv[1]  # 'MD04'
baseName =  sys.argv[2]  # "sat_0.62840"
print(env, baseName)

test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_agn_file = os.path.join(test_dir, baseName + '_agn.fits')

dir_2_eRO_catalog = os.path.join(test_dir, 'cat_AGN_' + baseName)
if os.path.isdir(dir_2_eRO_catalog) == False:
	os.system('mkdir -p ' + dir_2_eRO_catalog)


# eRosita flux limits
path_2_flux_limits = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "erosita",
    "flux_limits.fits")


# simulation setup
if env[:2] == "MD" : # env == "MD04" or env == "MD40" or env == "MD10" or env == "MD25"
	from astropy.cosmology import FlatLambdaCDM
	import astropy.units as u
	cosmoMD = FlatLambdaCDM(
		H0=67.77 * u.km / u.s / u.Mpc,
		Om0=0.307115)  # , Ob0=0.048206)
	h = 0.6777
	L_box = 1000.0 / h
	cosmo = cosmoMD
if env[:4] == "UNIT" : # == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
	from astropy.cosmology import FlatLambdaCDM
	import astropy.units as u
	cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
	h = 0.6774
	L_box = 1000.0 / h
	cosmo = cosmoUNIT

f3 = fits.open(path_2_agn_file)[1].data
ids_active = f3["ids_active"]
print('agn file opened', time.time() - t0)

f0 = fits.open(path_2_light_cone)[1].data[ids_active]
print('halo file opened', time.time() - t0)

f1 = fits.open(path_2_galaxy_file)[1].data[ids_active]
print('galaxy file opened', time.time() - t0)

f2 = fits.open(path_2_coordinate_file)[1].data[ids_active]
ra = f2['ra']
dec = f2['dec']
print('coordinate file opened', time.time() - t0)

# eROSITA flux limit
pix_ids = healpy.ang2pix(
    512,
    n.pi /
    2. -
    f2['g_lat'] *
    n.pi /
    180.,
    f2['g_lon'] *
    n.pi /
    180.,
    nest=True)
flux_lim_data = fits.open(path_2_flux_limits)  # [1].data
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3
flux_limit = flux_lim_data[1].data['flux_limit_eRASS8']
print('flux limit file opened', time.time() - t0)

# selection function here :
# ( FX_soft_attenuated > 10**(flux_limit[pix_ids]-2) ) & ( SDSS_r_AB_attenuated < 26.5 )
sf1 = (f3['FX_soft_attenuated'] > 0)
print('selection function applied', time.time() - t0)

HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - dec*n.pi/180. , ra*n.pi/180. , nest=True)

for HEALPIX_8_id in n.arange(healpy.nside2npix(8)):
	path_2_eRO_catalog = os.path.join(dir_2_eRO_catalog, str(HEALPIX_8_id).zfill(6) + '.fit')
	# print(path_2_eRO_catalog)
	sf = (sf1) & (HEALPIX_8 == HEALPIX_8_id)
	if len(sf.nonzero()[0]) > 0:
		t = Table()
		##
		# Coordinates
		##
		for col_name, unit_val in zip(f2.columns.names, f2.columns.units):
			t.add_column(Column(name=col_name, data=f2[col_name][sf], unit=unit_val,dtype=n.float32 ) )
		##
		# Galaxy
		##
		for col_name, unit_val in zip(f1.columns.names, f1.columns.units):
			t.add_column(Column(name='galaxy_'+col_name, data=f1[col_name][sf], unit=unit_val, dtype=n.float32 ) )
		##
		# HALO
		##
		for col_name, unit_val in zip(f0.columns.names, f0.columns.units):
			if col_name == 'Mvir' or col_name == 'M200c' or col_name == 'M500c':
				t.add_column(Column(name='HALO_'+col_name, data=f0[col_name][sf]/h, unit=unit_val, dtype=n.float32 ) )
			elif col_name == 'id' or col_name == 'pid' :
				t.add_column(Column(name='HALO_'+col_name, data=f0[col_name][sf], unit=unit_val, dtype=n.int64 ) )
			else:
				t.add_column(Column(name='HALO_'+col_name, data=f0[col_name][sf], unit=unit_val, dtype=n.float32 ) )
		##
		# AGN
		##
		for col_name, unit_val in zip(f3.columns.names, f3.columns.units):
			t.add_column(Column(name=col_name, data=f3[col_name][sf], unit=unit_val,dtype=n.float32 ) )

		# writes
		t.write(path_2_eRO_catalog, overwrite=True)

