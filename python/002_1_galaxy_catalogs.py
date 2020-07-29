"""
What it does
------------

Creates a GALAXY fits catalog for each healpix pixel of 13 deg2

Command to run
--------------

python3 002_1_galaxy_catalogs.py environmentVAR fileBasename

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib

"""
import sys, os, time, astropy #, h5py

if sys.version[0]=='2':
    import healpy
if sys.version[0]=='3':
    from astropy_healpix import healpy

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.io.fits as fits
from astropy.table import Table, Column
import numpy as n
print('CREATES GALAXY FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits

# import all pathes

env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "sat_0.62840"
print(env, baseName)

test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')

dir_2_eRO_catalog = os.path.join(test_dir, 'cat_GALAXY_' + baseName)
if os.path.isdir(dir_2_eRO_catalog) == False:
    os.system('mkdir -p ' + dir_2_eRO_catalog)

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


f0 = fits.open(path_2_light_cone)
#Mvir = f0[1].data['Mvir'] / h
#M500c = f0[1].data['M500c'] / h
#halo_id = f0[1].data['id']
#if os.path.basename(path_2_light_cone)[:3] == 'sat':
#halo_host_id = f0[1].data['pid']
#else:
#	halo_host_id = -n.ones_like(halo_id).astype('int')
#f0.close()
print('halo file opened', time.time() - t0)

f1 = fits.open(path_2_galaxy_file)
#galaxy_stellar_mass = f1[1].data['SMHMR_mass']
#galaxy_star_formation_rate = f1[1].data['star_formation_rate']
#galaxy_LX_hard = f1[1].data['LX_hard']
#galaxy_mag_r = f1[1].data['mag_r']
#f1.close()
print('galaxy file opened', time.time() - t0)

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra']
dec = f2[1].data['dec']
#zzr = f2[1].data['redshift_R']
#zzs = f2[1].data['redshift_S']
#dL_cm = f2[1].data['dL']
#galactic_NH = f2[1].data['nH']
#galactic_ebv = f2[1].data['ebv']
#g_lat = f2[1].data['g_lat']
#g_lon = f2[1].data['g_lon']
#ecl_lat = f2[1].data['ecl_lat']
#ecl_lon = f2[1].data['ecl_lon']
#N_galaxies = len(zzr)
#f2.close()
print('coordinate file opened', time.time() - t0)

HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - dec*n.pi/180. , ra*n.pi/180. , nest=True)

for HEALPIX_8_id in n.arange(healpy.nside2npix(8)):
	"""
	Loops over healpix pixels and writes the files to path_2_eRO_catalog
	"""
	path_2_eRO_catalog = os.path.join(dir_2_eRO_catalog, str(HEALPIX_8_id).zfill(6) + '.fit')
	# print(path_2_eRO_catalog)
	sf = (HEALPIX_8 == HEALPIX_8_id)
	if len(sf.nonzero()[0]) > 0:
		t = Table()
		##
		# Coordinates
		##
		for col_name, unit_val in zip(f2[1].data.columns.names, f2[1].data.columns.units):
			t.add_column(Column(name=col_name, data=f2[1].data[col_name][sf], unit=unit_val,dtype=n.float32 ) )
		##
		# Galaxy
		##
		for col_name, unit_val in zip(f1[1].data.columns.names, f1[1].data.columns.units):
			t.add_column(Column(name='galaxy_'+col_name, data=f1[1].data[col_name][sf], unit=unit_val, dtype=n.float32 ) )
		##
		# HALO
		##
		for col_name, unit_val in zip(f0[1].data.columns.names, f0[1].data.columns.units):
			if col_name == 'Mvir' or col_name == 'M200c' or col_name == 'M500c':
				t.add_column(Column(name='HALO_'+col_name, data=f0[1].data[col_name][sf]/h, unit=unit_val, dtype=n.float32 ) )
			elif col_name == 'id' or col_name == 'pid' :
				t.add_column(Column(name='HALO_'+col_name, data=f0[1].data[col_name][sf], unit=unit_val, dtype=n.int64 ) )
			else:
				t.add_column(Column(name='HALO_'+col_name, data=f0[1].data[col_name][sf], unit=unit_val, dtype=n.float32 ) )
		# writes
		t.write(path_2_eRO_catalog, overwrite=True)

