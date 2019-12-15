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


f1 = fits.open(path_2_light_cone)
Mvir = f1[1].data['Mvir'] / h
M500c = f1[1].data['M500c'] / h
halo_id = f1[1].data['id']
#if os.path.basename(path_2_light_cone)[:3] == 'sat':
halo_host_id = f1[1].data['pid']
#else:
#	halo_host_id = -n.ones_like(halo_id).astype('int')
f1.close()
print('halo file opened', time.time() - t0)

f1 = fits.open(path_2_galaxy_file)
galaxy_stellar_mass = f1[1].data['SMHMR_mass']
galaxy_star_formation_rate = f1[1].data['star_formation_rate']
galaxy_LX_hard = f1[1].data['LX_hard']
galaxy_mag_r = f1[1].data['mag_r']
f1.close()
print('galaxy file opened', time.time() - t0)

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra']
dec = f2[1].data['dec']
zzr = f2[1].data['redshift_R']
zzs = f2[1].data['redshift_S']
dL_cm = f2[1].data['dL']
galactic_NH = f2[1].data['nH']
galactic_ebv = f2[1].data['ebv']
g_lat = f2[1].data['g_lat']
g_lon = f2[1].data['g_lon']
ecl_lat = f2[1].data['ecl_lat']
ecl_lon = f2[1].data['ecl_lon']
N_galaxies = len(zzr)
f2.close()
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
		t.add_column(Column(name="ra", format='D', unit='deg', data=ra[sf])) 
		t.add_column(Column(name="dec", format='D', unit='deg', data=dec[sf])) 
		t.add_column(Column(name="g_lat", format='D', unit='deg', data=g_lat[sf])) 
		t.add_column(Column(name="g_lon", format='D', unit='deg', data=g_lon[sf])) 
		t.add_column(Column(name="ecl_lat", format='D', unit='deg', data=ecl_lat[sf])) 
		t.add_column(Column(name="ecl_lon", format='D', unit='deg', data=ecl_lon[sf]))   
		##
		# distances
		# extinction maps
		# Galaxy properties
		##
		t.add_column(Column(name="redshift_R", format='D', unit='real space', data=zzr[sf])) 
		t.add_column(Column(name="redshift_S", format='D', unit='redshift space', data=zzs[sf])) 
		t.add_column(Column(name="dL_cm", format='D', unit='cm', data=dL_cm[sf]))
		t.add_column(Column(name="galactic_NH", format='D', unit='cm-2', data=galactic_NH[sf]))
		t.add_column(Column(name="galactic_ebv", format='D', unit='mag', data=galactic_ebv[sf]))
		t.add_column(Column(name="galaxy_stellar_mass", format='D', unit='log10(M/[M_sun])', data=galaxy_stellar_mass[sf])) 
		t.add_column(Column(name="galaxy_star_formation_rate", format='D', unit='log10(SFR/[M_sun/year])', data=galaxy_star_formation_rate[sf]))
		t.add_column(Column(name="galaxy_LX_hard", format='D', unit='log10(LX (2-10keV)/[erg/s])', data=galaxy_LX_hard[sf])) 
		t.add_column(Column(name="galaxy_mag_r", format='D', unit='mag', data=galaxy_mag_r[sf])) 
		##
		# Dark matter halo
		##                                                                                           
		t.add_column(Column(name="HALO_M500c", format='D', unit='log10(M/[M_sun])', data=n.log10(M500c[sf]) )) 
		t.add_column(Column(name="HALO_Mvir", format='D', unit='log10(M/[M_sun])', data=n.log10(Mvir[sf]))) 
		t.add_column(Column(name="HALO_halo_id", format='K', unit='', data=halo_id[sf]))
		t.add_column(Column(name="HALO_pid", format='K', unit='', data=halo_host_id[sf]))
		# writes
		t.write(path_2_eRO_catalog, overwrite=True)
