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
from astropy.table import Table, Column
#import h5py
import numpy as n
print('CREATES FITS FILES')
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

f3 = fits.open(path_2_agn_file)
FX_hard = f3[1].data["FX_hard"]
FX_soft = f3[1].data["FX_soft"]
FX_soft_attenuated = f3[1].data["FX_soft_attenuated"]
LX_hard = f3[1].data["LX_hard"]
LX_soft = f3[1].data["LX_soft"]
SDSS_r_AB = f3[1].data["SDSS_r_AB"]
SDSS_r_AB_attenuated = f3[1].data["SDSS_r_AB_attenuated"]
agn_type = f3[1].data["agn_type"]
ids_active = f3[1].data["ids_active"]
logNH = f3[1].data["logNH"]
random = f3[1].data["random"]
f3.close()
print('agn file opened', time.time() - t0)

f1 = fits.open(path_2_light_cone)
Mvir = f1[1].data['Mvir'][ids_active] / h
M500c = f1[1].data['M500c'][ids_active] / h
halo_id = f1[1].data['id'][ids_active]
halo_host_id = f1[1].data['pid'][ids_active]
f1.close()
print('halo file opened', time.time() - t0)

f1 = fits.open(path_2_galaxy_file)
galaxy_stellar_mass = f1[1].data['SMHMR_mass'][ids_active]
galaxy_star_formation_rate = f1[1].data['star_formation_rate'][ids_active]
galaxy_LX_hard = f1[1].data['LX_hard'][ids_active]
galaxy_mag_r = f1[1].data['mag_r'][ids_active]
f1.close()
print('galaxy file opened', time.time() - t0)

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra'][ids_active]
dec = f2[1].data['dec'][ids_active]
zzr = f2[1].data['redshift_R'][ids_active]
zzs = f2[1].data['redshift_S'][ids_active]
dL_cm = f2[1].data['dL'][ids_active]
galactic_NH = f2[1].data['nH'][ids_active]
galactic_ebv = f2[1].data['ebv'][ids_active]
g_lat = f2[1].data['g_lat'][ids_active]
g_lon = f2[1].data['g_lon'][ids_active]
ecl_lat = f2[1].data['ecl_lat'][ids_active]
ecl_lon = f2[1].data['ecl_lon'][ids_active]
N_galaxies = len(zzr)
f2.close()
print('coordinate file opened', time.time() - t0)

pix_ids = healpy.ang2pix(
    512,
    n.pi /
    2. -
    g_lat *
    n.pi /
    180.,
    g_lon *
    n.pi /
    180.,
    nest=True)
flux_lim_data = fits.open(path_2_flux_limits)  # [1].data
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3
flux_limit = flux_lim_data[1].data['flux_limit_eRASS8']
print('flux limit file opened', time.time() - t0)

# selection function here :
# ( FX_soft_attenuated > 10**(flux_limit[pix_ids]-2) ) & ( SDSS_r_AB_attenuated < 26.5 )
sf1 = (FX_soft_attenuated > 0)
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
		t.add_column(Column(name="ra", format='D', unit='degree', data=ra[sf])) 
		t.add_column(Column(name="dec", format='D', unit='degree', data=dec[sf])) 
		t.add_column(Column(name="g_lat", format='D', unit='degree', data=g_lat[sf])) 
		t.add_column(Column(name="g_lon", format='D', unit='degree', data=g_lon[sf])) 
		t.add_column(Column(name="ecl_lat", format='D', unit='degree', data=ecl_lat[sf])) 
		t.add_column(Column(name="ecl_lon", format='D', unit='degree', data=ecl_lon[sf]))            # distances
		# extinction maps
		# Galaxy properties
		t.add_column(Column(name="redshift_R", format='D', unit='real space', data=zzr[sf])) 
		t.add_column(Column(name="redshift_S", format='D', unit='redshift space', data=zzs[sf])) 
		t.add_column(Column(name="dL_cm", format='D', unit='cm', data=dL_cm[sf])) 
		t.add_column(Column(name="galactic_NH", format='D', unit='cm-2', data=galactic_NH[sf])) 
		t.add_column(Column(name="galactic_ebv", format='D', unit='mag', data=galactic_ebv[sf])) 
		t.add_column(Column(name="galaxy_stellar_mass", format='D', unit='log10(M/[M_sun])', data=galaxy_stellar_mass[sf])) 
		t.add_column(Column(name="galaxy_star_formation_rate", format='D', unit='log10(SFR/[M_sun/year])', data=galaxy_star_formation_rate[sf])) 
		t.add_column(Column(name="galaxy_LX_hard", format='D', unit='log10(LX (2-10keV)/[erg/s])', data=galaxy_LX_hard[sf]))            # Dark matter halo
		t.add_column(Column(name="galaxy_mag_r", format='D', unit='mag', data=galaxy_mag_r[sf])) 
		##
		# Dark matter halo
		##                                                                                           
		t.add_column(Column(name="HALO_M500c", format='D', unit='log10(M/[M_sun])', data=n.log10(M500c[sf])))
		t.add_column(Column(name="HALO_Mvir", format='D', unit='log10(M/[M_sun])', data=n.log10(Mvir[sf])))
		t.add_column(Column(name="HALO_halo_id", format='K', unit='', data=halo_id[sf]))
		t.add_column(Column(name="HALO_pid", format='K', unit='', data=halo_host_id[sf])) 
		##
		# AGN properties
		##
		t.add_column(Column(name="AGN_LX_soft", format='D', unit='log10(Luminosity/[erg/s] 0.5-2 keV)', data=LX_soft[sf]))
		t.add_column(Column(name="AGN_FX_soft", format='D', unit='Flux/[erg/cm2/s] 0.5-2 keV', data=FX_soft_attenuated[sf])) 
		t.add_column(Column(name="AGN_LX_hard", format='D', unit='log10(Luminosity/[erg/s] 2-10 keV)', data=LX_hard[sf])) 
		t.add_column(Column(name="AGN_FX_hard", format='D', unit='Flux/[erg/cm2/s] 2-10 keV', data=FX_hard[sf])) 
		t.add_column(Column(name="AGN_SDSS_r_magnitude", format='D', unit='mag', data=SDSS_r_AB_attenuated[sf])) 
		t.add_column(Column(name="AGN_nH", format='D', unit='log10(nH/[cm-2])', data=logNH[sf]))
		t.add_column(Column(name="AGN_random_number", format='D', unit='', data=random[sf])) 
		t.add_column(Column(name="AGN_type", format='D', unit='X/opt type: 11, 12, 21, 22', data=agn_type[sf]))
		#print('fits columnes created', time.time()-t0)
		t.write(path_2_eRO_catalog, overwrite=True)
		print('written', path_2_eRO_catalog, time.time() - t0)
