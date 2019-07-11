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
print('CREATES FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')
import time 
t0 = time.time()
import numpy as n
import os, sys
#import astropy.io.fits as fits
import h5py
import astropy.io.fits as fits
from astropy_healpix import healpy 

# import all pathes 

env = sys.argv[1] # 'MD04'
baseName = sys.argv[2] # "sat_0.62840"
print(env, baseName)

test_dir = os.path.join( os.environ[env], 'hlists', 'fits' )

path_2_light_cone = os.path.join(test_dir, baseName+'.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName+'_coordinates.h5')
path_2_galaxy_file = os.path.join(test_dir, baseName+'_galaxy.h5')
path_2_agn_file = os.path.join(test_dir, baseName+'_agn.h5')

dir_2_eRO_catalog = os.path.join( test_dir, 'cat_AGN_'+baseName )
if os.path.isdir(dir_2_eRO_catalog)==False:
	os.system('mkdir -p '+dir_2_eRO_catalog)


# eRosita flux limits 
path_2_flux_limits = os.path.join( os.environ['GIT_AGN_MOCK'], "data", "erosita", "flux_limits.fits")


# simulation setup 
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
if env=="MD10" or env=="MD04":
	cosmoMD = FlatLambdaCDM(H0=67.77*u.km/u.s/u.Mpc, Om0=0.307115)#, Ob0=0.048206)
	h = 0.6777
	L_box = 1000.0 / h
	cosmo = cosmoMD
if env=="UNIT_fA1_DIR" or env=="UNIT_fA1i_DIR" or env=="UNIT_fA2_DIR":
	cosmoUNIT = FlatLambdaCDM(H0=67.74*u.km/u.s/u.Mpc, Om0=0.308900)
	h = 0.6774
	L_box = 1000.0 / h
	cosmo = cosmoUNIT

f3 = h5py.File(path_2_agn_file, 'r')
FX_hard              = f3["AGN/FX_hard"].value
FX_soft              = f3["AGN/FX_soft"].value
FX_soft_attenuated   = f3["AGN/FX_soft_attenuated"].value
LX_hard              = f3["AGN/LX_hard"].value
LX_soft              = f3["AGN/LX_soft"].value
SDSS_r_AB            = f3["AGN/SDSS_r_AB"].value
SDSS_r_AB_attenuated = f3["AGN/SDSS_r_AB_attenuated"].value
agn_type             = f3["AGN/agn_type"].value
ids_active           = f3["AGN/ids_active"].value
logNH                = f3["AGN/logNH"].value
random               = f3["AGN/random"].value
f3.close()
print('agn file opened', time.time()-t0)

f1 = fits.open(path_2_light_cone)
Mvir = f1[1].data['Mvir'][ids_active]/h 
M500c = f1[1].data['M500c'][ids_active]/h 
halo_id = f1[1].data['id'][ids_active]
if os.path.basename(path_2_light_cone)[:3]=='sat':
	halo_host_id = f1[1].data['pid'][ids_active]
else:
	halo_host_id = -n.ones_like(halo_id).astype('int')
f1.close()
print('halo file opened', time.time()-t0)

f1 = h5py.File(path_2_galaxy_file, 'r')
galaxy_stellar_mass = f1['galaxy/SMHMR_mass'].value[ids_active]
galaxy_star_formation_rate = f1['galaxy/star_formation_rate'].value[ids_active]
galaxy_LX_hard = f1['galaxy/LX_hard'].value[ids_active]
f1.close()
print('galaxy file opened', time.time()-t0)

f2 = h5py.File(path_2_coordinate_file, 'r')
ra = f2['/coordinates/ra'].value[ids_active]
dec = f2['/coordinates/dec'].value[ids_active]
zzr = f2['/coordinates/redshift_R'].value[ids_active]
zzs = f2['/coordinates/redshift_S'].value[ids_active]
dL_cm = f2['/coordinates/dL'].value[ids_active]
galactic_NH = f2['/coordinates/NH'].value[ids_active]
galactic_ebv = f2['/coordinates/ebv'].value[ids_active]
g_lat = f2['/coordinates/g_lat'].value[ids_active]
g_lon = f2['/coordinates/g_lon'].value[ids_active]
ecl_lat = f2['/coordinates/ecl_lat'].value[ids_active]
ecl_lon = f2['/coordinates/ecl_lon'].value[ids_active]
N_galaxies = len(zzr)
f2.close()
print('coordinate file opened', time.time()-t0)

pix_ids = healpy.ang2pix(512, n.pi/2.-g_lat*n.pi/180., g_lon*n.pi/180., nest=True)
flux_lim_data = fits.open(path_2_flux_limits)# [1].data
#flux_limit_eRASS3, flux_limit_eRASS8, flux_limit_SNR3  
flux_limit = flux_lim_data[1].data['flux_limit_eRASS8']
print('flux limit file opened', time.time()-t0)

# selection function here :
sf1 = (FX_soft_attenuated>0) # ( FX_soft_attenuated > 10**(flux_limit[pix_ids]-2) ) & ( SDSS_r_AB_attenuated < 26.5 )
print('selection function applied', time.time()-t0)

HEALPIX_32 = healpy.ang2pix(8,  dec*n.pi/180. + n.pi/2., ra*n.pi/180. ) 

for HEALPIX_32_id in n.arange(healpy.nside2npix(8)):
	path_2_eRO_catalog = os.path.join(dir_2_eRO_catalog, str(HEALPIX_32_id).zfill(6)+'.fit')
	#print(path_2_eRO_catalog)

	sf=(sf1)&(HEALPIX_32==HEALPIX_32_id)

	if len(sf.nonzero()[0])>0:
		hdu_cols  = fits.ColDefs([
		##
		# Coordinates
		##
		fits.Column(name="ra"              , format='D', unit='degree',         array = ra       [sf]   )
		,fits.Column(name= "dec"           , format='D', unit='degree',         array = dec      [sf]   )
		,fits.Column(name= "g_lat"         , format='D', unit='degree',         array = g_lat    [sf]   )
		,fits.Column(name= "g_lon"         , format='D', unit='degree',         array = g_lon    [sf]   )
		,fits.Column(name= "ecl_lat"       , format='D', unit='degree',         array = ecl_lat  [sf]   )
		,fits.Column(name= "ecl_lon"       , format='D', unit='degree',         array = ecl_lon  [sf]   )
		# distances
		,fits.Column(name= "redshift_R"    , format='D', unit='real space',     array = zzr      [sf]   )
		,fits.Column(name= "redshift_S"    , format='D', unit='redshift space', array = zzs      [sf]   )
		,fits.Column(name= "dL_cm"         , format='D', unit='cm',             array = dL_cm    [sf]   )
		# extinction maps
		,fits.Column(name= "galactic_NH"   , format='D', unit='cm-2',        array = galactic_NH   [sf]   )
		,fits.Column(name= "galactic_ebv"  , format='D', unit='mag',       array = galactic_ebv   [sf]   )
		##
		# Galaxy properties
		##
		,fits.Column(name= "galaxy_stellar_mass"  , format='D', unit='log10(M/[M_sun])', array = galaxy_stellar_mass [sf] )
		,fits.Column(name= "galaxy_star_formation_rate"  , format='D', unit='log10(SFR/[M_sun/year])', array = galaxy_star_formation_rate [sf] )
		,fits.Column(name= "galaxy_LX_hard"  , format='D', unit='log10(LX (2-10keV)/[erg/s])', array = galaxy_LX_hard [sf] )
		##
		# Dark matter halo 
		##
		,fits.Column(name=  "HALO_M500c"        , format='D', unit='log10(M/[M_sun])', array = n.log10(M500c[sf])              )
		,fits.Column(name=  "HALO_Mvir"         , format='D', unit='log10(M/[M_sun])', array = n.log10(Mvir[sf])               )
		,fits.Column(name=  "HALO_halo_id"      , format='K', unit='', array = halo_id[sf]               )
		,fits.Column(name=  "HALO_pid"          , format='K', unit='', array = halo_host_id[sf] )
		##
		# AGN properties
		##
		,fits.Column(name= "AGN_LX_soft"          , format='D', unit='log10(Luminosity/[erg/s] 0.5-2 keV)', array = LX_soft[sf]  )
		,fits.Column(name= "AGN_FX_soft"          , format='D', unit='Flux/[erg/cm2/s] 0.5-2 keV',   array = FX_soft_attenuated[sf]    )
		,fits.Column(name= "AGN_LX_hard"          , format='D', unit='log10(Luminosity/[erg/s] 2-10 keV)',  array = LX_hard[sf]   )
		,fits.Column(name= "AGN_FX_hard"          , format='D', unit='Flux/[erg/cm2/s] 2-10 keV',    array = FX_hard[sf]     )
		,fits.Column(name= "AGN_SDSS_r_magnitude" , format='D', unit='mag',                          array = SDSS_r_AB_attenuated[sf]          )
		,fits.Column(name= "AGN_nH"               , format='D', unit='log10(nH/[cm-2])',             array = logNH[sf]               )
		,fits.Column(name= "AGN_random_number"    , format='D', unit='',                             array = random[sf]       )
		,fits.Column(name= "AGN_type"             , format='D', unit='X/opt type: 11, 12, 21, 22',   array = agn_type[sf]           )

		])

		#print('fits columnes created', time.time()-t0)

		tb_hdu = fits.BinTableHDU.from_columns( hdu_cols )
		prihdr = fits.Header()
		prihdr['author'] = 'JC'
		prihdu = fits.PrimaryHDU(header=prihdr)
		thdulist = fits.HDUList([prihdu, tb_hdu])

		if os.path.isfile(path_2_eRO_catalog):
			os.system("rm "+path_2_eRO_catalog)
		thdulist.writeto(path_2_eRO_catalog)
		print('written', path_2_eRO_catalog, time.time()-t0)






