"""
Merges into a single fits catalog containing the 4FS input columns.


"""
import glob
import sys
from astropy_healpix import healpy
import os
from scipy.special import gammainc  # , gamma,  gammaincinv, gammaincc
from scipy.stats import scoreatpercentile
import pandas as pd  # external package
from scipy.special import erf
from astropy.coordinates import SkyCoord
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column, vstack
import astropy.units as u
import numpy as n
import extinction
print('CREATES 4FS FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

nl = lambda sel : len(sel.nonzero()[0])
path_2_output = os.path.join('/data42s/comparat/firefly/mocks/2020-03/QMOST/', 'S6_4MOST_R-228-234.fit')

def get_t123(env):
	print(env)
	path_2_deep = os.path.join(os.environ[env], 'AGN_DEEP_4MOST.fits')
	path_2_wide = os.path.join(os.environ[env], 'AGN_WIDE_4MOST.fits')
	path_2_ir   = os.path.join(os.environ[env], 'AGN_IR_4MOST.fits')
	t1 = Table.read(path_2_deep)   
	t2 = Table.read(path_2_wide)
	t3 = Table.read(path_2_ir)     
	if env=='MD10':
		area_ero      = lambda t : ( abs(t['g_lat']) > 10 ) & ( t['g_lon'] > 180 ) & ( t['DEC'] < 5 ) & (t['redshift_R']>0.3)
		area_ero_deep = lambda t : ( abs(t['ecl_lat']) > 80 ) & (t['redshift_R']>0.3)
		magl_lim       = lambda t, mag_lim : (t['MAG'] < mag_lim) & (t['redshift_R']>0.3)
	if env=='MD04':
		area_ero      = lambda t : ( abs(t['g_lat']) > 10 ) & ( t['g_lon'] > 180 ) & ( t['DEC'] < 5 ) & (t['redshift_R']<=0.3)
		area_ero_deep = lambda t : ( abs(t['ecl_lat']) > 80 ) & (t['redshift_R']<=0.3)
		magl_lim       = lambda t, mag_lim : (t['MAG'] < mag_lim) & (t['redshift_R']<=0.3)
		
	s1 = magl_lim(t1, 23.4) & area_ero( t1 ) & area_ero_deep( t1 ) 
	print('t1 DEEP', nl(s1))
	s2 = magl_lim(t2, 22.8) & area_ero( t2 )
	print('t2 WIDE', nl(s2))
	s3 = magl_lim(t3, 22.8) & area_ero( t3 )
	print('t3 IR', nl(s3))

	t1_b = Table(t1[s1])
	t2_b = Table(t2[s2])
	t3_b = Table(t3[s3])
	#t1_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
	#t2_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
	#t3_b.remove_columns(['Z','Mstar','SFR','EBV','K_mag_abs','rtot','rfib','ug','gr','ri','iz','zy','yj','jh','hks','g_lat','g_lon'])
	test = vstack((t1_b, t2_b, t3_b))
	test['MAG_TYPE'][:]="SDSS_r_AB"
	return test

t_md04 = get_t123('MD04')
t_md10 = get_t123('MD10')
table_S6 = vstack((t_md04, t_md10))

table_S6.write (path_2_output , overwrite=True)


