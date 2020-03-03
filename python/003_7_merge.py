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


env = 'MD10'  # sys.argv[1]

# simulation setup
if env == "MD10" or env == "MD04":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

root_dir = os.path.join(os.environ[env])

path_2_deep   = os.path.join(root_dir, 'AGN_DEEP_4MOST.fits')
path_2_ir   = os.path.join(root_dir, 'AGN_IR_4MOST.fits')
path_2_wide = os.path.join(root_dir, 'AGN_WIDE_4MOST.fits')

#sub_survey_names = n.array([ 'cluster_BCG', 'cluster_redGAL', 'filament_GAL'])

t1 = Table.read(path_2_deep)   
t2 = Table.read(path_2_wide)
t3 = Table.read(path_2_ir)     

area_ero      = lambda t : ( abs(t['g_lat']) > 10 ) & ( t['g_lon'] > 180 ) & ( t['DEC'] < 5 )
area_ero_deep = lambda t : ( abs(t['ecl_lat']) > 80 ) 
magl_lim       = lambda t, mag_lim : (t['MAG'] < mag_lim) 

nl = lambda sel : len(sel.nonzero()[0])
s1 = magl_lim(t1, 23.4) & area_ero( t1 ) & area_ero_deep( t1 ) 
print(nl(s1))
s2 = magl_lim(t2, 22.8) & area_ero( t2 )
print(nl(s2))
s3 = magl_lim(t3, 22.8) & area_ero( t3 )
print(nl(s3))

path_2_output = os.path.join(root_dir, 'S6_4MOST_R-228-234.fit')

t1_b = Table(t1[s1])
t2_b = Table(t2[s2])
t3_b = Table(t3[s3])

#t1_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
#t2_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
#t3_b.remove_columns(['Z','Mstar','SFR','EBV','K_mag_abs','rtot','rfib','ug','gr','ri','iz','zy','yj','jh','hks','g_lat','g_lon'])

test = vstack((t1_b, t2_b, t3_b))

test['MAG_TYPE'][:]="SDSS_r_AB"

test.write (path_2_output , overwrite=True)


