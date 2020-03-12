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


env = sys.argv[1]

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

path_2_clusterBCG   = os.path.join(root_dir, 'S5_BCG_4MOST.fit')
path_2_filament   = os.path.join(root_dir, 'FILAMENTS5_4MOST.fits')
path_2_clusterredGAL = os.path.join(root_dir, 'S5_CGAL_4MOST.fit')
path_2_output = os.path.join(root_dir, 'S5_4MOST.fit')

#sub_survey_names = n.array([ 'cluster_BCG', 'cluster_redGAL', 'filament_GAL'])

t1 = Table.read(path_2_clusterBCG)   
t2 = Table.read(path_2_clusterredGAL)
t3 = Table.read(path_2_filament)     

coords = SkyCoord(t3['RA'], t3['DEC'], unit='deg', frame='icrs')
bb_gal = coords.galactic.b.value
ll_gal = coords.galactic.l.value
t3.add_column(Column(name='g_lat', data=bb_gal, unit='deg'))
t3.add_column(Column(name='g_lon', data=ll_gal, unit='deg'))

area_ero = lambda dec, g_lat, g_lon : (abs(g_lat) > 10) & (g_lon > 180) & (dec < 20)
magl_lim = lambda mag, mag_lim : (mag < mag_lim) 
selection = lambda t, mag_lim : area_ero(t['DEC'], t['g_lat'], t['g_lon']) & magl_lim( t['MAG'], mag_lim )
nl = lambda sel : len(sel.nonzero()[0])
#t = t1
s1 = selection(t1, 24.5)
print(nl(s1))
s2 = selection(t2, 22.0)
print(nl(s2))
s3 = selection(t3, 20.5)
print(nl(s3))
path_2_output = os.path.join(root_dir, 'S5_4MOST_RFIB215.fit')

t1_b = Table(t1[s1])
t2_b = Table(t2[s2])
t3_b = Table(t3[s3])

t1_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
t2_b.remove_columns(['ra','dec','g_lat','g_lon','ecl_lat','ecl_lon','redshift_R','redshift_S','dL_cm','galactic_NH','galactic_ebv','galaxy_stellar_mass','galaxy_star_formation_rate','galaxy_LX_hard','is_quiescent','HALO_M200c','HALO_M500c','HALO_Mvir','HALO_Acc_Rate_1Tdyn','HALO_rs','HALO_rvir','HALO_vmax','Vrms','HALO_id','HALO_host_id','angular_distance_to_cluster_center_in_rvir','comoving_distance_to_cluster_in_rvir','redshift_R_distance_to_cluster','redshift_S_distance_to_cluster','galaxy_gr','galaxy_ri','galaxy_iz','richness','galaxy_mag_abs_r','galaxy_mag_r','sdss_g','sdss_r','sdss_i','sdss_z','sdss_g_err','sdss_r_err','sdss_i_err','sdss_z_err'])
t3_b.remove_columns(['Z','Mstar','SFR','EBV','K_mag_abs','rtot','rfib','ug','gr','ri','iz','zy','yj','jh','hks','g_lat','g_lon'])
t3_c=Table(t3_b)
t3_c.remove_columns('NAME')
t3_c.add_column(Column(name='NAME', data=t3_b['NAME']), index=0)

test = vstack((t1_b, t2_b, t3_c))

test['MAG_TYPE'][:]="SDSS_r_AB"

test.write (path_2_output , overwrite=True)


