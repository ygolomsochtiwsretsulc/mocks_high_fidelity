"""
What it does
------------

Take the complete cluster file ($ENV_eRO_CLU.fit) and identifies satellites halos around clusters (within 1 angular rvir and 1 rvir in 3D) in each galaxy shell. Creates a set of satellite files.
Uses BallTrees to find neighbours efficiently.

Create a satellite file for each shell/ftype.
Then concatenates (stilts/topcat) into a (very big) single file. $ENV_eRO_CLU_SAT.fit

References
----------

Command to run
--------------

python3 004_2_cluster_galaxies.py environmentVAR laptop

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/


laptop: if running on the server or on the laptop. Redirects properly to the stilts command

Dependencies
------------

topcat/stilts
import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib, sklearn


"""

import glob
import sys
from astropy_healpix import healpy
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from sklearn.neighbors import BallTree
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import h5py
import numpy as n
print('CREATES FITS FILE with galaxies around clusters')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits
# import all pathes
deg_to_rad = n.pi / 180.
#from astropy.table import Table,unique
#from math import radians, cos, sin, asin, sqrt, pi
#from sklearn.neighbors import DistanceMetric

env = sys.argv[1]
laptop = sys.argv[2]  # 'True'

if laptop == "True":
    stilts_cmd = 'java -jar /home/comparat/software/stilts.jar'
else:
    stilts_cmd = 'stilts'

test_dir = os.path.join(os.environ[env], 'hlists', 'fits')
tmp_dir = os.path.join(os.environ[env], 'hlists', 'fits', 'tmp')
if os.path.isdir(tmp_dir) == False:
    os.system('mkdir -p ' + tmp_dir)

# input cluster catalog
path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
# output cluster galaxy catalog
path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')


# simulation setup
if env == "MD10" or env == "MD04":
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

# input columns of interest
hdu_clu = fits.open(path_2_CLU_catalog)
id_CLU = hdu_clu[1].data['halo_id']
rvir_CLU = hdu_clu[1].data['HALO_rvir']
x_coord_CLU = hdu_clu[1].data['x']
y_coord_CLU = hdu_clu[1].data['y']
z_coord_CLU = hdu_clu[1].data['z']
ra_CLU = hdu_clu[1].data['ra']
dec_CLU = hdu_clu[1].data['dec']
zr_CLU = hdu_clu[1].data['redshift_R']
zs_CLU = hdu_clu[1].data['redshift_S']
hdu_clu.close()

# interpolation of angular diameter distance with redshift
zzs = n.arange(0, n.max(zr_CLU) + 0.1, 0.01)
itp_angular_size = interp1d(zzs, n.array(
    [cosmo.kpc_comoving_per_arcmin(zz_i).value for zz_i in zzs]))
search_rad_1_rvir_arcmin = rvir_CLU / itp_angular_size(zr_CLU)
search_rad_1_rvir_deg = search_rad_1_rvir_arcmin / 60.

# get all possible galaxy files
path_2_sat_galaxy_files = sorted(
    n.array(
        glob.glob(
            os.path.join(
                test_dir,
                '*_galaxy.h5'))))

baseNames = sorted(n.array([os.path.basename(path_2_sat_galaxy_file)[
                   :-10] for path_2_sat_galaxy_file in path_2_sat_galaxy_files]))

# loop over all these files to create fit catalogs of galaxies around
# clusters in 2 rvir angular


def get_data(
        path_2_light_cone,
        path_2_coordinate_file,
        path_2_galaxy_file,
        N_rvir_angular=1.):
    """
    parameters:
     - path_2_light_cone: light cone file
     - path_2_coordinate_file: coordinate file
     - path_2_galaxy_file: galaxy properties file
     - N_rvir_angular: search within X angular rvir, Default=2.

    output:
     - hdu_cols: hdu data unit containing the matched satellite file.

    Filter the galaxy catalog to every galaxy within 1_rvir angular

    Then retains galaxies within 1 r_vir (3D positions)

    Outputs the catalogue containing ( < 1 rvir angular) AND (< 1 r_vir) AND (r magnitude observed < 26.5) of a cluster.

    """
    f1 = h5py.File(path_2_galaxy_file, 'r')
    galaxy_stellar_mass = f1['galaxy/SMHMR_mass'][:]
    galaxy_star_formation_rate = f1['galaxy/star_formation_rate'][:]
    galaxy_LX_hard = f1['galaxy/LX_hard'][:]
    galaxy_mag_r = f1['galaxy/mag_r'][:]
    galaxy_mag_abs_r = f1['galaxy/mag_abs_r'][:]
    is_quiescent = f1['galaxy/is_quiescent'][:]
    f1.close()

    f2 = h5py.File(path_2_coordinate_file, 'r')
    ra = f2['/coordinates/ra'][:]
    dec = f2['/coordinates/dec'][:]
    zzr = f2['/coordinates/redshift_R'][:]
    zzs = f2['/coordinates/redshift_S'][:]
    dL_cm = f2['/coordinates/dL'][:]
    galactic_NH = f2['/coordinates/NH'][:]
    galactic_ebv = f2['/coordinates/ebv'][:]
    g_lat = f2['/coordinates/g_lat'][:]
    g_lon = f2['/coordinates/g_lon'][:]
    ecl_lat = f2['/coordinates/ecl_lat'][:]
    ecl_lon = f2['/coordinates/ecl_lon'][:]
    N_galaxies = len(zzr)
    f2.close()

    f1 = fits.open(path_2_light_cone)
    halo_id = f1[1].data['id']
    x = f1[1].data['x']
    y = f1[1].data['y']
    z = f1[1].data['z']

    # angular search radius
    rsearch_rad = search_rad_1_rvir_deg * N_rvir_angular * n.pi / 180.

    # coordinates of the galaxies
    obj_coord_rad = deg_to_rad * n.array([dec, ra]).T
    Tree_obj_Haversine = BallTree(obj_coord_rad, metric='haversine')
    # Looks at individual clusters
    # retrieve targets within a circle of N angular rvir
    # returns the IDS of objects within the search radius for each tile
    idx_arrays, distances = Tree_obj_Haversine.query_radius(
        deg_to_rad * n.transpose([dec_CLU, ra_CLU]), r=rsearch_rad, return_distance=True)

    # computes the number of targets in each tile
    N_in_test = n.array([len(t_i) for t_i in idx_arrays])

    CLU_ids = n.hstack(
        (n.array([n.ones(el) * ii for ii, el in enumerate(N_in_test)]))).astype('int')

    sf = n.hstack((idx_arrays))

    angular_distance_to_cluster = n.hstack(
        (distances)) / (search_rad_1_rvir_deg[CLU_ids] * n.pi / 180.)
    r_o_rvir = ((x_coord_CLU[CLU_ids] - x[sf])**2. + (y_coord_CLU[CLU_ids] - y[sf])
                ** 2. + (z_coord_CLU[CLU_ids] - z[sf])**2.)**(0.5) / (rvir_CLU[CLU_ids] / 1000.)
    redshift_distance_R = zzr[sf] - zr_CLU[CLU_ids]
    redshift_distance_S = zzs[sf] - zs_CLU[CLU_ids]

    s_R = (r_o_rvir < 1) & (galaxy_mag_r[sf] < 26.5)

    hdu_cols = fits.ColDefs([
        ##
        # Coordinates
        ##
        fits.Column(name="ra", format='D', unit='degree', array=ra[sf[s_R]]), fits.Column(name="dec", format='D', unit='degree', array=dec[sf[s_R]]), fits.Column(name="g_lat", format='D', unit='degree', array=g_lat[sf[s_R]]), fits.Column(
            name="g_lon", format='D', unit='degree', array=g_lon[sf[s_R]]), fits.Column(name="ecl_lat", format='D', unit='degree', array=ecl_lat[sf[s_R]]), fits.Column(name="ecl_lon", format='D', unit='degree', array=ecl_lon[sf[s_R]])        # distances
        # extinction maps
        # Galaxy properties
        , fits.Column(name="redshift_R", format='D', unit='real space', array=zzr[sf[s_R]]), fits.Column(name="redshift_S", format='D', unit='redshift space', array=zzs[sf[s_R]]), fits.Column(name="dL_cm", format='D', unit='cm', array=dL_cm[sf[s_R]]), fits.Column(name="galactic_NH", format='D', unit='cm', array=galactic_NH[sf[s_R]]), fits.Column(name="galactic_ebv", format='D', unit='cm', array=galactic_ebv[sf[s_R]]), fits.Column(name="galaxy_stellar_mass", format='D', unit='M_sun', array=galaxy_stellar_mass[sf[s_R]]), fits.Column(name="galaxy_star_formation_rate", format='D', unit='M_sun', array=galaxy_star_formation_rate[sf[s_R]]), fits.Column(name="galaxy_LX_hard", format='D', unit='M_sun', array=galaxy_LX_hard[sf[s_R]]), fits.Column(name="galaxy_mag_r", format='D', unit='', array=galaxy_mag_r[sf[s_R]]), fits.Column(name="galaxy_mag_abs_r", format='D', unit='', array=galaxy_mag_abs_r[sf[s_R]]), fits.Column(name="is_quiescent", format='L', unit='', array=is_quiescent[sf[s_R]])        # Dark matter halo
        ##
        # ,fits.Column(name=  "scale_of_last_MM"   , format='D', unit='',     array = f1[1].data["scale_of_last_MM"   ][sf[s_R]] )
        , fits.Column(name="HALO_M200c", format='D', unit='log10(M_sun)', array=f1[1].data["M200c"][sf[s_R]] / h), fits.Column(name="HALO_M500c", format='D', unit='log10(M_sun)', array=f1[1].data["M500c"][sf[s_R]] / h), fits.Column(name="HALO_Mvir", format='D', unit='log10(M_sun)', array=f1[1].data["Mvir"][sf[s_R]] / h), fits.Column(name="HALO_Acc_Rate_1Tdyn", format='D', unit='Msun/yr', array=f1[1].data["Acc_Rate_1_Tdyn"][sf[s_R]]), fits.Column(name="HALO_rs", format='D', unit='kpc', array=f1[1].data["rs"][sf[s_R]]), fits.Column(name="HALO_rvir", format='D', unit='kpc', array=f1[1].data["rvir"][sf[s_R]]), fits.Column(name="HALO_vmax", format='D', unit='km/s', array=f1[1].data["vmax"][sf[s_R]]), fits.Column(name="Vrms", format='D', unit='km/s', array=(f1[1].data["vx"][sf[s_R]]**2 + f1[1].data["vy"][sf[s_R]]**2 + f1[1].data["vz"][sf[s_R]]**2)**0.5)        #
        , fits.Column(name="HALO_id", format='K', unit='', array=halo_id[sf[s_R]]), fits.Column(name="HALO_host_id", format='K', unit='', array=CLU_ids[s_R]), fits.Column(name="angular_distance_to_cluster_center_in_rvir", format='D', unit='theta/theta_vir', array=angular_distance_to_cluster[s_R]), fits.Column(name="comoving_distance_to_cluster_in_rvir", format='D', unit='r/r_vir', array=r_o_rvir[s_R]), fits.Column(name="redshift_R_distance_to_cluster", format='D', unit='z_R-z_R_cluster', array=redshift_distance_R[s_R]), fits.Column(name="redshift_S_distance_to_cluster", format='D', unit='z_S-z_S_cluster', array=redshift_distance_S[s_R])
    ])
    f1.close()
    return hdu_cols


# LOOPS OVER ALL FILES AND CREATE THE TEMPORARY FILE
hdu_col_array = []
for baseName in baseNames[::-1]:
    path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
    path_2_coordinate_file = os.path.join(
        test_dir, baseName + '_coordinates.h5')
    path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.h5')
    path_2_out_file = os.path.join(
        tmp_dir, baseName + '_galaxiesAroundClusters.fit')

    hdu_col = get_data(
        path_2_light_cone,
        path_2_coordinate_file,
        path_2_galaxy_file)
    tb_hdu = fits.BinTableHDU.from_columns(hdu_col)
    prihdr = fits.Header()
    prihdr['author'] = 'JC'
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tb_hdu])
    if os.path.isfile(path_2_out_file):
        os.system("rm " + path_2_out_file)
    thdulist.writeto(path_2_out_file)
    print('written', path_2_out_file, time.time() - t0)


# concatenates all temporary files into a single fits catalog for the clusters
os.chdir(tmp_dir)
c1 = 'ls *_galaxiesAroundClusters.fit > fit_list_galaxiesAroundClusters.list'
print(c1)
os.system(c1)
c2 = stilts_cmd + """ tcat in=@fit_list_galaxiesAroundClusters.list ifmt=fits omode=out ofmt=fits out=""" + path_2_CLU_SAT_catalog
print(c2)
os.system(c2)
os.system('rm -rf ' + tmp_dir)
