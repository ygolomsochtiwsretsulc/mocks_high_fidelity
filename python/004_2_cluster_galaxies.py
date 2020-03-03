"""
What it does
------------

Take the complete cluster file ($ENV_eRO_CLU.fit) and identifies satellites halos around clusters (within 1 angular rvir and 1 rvir in 3D) in each galaxy shell. Creates a set of satellite files.
Uses BallTrees to find neighbours efficiently.

Create a satellite file for each shell/ftype.
Then concatenates (stilts/topcat) into a (very big) single file. $ENV_eRO_CLU_SAT.fit


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
from astropy.table import Table, Column
	
print('CREATES FITS FILE with galaxies around clusters')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

env = sys.argv[1]
baseName = sys.argv[2]
delta_crit = sys.argv[3]
#env="MD10" 
#baseName="all_0.89510"
#delta_crit = '200c'
#delta_crit = '500c'
#delta_crit = 'vir'
#delta_crit = '2rvir'

z_snap = 1./float(baseName.split('_')[1])-1.
aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
print(env, baseName, delta_crit)

test_dir = os.path.join(os.environ[env], 'fits')

# input cluster catalog
path_2_CLU_catalog = os.path.join(os.environ[env], env + '_eRO_CLU.fit')

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

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_out_file = os.path.join(test_dir, baseName + '_galaxiesAroundClusters.fit')

f1 = Table.read(path_2_galaxy_file)

f2 = Table.read(path_2_coordinate_file)
zzr = f2['redshift_R']

shell_zmin = n.min(zzr)
shell_zmax = n.max(zzr)

# input columns of interest
hdu_clu_i = Table.read(path_2_CLU_catalog)
# selecte redshift range around the shell
s_z = (hdu_clu_i['redshift_R']>shell_zmin-0.05) & (hdu_clu_i['redshift_R']<shell_zmax+0.05)
# define sub set
hdu_clu = hdu_clu_i[s_z]
rvir_CLU = hdu_clu['HALO_Rvir']
x_coord_CLU = hdu_clu['HALO_x']
y_coord_CLU = hdu_clu['HALO_y']
z_coord_CLU = hdu_clu['HALO_z']


omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)

HOST_HALO_Mvir = hdu_clu['HALO_Mvir'] / h
HOST_HALO_Rvir = hdu_clu['HALO_Rvir']
HOST_HALO_M500c = hdu_clu['HALO_M500c'] / h
HOST_HALO_R500c = (DeltaVir_bn98(z_snap)/500. * HOST_HALO_M500c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir
HOST_HALO_M200c = hdu_clu['HALO_M200c'] / h
HOST_HALO_R200c = (DeltaVir_bn98(z_snap)/200. * HOST_HALO_M200c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir

if delta_crit == '200c' :
	frac_rvir = HOST_HALO_R200c/HOST_HALO_Rvir
	RADIUS = HOST_HALO_R200c
if delta_crit == '500c' :
	frac_rvir = HOST_HALO_R500c/HOST_HALO_Rvir
	RADIUS = HOST_HALO_R500c
if delta_crit == 'vir' :
	frac_rvir = HOST_HALO_Rvir/HOST_HALO_Rvir
	RADIUS = hdu_clu_bin['HOST_HALO_Rvir']
if delta_crit == '2rvir' :
	frac_rvir = 2.*n.ones_like(HOST_HALO_R200c)
	RADIUS = 2.*hdu_clu_bin['HOST_HALO_Rvir']


f3 = Table.read(path_2_light_cone)
halo_id = f3['id']
halo_pid = f3['pid']
x = f3['x']
y = f3['y']
z = f3['z']

# 3D search
Tree_obj = BallTree(n.transpose([x, y, z]))
idx_arrays, distances = Tree_obj.query_radius(n.transpose([x_coord_CLU, y_coord_CLU, z_coord_CLU]), r = RADIUS/1000. , return_distance=True)
N_in_test = n.array([len(t_i) for t_i in idx_arrays])

N_in_test = n.array([len(t_i) for t_i in idx_arrays])
# cluster ID: line in the cluster file for each satellite
# for hdu_clu
CLU_ids = n.hstack((n.array([n.ones(el) * ii for ii, el in enumerate(N_in_test)]))).astype('int')
# IDs of the satellite in the galaxy and coordinate file
# for f1, f2, f3
sf = n.hstack((idx_arrays))

# distance to cluster_center
#r_o_rvir = n.hstack((distances)) / (rvir_CLU[CLU_ids] / 1000.)


t = Table()

for col_name in f1.columns.keys():
	print(col_name,)
	if col_name=='LX_hard':
		print('--not added',)
	else:
		t.add_column(Column(name=col_name, data=f1[col_name][sf] ) )

for col_name in f2.columns.keys():
	print(col_name,)
	if col_name=='LX_hard':
		print('--not added',)
	else:
		t.add_column(Column(name=col_name, data=f2[col_name][sf] ) )

for col_name in f3.columns.keys():
	print(col_name,)
	if col_name in n.array(['M200c', 'M500c', 'Xoff', 'b_to_a_500c', 'c_to_a_500c', 'scale_of_last_MM', 'Acc_Rate_1_Tdyn']):
		print('--not added',)
	else:
		t.add_column(Column(name=col_name, data=f3[col_name][sf] ) )


for col_name in hdu_clu.columns.keys():
	print(col_name,)
	if col_name in n.array(['g_lon', 'g_lat', 'ecl_lon', 'ecl_lat', 'dL', 'nH', 'ebv', 'HALO_lineID','galaxy_SMHMR_mass', 'HALO_c_to_a_500c',  'galaxy_star_formation_rate', 'galaxy_is_quiescent',  'galaxy_LX_hard', 'galaxy_mag_abs_r',  'galaxy_mag_r']):
		print('--not added',)
	else:
		t.add_column(Column(name='HOST_'+col_name, data=hdu_clu[col_name][CLU_ids] ) )


t.write(path_2_out_file, overwrite=True)
print('written to ',path_2_out_file, 'after ', time.time()-t0, ' seconds')
