"""
What it does
------------

Creates a light cone with every halo with mass > XX

"""
from astropy.table import Table, Column
from astropy_healpix import healpy
import sys
import os
import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import h5py
import numpy as n
from lib_cluster import write_img, create_matrix
print('Creates the light cone files per shell')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


# import all pathes

env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "all_0.62840"
z_snap = 1./float(baseName.split('_')[1])-1.
print(env, baseName,z_snap)
make_figure = True
make_figure = False

# initializes pathes to files
test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_CLU_file = os.path.join(test_dir, baseName + '_LC.fits')

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

omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)


f1 = fits.open(path_2_light_cone)
N_obj = len(f1[1].data['M500c'])
cluster = ( n.log10(f1[1].data['M500c'] / h) > 13 ) & ( f1[1].data['pid'] == -1 )

Mvir = f1[1].data['Mvir'][cluster] / h
Rvir = f1[1].data['Rvir'][cluster]
M500c = f1[1].data['M500c'][cluster] / h
R500c = (DeltaVir_bn98(z_snap)/500. * M500c / Mvir)**(1./3.)*Rvir
logM500c = n.log10(M500c)
log_vmax = n.log10(f1[1].data['vmax'][cluster])
scale_of_last_MM = f1[1].data['scale_of_last_MM'][cluster]
N_clu = len(Mvir)
print(N_clu, 'clusters')
Xoff = f1[1].data['Xoff'][cluster] / f1[1].data['Rvir'][cluster]
b_to_a_500c = f1[1].data['b_to_a_500c'][cluster]
#f1.close()
halo_lineID = n.arange(N_obj)[cluster]
ids_cluster = (n.arange(N_obj)[cluster]*1e12 + f1[1].data['id'][cluster] ).astype('int')
ids_cluster_str = n.array([str(el).zfill(22) for el in ids_cluster])

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra'][cluster]
dec = f2[1].data['dec'][cluster]
zz = f2[1].data['redshift_R'][cluster]
zzs = f2[1].data['redshift_S'][cluster]
dL_cm = f2[1].data['dL'][cluster]
galactic_NH = f2[1].data['nH'][cluster]
galactic_ebv = f2[1].data['ebv'][cluster]
g_lat = f2[1].data['g_lat'][cluster]
g_lon = f2[1].data['g_lon'][cluster]
ecl_lat = f2[1].data['ecl_lat'][cluster]
ecl_lon = f2[1].data['ecl_lon'][cluster]
#f2.close()
print('coordinate file opened', time.time() - t0)

# cosmological volume
zmin = n.min(zz)
zmax = n.max(zz)
z_mean = 0.5 * (zmin + zmax)
print(zmin, '<z<', zmax)
vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value)
DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
print('volume', vol, 'Mpc3', time.time() - t0)

f3 = fits.open(path_2_galaxy_file)
mass = f3[1].data['SMHMR_mass'][cluster]  # log of the stellar mass
# f3.close()

t = Table()

for col_name, unit_val in zip(f2[1].data.columns.names, f2[1].data.columns.units):
	t.add_column(Column(name=col_name, data=f2[1].data[col_name][cluster], unit=unit_val,dtype=n.float32 ) )

for col_name in f1[1].data.columns.names:
	if col_name == 'Mvir' or col_name == 'M200c' or col_name == 'M500c':
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster]/h, unit='', dtype=n.float32 ) )
	else:
		t.add_column(Column(name='HALO_'+col_name, data=f1[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.add_column(Column(name='HALO_lineID', data=halo_lineID.astype('int'), unit='', dtype=n.int64 ) )

for col_name, unit_val in zip(f3[1].data.columns.names, f3[1].data.columns.units):
	t.add_column(Column(name='galaxy_'+col_name, data=f3[1].data[col_name][cluster], unit='', dtype=n.float32 ) )

t.write(path_2_CLU_file, overwrite=True)
