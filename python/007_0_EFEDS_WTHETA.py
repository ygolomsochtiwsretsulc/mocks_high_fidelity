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
path_2_out_file = os.path.join(test_dir, baseName + '_EFEDSwtheta.fits')

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

f2 = fits.open(path_2_coordinate_file)
ra = f2[1].data['ra']
dec = f2[1].data['dec']
zz = f2[1].data['redshift_R']
cluster = ( zz < 6 ) & (abs(dec)<10) & (ra > 130) & (ra < 210 )
print(len(ra[cluster]), 'points')

t = Table()

f1 = fits.open(path_2_light_cone)
t.add_column(Column(name='pid' , data=f1[1].data['pid'][cluster] ) )
t.add_column(Column(name='Mvir' , data=f1[1].data['Mvir'][cluster] / h ) )
t.add_column(Column(name='Rvir' , data=f1[1].data['Rvir'][cluster] ) )
t.add_column(Column(name='vmax' , data=f1[1].data['vmax'][cluster] ) )

t.add_column(Column(name='ra' , data= f2[1].data['ra'][cluster]                  ) )
t.add_column(Column(name='dec' , data=  f2[1].data['dec'][cluster]               ) )
t.add_column(Column(name='redshift_R' , data = f2[1].data['redshift_R'][cluster] ) )
t.add_column(Column(name='redshift_S' , data = f2[1].data['redshift_S'][cluster] ) )
t.add_column(Column(name='dL' , data=  f2[1].data['dL'][cluster]                 ) )
t.add_column(Column(name='nH' , data=  f2[1].data['nH'][cluster]                 ) )
t.add_column(Column(name='ebv' , data=  f2[1].data['ebv'][cluster]               ) )

f3 = fits.open(path_2_galaxy_file)
t.add_column(Column(name='SMHMR_mass' , data= f3[1].data['SMHMR_mass'][cluster] ) ) 
f1.close()
f1.close()
f2.close()
f3.close()

t.write(path_2_out_file, overwrite=True)
